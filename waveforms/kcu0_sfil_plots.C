// kcu0_sfil_plots.C: the KCU0 SFIL channels of run 3062 from data/kcu0_3062.root: raw waveforms, pedestal-subtracted
// waveforms, TOA code and TOT code versus sample, North and South end of every cell of SFIL1-4, with event-quality
// and paddle-trigger flags.
//
//   root -l -b -q 'kcu0_sfil_plots.C+(200000)'                  first 200k triggers, default flags
//   root -l -b -q 'kcu0_sfil_plots.C+(200000, 100000)'          200k triggers starting at trigger 100000
//   root -l 'kcu0_sfil_plots.C+'   then   Flags f; f.requirePaddle = 6; kcu0_sfil_plots(-1, 0, f);   all triggers, paddle P6 only
//
// Output: kcu0_sfil_plots.root (histograms hRaw/hPed/hTOA/hTOT_<layer>_<S|N>_c<cell> and _P3..P6, paddle dips and
// multiplicity, the cut flow hFlags, the error counters hErrors and hErrorsPerLink) and one figure per quantity with
// every channel laid out like the setup: kcu0_all_<raw|pedsub|toa|tot>.png, 12 columns (paddles P4 P3 P6 P5 along the
// module in front of SFIL1, then SFIL1-4 with North | South ends) x 4 rows (hookup cells 1 top .. 4 bottom, beam row = 3),
// paddle fire fractions and multiplicity in the caption. perLayer = true adds the 2 x 4 figures per layer.
#include "TFile.h"
#include "TTree.h"
#include "TH2F.h"
#include "TH1D.h"
#include "TCanvas.h"
#include "TPad.h"
#include "TStyle.h"
#include "TLatex.h"
#include "TString.h"
#include <algorithm>
#include <cstdio>

struct Flags
{
  bool requireComplete = true;    // 20 samples with all 8 links (header.complete)
  bool rejectCRC = true;          // any line of the trigger with a failed CRC-32        (header.crc_ok[s][l])
  bool rejectHamming = true;      // any line with H1-H3 Hamming flags set              (header.ham[s][l])
  bool rejectBackward = true;     // trigger out of time order                           (header.backward)
  bool rejectOverflow = true;     // pool offered more than 20 samples                   (header.nsamp_overflow)
  bool sfilLinksOnly = false;     // apply the CRC / Hamming cuts to the SFIL links (ASIC0, ASIC2) only, not to the paddles / ASIC3
  int  requirePaddle = -1;         // 0 = no paddle condition, -1 = any paddle fired, 3/4/5/6 = that paddle fired (dip > paddleDip)
  bool rejectMultiPaddle = true; // require exactly one paddle
  double paddleDip = 40.;         // paddle hit: pedestal (samples 0-3) minus the minimum of samples 8-16, ADC
  int  pedFirst = 1, pedLast = 3; // pedestal = mean of samples pedFirst..pedLast (sample 0 is left out)
  bool perLayer = true;          // also write the per-layer 2 x 4 figures
  // waveform-level cleaning (applied per channel to the pedestal-subtracted, TOA and TOT histograms; the raw ones take everything)
  double minAmp = 100.;             // keep a channel's waveform only if max(samples 4-12) - pedestal > minAmp (0 = off)
  double maxPedSpread = 30.;       // drop a waveform whose pedestal samples spread (max - min) by more than this (0 = off): pile-up tails, baseline jumps
  bool rejectSaturated = false;   // drop waveforms with a sample at 1023
  bool toaPulseOnly = true;      // TOA / TOT words only in the crossing sample and the next two (and only when the channel passes minAmp)
  bool normalise = false;         // also fill hNorm_*: (ADC - pedestal) / amplitude vs sample for the kept pulses
};

static const int SFIL_S[4][4] = {{2, 6, 10, 14}, {21, 25, 29, 33}, {63, 59, 52, 48}, {43, 37, 0, 4}};
static const int SFIL_N[4][4] = {{146, 150, 154, 158}, {165, 169, 173, 177}, {207, 203, 196, 192}, {187, 186, 144, 148}};
static int tot_decompress(int v) { return v < 512 ? v : 8 * (v - 512) + 4; }

void kcu0_sfil_plots(
  Long64_t nmax = 400000, 
  Long64_t start = 200000, 
  Flags f = Flags(), 
  const char *file = "../data/kcu0_3062.root", 
  const char *out = "kcu0_sfil_plots.root"
)
{
  TFile *fin = TFile::Open(file); if (!fin || fin->IsZombie()) { printf("cannot open %s\n", file); return; }
  TTree *h = (TTree *) fin->Get("header"), *pad = (TTree *) fin->Get("paddles"), *sf[4];
  for (int L = 0; L < 4; L++) sf[L] = (TTree *) fin->Get(Form("sfil%d", L + 1));
  // header flags
  Bool_t complete, backward, overflow, crc_ok[20][8], present[20][8]; UChar_t ham[20][8]; Int_t nlines;
  h->SetBranchStatus("*", 0); for (const char *b : {"complete", "backward", "nsamp_overflow", "crc_ok", "ham", "nlines", "present"}) h->SetBranchStatus(b, 1);
  h->SetBranchAddress("complete", &complete); h->SetBranchAddress("backward", &backward); h->SetBranchAddress("nsamp_overflow", &overflow); h->SetBranchAddress("crc_ok", crc_ok); h->SetBranchAddress("ham", ham); h->SetBranchAddress("nlines", &nlines); h->SetBranchAddress("present", present);
  // paddles
  Short_t P[4][20], Ptoa[4][20], Ptot[4][20]; const char *pn[4] = {"P3", "P4", "P5", "P6"};
  for (int k = 0; k < 4; k++) { pad->SetBranchAddress(Form("%s_adc", pn[k]), P[k]); pad->SetBranchAddress(Form("%s_toa", pn[k]), Ptoa[k]); pad->SetBranchAddress(Form("%s_tot", pn[k]), Ptot[k]); }
  // SFIL channels: S[layer][cell][sample], N[layer][cell][sample] for adc / toa / tot
  Short_t Sadc[4][4][20], Nadc[4][4][20], Stoa[4][4][20], Ntoa[4][4][20], Stot[4][4][20], Ntot[4][4][20];
  for (int L = 0; L < 4; L++) for (int c = 0; c < 4; c++)
  {
    sf[L]->SetBranchAddress(Form("Sc%d_adc", c + 1), Sadc[L][c]); sf[L]->SetBranchAddress(Form("Nc%d_adc", c + 1), Nadc[L][c]);
    sf[L]->SetBranchAddress(Form("Sc%d_toa", c + 1), Stoa[L][c]); sf[L]->SetBranchAddress(Form("Nc%d_toa", c + 1), Ntoa[L][c]);
    sf[L]->SetBranchAddress(Form("Sc%d_tot", c + 1), Stot[L][c]); sf[L]->SetBranchAddress(Form("Nc%d_tot", c + 1), Ntot[L][c]);
  }
  // histograms
  TFile *fout = new TFile(out, "RECREATE");
  TH2F *hRaw[4][2][4], *hPed[4][2][4], *hTOA[4][2][4], *hTOT[4][2][4], *hNorm[4][2][4]; const char *end[2] = {"S", "N"};
  for (int L = 0; L < 4; L++) for (int e = 0; e < 2; e++) for (int c = 0; c < 4; c++)
  {
    const int ch = e == 0 ? SFIL_S[L][c] : SFIL_N[L][c];
    hRaw[L][e][c] = new TH2F(Form("hRaw_%d_%s_c%d", L + 1, end[e], c + 1), Form("SFIL%d %s cell %d (ch %d) raw;sample (25 ns);ADC", L + 1, end[e], c + 1, ch), 20, -0.5, 19.5, 256, -0.5, 1023.5);
    hPed[L][e][c] = new TH2F(Form("hPed_%d_%s_c%d", L + 1, end[e], c + 1), Form("SFIL%d %s cell %d (ch %d) pedestal-subtracted;sample (25 ns);ADC - pedestal (samples %d-%d)", L + 1, end[e], c + 1, ch, f.pedFirst, f.pedLast), 20, -0.5, 19.5, 280, -100, 1020);
    hTOA[L][e][c] = new TH2F(Form("hTOA_%d_%s_c%d", L + 1, end[e], c + 1), Form("SFIL%d %s cell %d (ch %d) TOA;sample (25 ns);TOA code", L + 1, end[e], c + 1, ch), 20, -0.5, 19.5, 256, -0.5, 1023.5);
    hTOT[L][e][c] = new TH2F(Form("hTOT_%d_%s_c%d", L + 1, end[e], c + 1), Form("SFIL%d %s cell %d (ch %d) TOT;sample (25 ns);TOT (decompressed)", L + 1, end[e], c + 1, ch), 20, -0.5, 19.5, 256, -0.5, 4095.5);
    hNorm[L][e][c] = new TH2F(Form("hNorm_%d_%s_c%d", L + 1, end[e], c + 1), Form("SFIL%d %s cell %d (ch %d) normalised;sample (25 ns);(ADC - pedestal) / amplitude", L + 1, end[e], c + 1, ch), 20, -0.5, 19.5, 130, -0.2, 1.1);
  }
  TH2F *hRawP[4], *hPedP[4], *hTOAP[4], *hTOTP[4]; TH1D *hDip[4]; const int PCH[4] = {84, 86, 81, 97};
  for (int k = 0; k < 4; k++)
    {
      hRawP[k] = new TH2F(Form("hRaw_%s", pn[k]), Form("paddle %s (ch %d) raw;sample (25 ns);ADC", pn[k], PCH[k]), 20, -0.5, 19.5, 256, -0.5, 1023.5);
      hPedP[k] = new TH2F(Form("hPed_%s", pn[k]), Form("paddle %s (ch %d) pedestal-subtracted;sample (25 ns);ADC - pedestal", pn[k], PCH[k]), 20, -0.5, 19.5, 300, -150, 150);
      hTOAP[k] = new TH2F(Form("hTOA_%s", pn[k]), Form("paddle %s (ch %d) TOA;sample (25 ns);TOA code", pn[k], PCH[k]), 20, -0.5, 19.5, 256, -0.5, 1023.5);
      hTOTP[k] = new TH2F(Form("hTOT_%s", pn[k]), Form("paddle %s (ch %d) TOT;sample (25 ns);TOT (decompressed)", pn[k], PCH[k]), 20, -0.5, 19.5, 256, -0.5, 4095.5);
      hDip[k] = new TH1D(Form("hDip_%s", pn[k]), Form("paddle %s dip = pedestal - min(samples 8-16);ADC;triggers", pn[k]), 220, -20, 200);
    }
  TH1D *hMult = new TH1D("hPaddleMult", "paddles fired per kept trigger;paddles;triggers", 5, -0.5, 4.5);
  TH1D *hErr = new TH1D("hErrors", "error counters over all triggers looked at;;count", 10, 0, 10);
  const char *en[10] = {"triggers", "incomplete", "CRC-failed lines", "triggers with CRC fail", "Hamming lines", "triggers with Hamming", "missing lines", "triggers with missing lines", "out of order", "overflow"};
  for (int i = 0; i < 10; i++) hErr->GetXaxis()->SetBinLabel(i + 1, en[i]);
  TH1D *hErrLink = new TH1D("hErrorsPerLink", "CRC-failed (bins 0-7) and Hamming-flagged (bins 8-15) lines per link = 2 ASIC + half;;lines", 16, 0, 16);
  for (int l = 0; l < 8; l++) { hErrLink->GetXaxis()->SetBinLabel(l + 1, Form("CRC A%d.%d", l / 2, l % 2)); hErrLink->GetXaxis()->SetBinLabel(l + 9, Form("Ham A%d.%d", l / 2, l % 2)); }
  TH1D *hFlags = new TH1D("hFlags", "triggers passing each cut in turn;;triggers", 9, 0, 9);
  const char *fl[9] = {"all", "complete", "no CRC fail", "no Hamming", "in order", "no overflow", "paddle", "single paddle", "kept"}; for (int i = 0; i < 9; i++) hFlags->GetXaxis()->SetBinLabel(i + 1, fl[i]);

  const Long64_t total = h->GetEntries(); if (start < 0) start = 0; if (start >= total) { printf("start %lld is beyond the %lld triggers\n", start, total); return; }
  const Long64_t stop = nmax < 0 ? total : std::min(start + nmax, total), n = stop - start; Long64_t kept = 0;
  printf("triggers %lld to %lld of %lld\n", start, stop - 1, total);
  for (Long64_t i = start; i < stop; i++)
    {
      h->GetEntry(i); hFlags->Fill(0); hErr->Fill(0);
      // error counters, before any cut
      { int ncrc = 0, nham = 0; for (int s = 0; s < 20; s++) for (int l = 0; l < 8; l++) { if (!present[s][l]) continue; if (!crc_ok[s][l]) { ncrc++; hErrLink->Fill(l); } if (ham[s][l]) { nham++; hErrLink->Fill(8 + l); } }   // a missing line is counted as missing, not as a CRC failure
        if (!complete) hErr->Fill(1); hErr->Fill(2, ncrc); if (ncrc) hErr->Fill(3); hErr->Fill(4, nham); if (nham) hErr->Fill(5);
        const int miss = 160 - nlines; if (miss > 0) { hErr->Fill(6, miss); hErr->Fill(7); } if (backward) hErr->Fill(8); if (overflow) hErr->Fill(9); }
      if (f.requireComplete && !complete) continue; hFlags->Fill(1);
      bool badcrc = false, badham = false;
      for (int s = 0; s < 20; s++) for (int l = 0; l < 8; l++) { if (!present[s][l]) continue; if (f.sfilLinksOnly && (l == 2 || l == 3 || l == 6 || l == 7)) continue; if (!crc_ok[s][l]) badcrc = true; if (ham[s][l]) badham = true; }
      if (f.rejectCRC && badcrc) continue; hFlags->Fill(2);
      if (f.rejectHamming && badham) continue; hFlags->Fill(3);
      if (f.rejectBackward && backward) continue; hFlags->Fill(4);
      if (f.rejectOverflow && overflow) continue; hFlags->Fill(5);
      pad->GetEntry(i); int npad = 0; bool fired[4]; double pped[4];
      for (int k = 0; k < 4; k++) { double ped = 0; for (int s = f.pedFirst; s <= f.pedLast; s++) ped += P[k][s]; ped /= (f.pedLast - f.pedFirst + 1); pped[k] = ped; int mn = 1023; for (int s = 8; s <= 16; s++) if (P[k][s] >= 0 && P[k][s] < mn) mn = P[k][s]; fired[k] = ped - mn > f.paddleDip; if (fired[k]) npad++; }
      if (f.requirePaddle == -1 && npad == 0) continue;
      if (f.requirePaddle >= 3 && !fired[f.requirePaddle - 3]) continue; hFlags->Fill(6);
      if (f.rejectMultiPaddle && npad != 1) continue; hFlags->Fill(7);
      hMult->Fill(npad);
      for (int k = 0; k < 4; k++)
        {
          if (P[k][0] < 0) continue;
          int mn = 1023; for (int s = 8; s <= 16; s++) if (P[k][s] >= 0 && P[k][s] < mn) mn = P[k][s]; hDip[k]->Fill(pped[k] - mn);
          for (int s = 0; s < 20; s++) { if (P[k][s] < 0) continue; hRawP[k]->Fill(s, P[k][s]); hPedP[k]->Fill(s, P[k][s] - pped[k]); if (Ptoa[k][s] > 0) hTOAP[k]->Fill(s, Ptoa[k][s]); if (Ptot[k][s] > 0) hTOTP[k]->Fill(s, tot_decompress(Ptot[k][s])); }
        }
      hFlags->Fill(8); kept++;
      for (int L = 0; L < 4; L++)
        {
          sf[L]->GetEntry(i);
          for (int e = 0; e < 2; e++) for (int c = 0; c < 4; c++)
            {
              const Short_t *adc = e == 0 ? Sadc[L][c] : Nadc[L][c], *toa = e == 0 ? Stoa[L][c] : Ntoa[L][c], *tot = e == 0 ? Stot[L][c] : Ntot[L][c];
              if (adc[0] < 0) continue;   // line missing in sample 0 (only possible when requireComplete is off)
              double ped = 0; int pmin = 1023, pmax = 0; for (int s = f.pedFirst; s <= f.pedLast; s++) { ped += adc[s]; pmin = std::min(pmin, (int) adc[s]); pmax = std::max(pmax, (int) adc[s]); } ped /= (f.pedLast - f.pedFirst + 1);
              for (int s = 0; s < 20; s++) if (adc[s] >= 0) hRaw[L][e][c]->Fill(s, adc[s]);
              // cleaning
              double amp = 0; int sat = 0, sx = -1; for (int s = 4; s < 20; s++) { if (adc[s] == 1023) sat++; if (s <= 12 && adc[s] - ped > amp) amp = adc[s] - ped; }
              if (f.minAmp > 0 && amp <= f.minAmp) continue;
              if (f.maxPedSpread > 0 && pmax - pmin > f.maxPedSpread) continue;
              if (f.rejectSaturated && sat) continue;
              if (amp > 0) for (int s = 1; s < 20; s++) if (adc[s] - ped >= 0.2 * amp) { sx = s; break; }   // crossing sample: first above 20% of the amplitude
              for (int s = 0; s < 20; s++)
                {
                  if (adc[s] < 0) continue;
                  hPed[L][e][c]->Fill(s, adc[s] - ped);
                  if (f.normalise && amp > 0) hNorm[L][e][c]->Fill(s, (adc[s] - ped) / amp);
                  const bool inwin = !f.toaPulseOnly || (sx >= 0 && s >= sx && s <= sx + 2 && (f.minAmp > 0 || amp > 50));
                  if (toa[s] > 0 && inwin) hTOA[L][e][c]->Fill(s, toa[s]);
                  if (tot[s] > 0 && inwin) hTOT[L][e][c]->Fill(s, tot_decompress(tot[s]));
                }
            }
        }
    }
  printf("error counters over %lld triggers:\n", n);
  for (int i = 0; i < 10; i++) printf("  %-28s %10.0f\n", en[i], hErr->GetBinContent(i + 1));
  printf("  per link (ASIC.half)  CRC-failed lines:"); for (int l = 0; l < 8; l++) printf(" A%d.%d:%.0f", l / 2, l % 2, hErrLink->GetBinContent(l + 1)); printf("\n  per link (ASIC.half)  Hamming lines:   "); for (int l = 0; l < 8; l++) printf(" A%d.%d:%.0f", l / 2, l % 2, hErrLink->GetBinContent(l + 9)); printf("\n");
  printf("cut flow: %lld of %lld triggers kept\n", kept, n);
  for (int i = 0; i < 9; i++) printf("  %-14s %10.0f\n", fl[i], hFlags->GetBinContent(i + 1));
  // figures
  gStyle->SetOptStat(0); gStyle->SetOptTitle(0); gStyle->SetPalette(57); gStyle->SetPadLeftMargin(0.16); gStyle->SetPadBottomMargin(0.15); gStyle->SetPadRightMargin(0.04); gStyle->SetPadTopMargin(0.1); gStyle->SetTitleSize(0.06, "xy"); gStyle->SetLabelSize(0.05, "xy"); gStyle->SetTitleOffset(1.2, "y"); gStyle->SetNdivisions(505, "xy");
  const int NQ = f.normalise ? 5 : 4; const char *qn[5] = {"raw", "pedsub", "toa", "tot", "norm"};
  TString clean = ""; if (f.minAmp > 0) clean += Form("; amplitude > %.0f", f.minAmp); if (f.maxPedSpread > 0) clean += Form("; pedestal spread < %.0f", f.maxPedSpread); if (f.rejectSaturated) clean += "; no saturated sample"; if (f.toaPulseOnly) clean += "; TOA/TOT words in the crossing sample +0..+2 only";
  const TString qdesc[5] = {"raw ADC vs sample", TString::Format("ADC minus pedestal (samples %d-%d) vs sample", f.pedFirst, f.pedLast) + clean, "TOA code vs sample (non-zero words)" + clean, "decompressed TOT vs sample (non-zero words)" + clean, "(ADC - pedestal) / amplitude vs sample" + clean};
  TString padinfo = "paddles fired per kept trigger: ";
  for (int k = 0; k < 4; k++) padinfo += Form("%s %.1f%%  ", pn[k], 100. * hDip[k]->Integral(hDip[k]->FindBin(f.paddleDip), hDip[k]->GetNbinsX() + 1) / std::max(1LL, kept));
  padinfo += Form("| multiplicity 0/1/2/3+: %.1f/%.1f/%.1f/%.1f%%", 100. * hMult->GetBinContent(1) / std::max(1LL, kept), 100. * hMult->GetBinContent(2) / std::max(1LL, kept), 100. * hMult->GetBinContent(3) / std::max(1LL, kept), 100. * (hMult->GetBinContent(4) + hMult->GetBinContent(5)) / std::max(1LL, kept));
  auto panel = [&](TH2F *hh, const char *title) { gPad->SetLogz(); hh->Draw("col"); TLatex t; t.SetNDC(); t.SetTextSize(0.085); t.DrawLatex(0.18, 0.92, title); };
  // one figure per quantity laid out like the setup: beam from the left; the four paddles along the module (North to
  // South: P4 P3 P6 P5) in front of SFIL1 at the beam row; then SFIL1..4 in beam order, each with its North end on
  // the left and its South end on the right; rows = cells 1 (top) .. 4 (bottom) of the hookup table (beam row = cell 3).
  const int PORD[4] = {1, 0, 3, 2};   // paddle panels North -> South: P4, P3, P6, P5 (k index into pn / PCH)
  for (int q = 0; q < NQ; q++)
    {
      TCanvas *c = new TCanvas(Form("call_%s", qn[q]), "", 3600, 1450); TPad *band = new TPad("band", "", 0, 0.93, 1, 1); band->Draw(); TPad *grid = new TPad("grid", "", 0, 0, 1, 0.93); grid->Draw(); grid->Divide(12, 4, 0.002, 0.004);
      for (int k = 0; k < 4; k++) { grid->cd(2 * 12 + k + 1); const int kk = PORD[k]; TH2F *hh = q == 0 ? hRawP[kk] : q == 1 ? hPedP[kk] : q == 2 ? hTOAP[kk] : q == 3 ? hTOTP[kk] : nullptr; if (hh) panel(hh, Form("paddle %s  (ch %d)", pn[kk], PCH[kk])); }
            for (int L = 0; L < 4; L++) for (int e = 0; e < 2; e++) for (int cell = 0; cell < 4; cell++)
        {
          const int col = 4 + 2 * L + (e == 0 ? 1 : 0);   // North end left, South end right
          grid->cd(cell * 12 + col + 1);
          TH2F *hh = q == 0 ? hRaw[L][e][cell] : q == 1 ? hPed[L][e][cell] : q == 2 ? hTOA[L][e][cell] : q == 3 ? hTOT[L][e][cell] : hNorm[L][e][cell];
          panel(hh, Form("SFIL%d %s cell %d  (ch %d)", L + 1, e == 0 ? "S" : "N", cell + 1, e == 0 ? SFIL_S[L][cell] : SFIL_N[L][cell]));
        }
      band->cd(); TLatex t; t.SetNDC(); t.SetTextSize(0.2); t.SetTextAlign(12);
      t.DrawLatex(0.003, 0.85, Form("run 3062, triggers %lld-%lld, %lld kept of %lld; %s", start, stop - 1, kept, n, qdesc[q].Data()));
      t.DrawLatex(0.003, 0.62, "layout: beam from the left; paddles along the module North to South (P4 P3 P6 P5) in front of SFIL1 at the beam row; SFIL1-4 in beam order with the North end left and the South end right; rows = hookup cells 1 (top) to 4 (bottom), beam row = cell 3, cells 1+2 are one row on v2 boards");
      t.DrawLatex(0.003, 0.4, padinfo);
      t.SetTextSize(0.26); t.SetTextAlign(22); t.SetTextColor(kBlue + 2);
      const char *hdr[5] = {"trigger paddles  (N #rightarrow S)", "SFIL1  (N | S)", "SFIL2  (N | S)", "SFIL3  (N | S)", "SFIL4  (N | S)"}; const double hx[5] = {2. / 12, 5. / 12, 7. / 12, 9. / 12, 11. / 12};
      for (int b = 0; b < 5; b++) t.DrawLatex(hx[b], 0.12, hdr[b]);
      c->SaveAs(Form("kcu0_all_%s.png", qn[q]));
    }
  if (f.perLayer)
    for (int L = 0; L < 4; L++) for (int q = 0; q < 4; q++)
      {
        TCanvas *c = new TCanvas(Form("c_%d_%s", L + 1, qn[q]), "", 1600, 830); TPad *band = new TPad("band", "", 0, 0.964, 1, 1); band->Draw(); TPad *grid = new TPad("grid", "", 0, 0, 1, 0.964); grid->Draw(); grid->Divide(4, 2, 0.003, 0.003);
        for (int e = 0; e < 2; e++) for (int cell = 0; cell < 4; cell++)
          {
            grid->cd(4 * e + cell + 1);
            TH2F *hh = q == 0 ? hRaw[L][e][cell] : q == 1 ? hPed[L][e][cell] : q == 2 ? hTOA[L][e][cell] : hTOT[L][e][cell];
            panel(hh, Form("SFIL%d %s cell %d  (ch %d)", L + 1, e == 0 ? "South" : "North", cell + 1, e == 0 ? SFIL_S[L][cell] : SFIL_N[L][cell]));
          }
        band->cd(); TLatex t; t.SetNDC(); t.SetTextSize(0.5); t.SetTextAlign(12); t.DrawLatex(0.01, 0.5, Form("run 3062, triggers %lld-%lld, %lld kept of %lld; %s", start, stop - 1, kept, n, qdesc[q].Data()));
        c->SaveAs(Form("kcu0_sfil%d_%s.png", L + 1, qn[q]));
      }
  fout->cd(); fout->Write(); printf("wrote %s and kcu0_all_<raw|pedsub|toa|tot>.png%s\n", out, f.perLayer ? " and kcu0_sfil<layer>_<quantity>.png" : "");
}
