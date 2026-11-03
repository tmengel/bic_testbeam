// kcu1_sfil5_raw.C: raw-waveform, TOA and TOT galleries of both SFIL5 ends on KCU1 (run 3062), from kcu1_3062.root.
//
//   root -l -b -q 'kcu1_sfil5_raw.C+("../data/kcu1_3062.root", 100000)'
//
// SFIL5 South (tree sfil5, ASIC2): 4 arrays A B C D of 16 single-SiPM channels, drawn as 4 rows x 16 columns (rows C, D, B, A top to bottom,
// their order in height); SFIL5 North (tree sfil5north, ASIC0 ch 0-71 and ASIC1 ch 72-143): 72 channels per ASIC, 6 rows x 12 columns.
// Per channel: pedestal-subtracted waveform (sample x ADC), TOA code x sample, TOT (decompressed) x sample; all at log z.
//   kcu1_south_{wf,toa,tot}.png   kcu1_north_asic0_{wf,toa,tot}.png   kcu1_north_asic1_{wf,toa,tot}.png
//   kcu1_south_toa_time.png: for pulses > 100 ADC, the TOA time t = 25 x sample + TOA_ns of the first TOA word minus 25 x (peak sample):
//                            the pulse-timed TOA sits in a narrow band; words on pedestal samples spread over the window.
// Cuts: complete triggers, no backward/overflow timestamps; CRC + Hamming on ASIC2 for the South end; Hamming only on ASIC0/1 for the North end (its CRC fails on most triggers).
#include "TFile.h"
#include "TTree.h"
#include "TH1D.h"
#include "TH2F.h"
#include "TCanvas.h"
#include "TPad.h"
#include "TStyle.h"
#include "TLatex.h"
#include "TString.h"
#include <cmath>
#include <cstdio>

static int totdec(int v) { return v < 512 ? v : 8 * (v - 512) + 4; }
static double toans(int v) { return (v & 7) * 0.025 + ((v >> 3) & 31) * 0.2 + ((v >> 8) & 3) * 6.25; }

void kcu1_sfil5_raw(const char *rawfile = "../data/kcu1_3062.root", Long64_t nraw = 100000)
{
  gStyle->SetOptStat(0); gStyle->SetOptTitle(0); gStyle->SetPalette(57); gStyle->SetPadLeftMargin(0.16); gStyle->SetPadBottomMargin(0.15); gStyle->SetPadRightMargin(0.04); gStyle->SetPadTopMargin(0.1); gStyle->SetTitleSize(0.065, "xy"); gStyle->SetLabelSize(0.055, "xy"); gStyle->SetNdivisions(505, "xy"); gStyle->SetTitleOffset(1.1, "y");
  TFile *fr = TFile::Open(rawfile); if (!fr || fr->IsZombie()) { printf("cannot open %s\n", rawfile); return; }
  TTree *h = (TTree *) fr->Get("header"), *t5 = (TTree *) fr->Get("sfil5"), *tn = (TTree *) fr->Get("sfil5north");
  Bool_t complete, backward, overflow, crc_ok[20][8], present[20][8]; UChar_t ham[20][8];
  h->SetBranchStatus("*", 0); for (const char *b : {"complete", "backward", "nsamp_overflow", "crc_ok", "ham", "present"}) h->SetBranchStatus(b, 1);
  h->SetBranchAddress("complete", &complete); h->SetBranchAddress("backward", &backward); h->SetBranchAddress("nsamp_overflow", &overflow); h->SetBranchAddress("crc_ok", crc_ok); h->SetBranchAddress("ham", ham); h->SetBranchAddress("present", present);
  const char *an[4] = {"A", "B", "C", "D"}; Int_t CH[4][16]; Short_t W[4][16][20], TA[4][16][20], TT[4][16][20];
  for (int a = 0; a < 4; a++) { t5->SetBranchAddress(Form("%s_ch", an[a]), CH[a]); t5->SetBranchAddress(Form("%s_adc", an[a]), W[a]); t5->SetBranchAddress(Form("%s_toa", an[a]), TA[a]); t5->SetBranchAddress(Form("%s_tot", an[a]), TT[a]); }
  Int_t NCH[144]; Short_t NW[144][20], NTA[144][20], NTT[144][20]; tn->SetBranchAddress("ch", NCH); tn->SetBranchAddress("adc", NW); tn->SetBranchAddress("toa", NTA); tn->SetBranchAddress("tot", NTT);
  TH2F *sW[64], *sA[64], *sT[64], *nW[144], *nA[144], *nT[144]; TH1D *sTime[64];
  for (int k = 0; k < 64; k++) { sW[k] = new TH2F(Form("sW%d", k), ";sample;ADC - pedestal", 20, -0.5, 19.5, 120, -100, 1100); sA[k] = new TH2F(Form("sA%d", k), ";sample;TOA code", 20, -0.5, 19.5, 128, 0, 1024); sT[k] = new TH2F(Form("sT%d", k), ";sample;TOT (decompressed)", 20, -0.5, 19.5, 128, 0, 4096); sTime[k] = new TH1D(Form("sTime%d", k), ";t_{TOA} - 25 #times peak sample [ns];pulses", 120, -300, 300); }
  for (int k = 0; k < 144; k++) { nW[k] = new TH2F(Form("nW%d", k), ";sample;ADC - pedestal", 20, -0.5, 19.5, 120, -100, 1100); nA[k] = new TH2F(Form("nA%d", k), ";sample;TOA code", 20, -0.5, 19.5, 128, 0, 1024); nT[k] = new TH2F(Form("nT%d", k), ";sample;TOT (decompressed)", 20, -0.5, 19.5, 128, 0, 4096); }
  Long64_t usedS = 0, usedN = 0; const Long64_t stop = std::min(nraw, h->GetEntries());
  for (Long64_t i = 0; i < stop; i++)
    {
      h->GetEntry(i); if (!complete || backward || overflow) continue;
      bool badS = false, badN = false;
      for (int s = 0; s < 20; s++) { for (int l : {4, 5}) if (present[s][l] && (!crc_ok[s][l] || ham[s][l])) badS = true; for (int l : {0, 1, 2, 3}) if (present[s][l] && ham[s][l]) badN = true; }   // North: no CRC requirement (half 1 of ASIC0/1 fails CRC on 52 / 12 % of lines, which would remove nearly every trigger)
      if (!badS)
        {
          t5->GetEntry(i); usedS++;
          for (int k = 0; k < 64; k++)
            {
              const Short_t *w = W[k / 16][k % 16], *ta = TA[k / 16][k % 16], *tt = TT[k / 16][k % 16]; if (w[0] < 0) continue;
              const double ped = (w[1] + w[2] + w[3]) / 3.; double a = -1e9; int ipk = -1; int fs = -1, fcode = 0;
              for (int s = 0; s < 20; s++) { sW[k]->Fill(s, w[s] - ped); if (ta[s] > 0) { sA[k]->Fill(s, ta[s]); if (fs < 0) { fs = s; fcode = ta[s]; } } if (tt[s] > 0) sT[k]->Fill(s, totdec(tt[s])); if (s >= 4 && s <= 12 && w[s] - ped > a) { a = w[s] - ped; ipk = s; } }
              if (a > 100 && fs >= 0) sTime[k]->Fill(25. * fs + toans(fcode) - 25. * ipk);
            }
        }
      if (!badN)
        {
          tn->GetEntry(i); usedN++;
          for (int k = 0; k < 144; k++)
            {
              const Short_t *w = NW[k], *ta = NTA[k], *tt = NTT[k]; if (w[0] < 0) continue; const double ped = (w[1] + w[2] + w[3]) / 3.;
              for (int s = 0; s < 20; s++) { nW[k]->Fill(s, w[s] - ped); if (ta[s] > 0) nA[k]->Fill(s, ta[s]); if (tt[s] > 0) nT[k]->Fill(s, totdec(tt[s])); }
            }
        }
    }
  t5->GetEntry(0); tn->GetEntry(0);
  const int rowOf[4] = {3, 2, 0, 1};
  auto gridS = [&](const char *name, auto draw) {
    TCanvas *c = new TCanvas(name, "", 4000, 1200); TPad *band = new TPad("band", "", 0, 0.95, 1, 1); band->Draw(); TPad *g = new TPad("grid", "", 0, 0, 1, 0.95); g->Draw(); g->Divide(16, 4, 0.0015, 0.004);
    for (int a = 0; a < 4; a++) for (int i = 0; i < 16; i++) { const int k = 16 * a + i; g->cd(rowOf[a] * 16 + i + 1); draw(k); TLatex t; t.SetNDC(); t.SetTextSize(0.09); t.DrawLatex(0.2, 0.91, Form("%s%d ch %d", an[a], i, CH[a][i])); }
    return std::make_pair(c, band); };
  auto gridN = [&](const char *name, int asic, auto draw) {
    TCanvas *c = new TCanvas(name, "", 3600, 1900); TPad *band = new TPad("band", "", 0, 0.96, 1, 1); band->Draw(); TPad *g = new TPad("grid", "", 0, 0, 1, 0.96); g->Draw(); g->Divide(12, 6, 0.0015, 0.003);
    for (int j = 0; j < 72; j++) { const int k = 72 * asic + j; g->cd(j + 1); draw(k); TLatex t; t.SetNDC(); t.SetTextSize(0.09); t.DrawLatex(0.2, 0.91, Form("ch %d", NCH[k])); }
    return std::make_pair(c, band); };
  auto cap = [&](TPad *band, const char *what, double size = 0.32) { band->cd(); TLatex t; t.SetNDC(); t.SetTextSize(size); t.SetTextAlign(12); t.DrawLatex(0.003, 0.5, what); };
  { auto cb = gridS("cSW", [&](int k) { gPad->SetLogz(); sW[k]->Draw("col"); }); cap(cb.second, Form("run 3062, KCU1 SFIL5 South (ASIC2, 4 arrays x 16 single SiPMs; rows top to bottom C, D, B, A): pedestal-subtracted waveforms of %lld clean triggers, log z", usedS)); cb.first->SaveAs("kcu1_south_wf.png"); }
  { auto cb = gridS("cSA", [&](int k) { gPad->SetLogz(); sA[k]->Draw("col"); }); cap(cb.second, Form("run 3062, KCU1 SFIL5 South: TOA code vs sample (words > 0 only), %lld clean triggers; a pulse-timed TOA appears one sample after the crossing", usedS)); cb.first->SaveAs("kcu1_south_toa.png"); }
  { auto cb = gridS("cST", [&](int k) { gPad->SetLogz(); sT[k]->Draw("col"); }); cap(cb.second, Form("run 3062, KCU1 SFIL5 South: decompressed TOT vs sample (words > 0 only), %lld clean triggers; TOT fires on pulses near and above the ceiling", usedS)); cb.first->SaveAs("kcu1_south_tot.png"); }
  { auto cb = gridS("cSt", [&](int k) { sTime[k]->SetLineWidth(2); sTime[k]->Draw(); }); cap(cb.second, "run 3062, KCU1 SFIL5 South: TOA time of the first TOA word minus 25 ns x peak sample, pulses > 100 ADC (forward convention t = 25 s + TOA_ns); the narrow band is the pulse-timed TOA, flat parts are words on pedestal samples"); cb.first->SaveAs("kcu1_south_toa_time.png"); }
  for (int asic = 0; asic < 2; asic++)
    {
      { auto cb = gridN(Form("cNW%d", asic), asic, [&](int k) { gPad->SetLogz(); nW[k]->Draw("col"); }); cap(cb.second, Form("run 3062, KCU1 SFIL5 North, ASIC%d (channels %d-%d): pedestal-subtracted waveforms of %lld triggers (complete, no Hamming flag on ASIC0/1; CRC NOT required here: half 1 of both ASICs fails it on 52 / 12 %% of lines); half 0 (ch %d-%d) is railed at 1023 in every trigger", asic, 72 * asic, 72 * asic + 71, usedN, 72 * asic, 72 * asic + 35), 0.3); cb.first->SaveAs(Form("kcu1_north_asic%d_wf.png", asic)); }
      { auto cb = gridN(Form("cNA%d", asic), asic, [&](int k) { gPad->SetLogz(); nA[k]->Draw("col"); }); cap(cb.second, Form("run 3062, KCU1 SFIL5 North, ASIC%d: TOA code vs sample (words > 0 only)", asic), 0.3); cb.first->SaveAs(Form("kcu1_north_asic%d_toa.png", asic)); }
      { auto cb = gridN(Form("cNT%d", asic), asic, [&](int k) { gPad->SetLogz(); nT[k]->Draw("col"); }); cap(cb.second, Form("run 3062, KCU1 SFIL5 North, ASIC%d: decompressed TOT vs sample (words > 0 only)", asic), 0.3); cb.first->SaveAs(Form("kcu1_north_asic%d_tot.png", asic)); }
    }
  printf("South: %lld clean triggers, North: %lld; wrote kcu1_south_{wf,toa,tot,toa_time}.png and kcu1_north_asic{0,1}_{wf,toa,tot}.png\n", usedS, usedN);
}
