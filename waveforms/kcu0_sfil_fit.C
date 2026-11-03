// kcu0_sfil_fit.C: pedestal-subtracted pulse fits of the KCU0 SFIL channels (run 3062, and any run via the file argument) and correlations of the fitted amplitudes.
//
//   root -l -b -q 'kcu0_sfil_fit.C+(20000, 200000)'      20k triggers from trigger 200000
//
// Per trigger (quality flags as in kcu0_sfil_plots.C) and per unmasked SFIL channel:
//   pedestal = mean of samples pedFirst..pedLast; amplitude = max(samples 4-12) - pedestal;
//   keep if amplitude > minAmp and the peak sample is in [peakLo, peakHi] (default 7-9), else the waveform is tossed;
//   fit a per-channel pulse TEMPLATE (built in a first pass from clean, unsaturated pulses of that channel, aligned on the
//   20% crossing time) with two free parameters, amplitude A and time t0, by a chi2 scan over t0 with A solved
//   analytically; saturated samples are masked (the fit extrapolates through them); sigma per sample =
//   sqrt(noise^2 + (shapeErr x model)^2); keep if chi2/ndf < maxChi2. Bad towers (default: the noisy SFIL1 North
//   column) are masked.
// Output: kcu0_sfil_fit.root with a tree "fits" (per trigger: ok/amp/afit/t0/chi2/ped per channel index k = 8 x layer +
// 4 x end + cell, end 0 = South, 1 = North), the templates (hTemplate_k) and histograms, plus figures kcu0_fit_*.png.
#include "TFile.h"
#include "TTree.h"
#include "TH1D.h"
#include "TH2F.h"
#include "TCanvas.h"
#include "TPad.h"
#include "TStyle.h"
#include "TLatex.h"
#include "TString.h"
#include "TSystem.h"
#include "TMath.h"
#include <set>
#include <vector>
#include <algorithm>
#include <cmath>
#include <cstdio>

struct FitFlags
{
  bool requireComplete = true, rejectCRC = true, rejectHamming = true, rejectBackward = true, rejectOverflow = true, sfilLinksOnly = true;
  int requirePaddle = -1; bool rejectMultiPaddle = true; double paddleDip = 40.;
  int pedFirst = 1, pedLast = 3;
  double minAmp = 60.;            // ADC above pedestal, else the waveform is tossed
  int peakLo = 7, peakHi = 9;     // peak sample window
  double maxChi2 = 8.;            // chi2 / ndf
  double noise = 5.;              // ADC, per-sample noise for the chi2
  double shapeErr = 0.05;         // relative shape uncertainty of the template, added in quadrature per sample (the tail fraction varies pulse to pulse)
  bool fitSaturated = true;       // fit saturated waveforms with the 1023 samples masked (the fit then extrapolates the amplitude)
  Long64_t nTemplate = 100000;    // triggers used in the first pass to build the templates (fast, no fits)
  int minTemplatePulses = 100;    // a channel needs at least this many clean pulses for a template, else it is skipped
  const char *templateFile = "";  // fit output of another run: channels with too few clean pulses take its hTemplate_k (low-rate runs 3065, 3069)
  std::set<int> mask = {146, 150, 154, 158};   // towers to ignore: the noisy SFIL1 North column
};
static const int SFIL_S[4][4] = {{2, 6, 10, 14}, {21, 25, 29, 33}, {63, 59, 52, 48}, {43, 37, 0, 4}};
static const int SFIL_N[4][4] = {{146, 150, 154, 158}, {165, 169, 173, 177}, {207, 203, 196, 192}, {187, 186, 144, 148}};
// 20% crossing time in ns (linear interpolation between the last sample below and the first at or above 0.2 x amplitude)
static double crossing(const double *y, double a, int ipk) { for (int s = 1; s <= ipk; s++) if (y[s] >= 0.2 * a) return 25. * (s - 1 + (0.2 * a - y[s - 1]) / std::max(1e-6, y[s] - y[s - 1])); return -1; }
// template: mean normalised pulse versus time relative to the crossing, 1 ns bins from -100 to +400 ns
struct Template { std::vector<double> sum, n; Template() : sum(500, 0), n(500, 0) {} double eval(double u) const { const double b = u + 100.; if (b < 0 || b >= 499) return 0; const int i = int(b); const double f = b - i; const double v0 = n[i] ? sum[i] / n[i] : 0, v1 = n[i + 1] ? sum[i + 1] / n[i + 1] : v0; return v0 + f * (v1 - v0); } double count() const { double m = 0; for (double v : n) m = std::max(m, v); return m; } };

void kcu0_sfil_fit(Long64_t nmax = 20000, Long64_t start = 0, FitFlags f = FitFlags(), const char *file = "../data/kcu0_3062.root", const char *out = "kcu0_sfil_fit.root")
{
  TFile *fin = TFile::Open(file); if (!fin || fin->IsZombie()) { printf("cannot open %s\n", file); return; }
  TTree *h = (TTree *) fin->Get("header"), *pad = (TTree *) fin->Get("paddles"), *sf[4]; for (int L = 0; L < 4; L++) sf[L] = (TTree *) fin->Get(Form("sfil%d", L + 1));
  Bool_t complete, backward, overflow, crc_ok[20][8], present[20][8]; UChar_t ham[20][8];
  h->SetBranchStatus("*", 0); for (const char *b : {"complete", "backward", "nsamp_overflow", "crc_ok", "ham", "present"}) h->SetBranchStatus(b, 1);
  h->SetBranchAddress("complete", &complete); h->SetBranchAddress("backward", &backward); h->SetBranchAddress("nsamp_overflow", &overflow); h->SetBranchAddress("crc_ok", crc_ok); h->SetBranchAddress("ham", ham); h->SetBranchAddress("present", present);
  Short_t P[4][20]; const char *pn[4] = {"P3", "P4", "P5", "P6"}; for (int k = 0; k < 4; k++) pad->SetBranchAddress(Form("%s_adc", pn[k]), P[k]);
  Short_t W[4][2][4][20];   // [layer][end][cell][sample]
  for (int L = 0; L < 4; L++) for (int c = 0; c < 4; c++) { sf[L]->SetBranchAddress(Form("Sc%d_adc", c + 1), W[L][0][c]); sf[L]->SetBranchAddress(Form("Nc%d_adc", c + 1), W[L][1][c]); }

  TFile *fout = new TFile(out, "RECREATE");
  // tree
  Long64_t trig; Int_t padid; Bool_t ok[32]; Float_t amp[32], afit[32], t0[32], chi2[32], pedv[32], taur[32], tauf[32], pw[32]; Int_t chan[32], masked[32];
  TTree *T = new TTree("fits", "pulse fits per trigger; index k = 8 layer + 4 end + cell (layer 0-3, end 0 South 1 North, cell 0-3)");
  T->Branch("trig", &trig, "trig/L"); T->Branch("pad", &padid, "pad/I"); T->Branch("chan", chan, "chan[32]/I"); T->Branch("masked", masked, "masked[32]/I");
  T->Branch("ok", ok, "ok[32]/O"); T->Branch("ped", pedv, "ped[32]/F"); T->Branch("amp", amp, "amp[32]/F"); T->Branch("afit", afit, "afit[32]/F"); T->Branch("t0", t0, "t0[32]/F"); T->Branch("chi2", chi2, "chi2[32]/F");
  for (int L = 0; L < 4; L++) for (int e = 0; e < 2; e++) for (int c = 0; c < 4; c++) { const int k = 8 * L + 4 * e + c, ch = e == 0 ? SFIL_S[L][c] : SFIL_N[L][c]; chan[k] = ch; masked[k] = f.mask.count(ch) ? 1 : 0; }
  // histograms
  TH1D *hA[32], *hChi[32], *hT0[32], *hRej = new TH1D("hReasons", "waveforms per outcome;;waveforms", 6, 0, 6); TH2F *hAA[32];
  const char *rn[6] = {"looked at", "amp < minAmp", "peak outside window", "fit failed", "chi2 too large", "kept"}; for (int i = 0; i < 6; i++) hRej->GetXaxis()->SetBinLabel(i + 1, rn[i]);
  for (int k = 0; k < 32; k++)
    {
      const int L = k / 8, e = (k / 4) % 2, c = k % 4;
      hA[k] = new TH1D(Form("hAfit_%d", k), Form("SFIL%d %s cell %d (ch %d) fitted amplitude;A_{fit} [ADC];waveforms", L + 1, e ? "N" : "S", c + 1, chan[k]), 150, 0, 3000);
      hChi[k] = new TH1D(Form("hChi2_%d", k), Form("SFIL%d %s cell %d chi2/ndf;#chi^{2}/ndf;waveforms", L + 1, e ? "N" : "S", c + 1), 100, 0, 20);
      hT0[k] = new TH1D(Form("hT0_%d", k), Form("SFIL%d %s cell %d t0;t_{0} [ns];waveforms", L + 1, e ? "N" : "S", c + 1), 100, 50, 200);
      hAA[k] = new TH2F(Form("hAfitVsAmp_%d", k), Form("SFIL%d %s cell %d;max - pedestal [ADC];A_{fit} [ADC]", L + 1, e ? "N" : "S", c + 1), 100, 0, 1100, 100, 0, 3000);
    }
  TH2F *hSN[4] = {nullptr}; for (int L = 0; L < 4; L++) hSN[L] = new TH2F(Form("hSN_%d", L + 1), Form("SFIL%d beam row: North vs South fitted amplitude;A_{fit} South [ADC];A_{fit} North [ADC]", L + 1), 100, 0, 3000, 100, 0, 3000);
  TH2F *hLL[3]; for (int L = 0; L < 3; L++) hLL[L] = new TH2F(Form("hLL_%d%d", L + 1, L + 2), Form("beam row, #sqrt{S#timesN} fitted: SFIL%d vs SFIL%d;SFIL%d [ADC];SFIL%d [ADC]", L + 2, L + 1, L + 1, L + 2), 100, 0, 3000, 100, 0, 3000);
  TH2F *hAsym = new TH2F("hAsym", "beam row (N - S) / (N + S) of the fitted amplitudes per layer;layer;(N - S) / (N + S)", 4, 0.5, 4.5, 100, -1, 1);
  TH1D *hSum = new TH1D("hSum", "sum over SFIL1-4 beam row of #sqrt{S#timesN} fitted;ADC;triggers", 150, 0, 9000);

  const Long64_t total = h->GetEntries(); if (start < 0) start = 0; const Long64_t stop = nmax < 0 ? total : std::min(start + nmax, total);
  Long64_t kept = 0; Template tpl[32]; bool hasTpl[32] = {false}, borrowed[32] = {false};
  for (int pass = 0; pass < 2; pass++)
  {
  if (pass == 1)
    {
      // a channel without enough clean pulses takes the template of templateFile if given (marked f), else (off-row cells) borrows
      // the template of the beam-row cell of the same layer and end
      bool fromFile[32] = {false};
      if (f.templateFile && *f.templateFile)
        {
          TFile *ftpl = TFile::Open(f.templateFile); if (!ftpl || ftpl->IsZombie()) { printf("cannot open template file %s\n", f.templateFile); return; }
          for (int k = 0; k < 32; k++) if (tpl[k].count() < f.minTemplatePulses) { TH1D *ht = (TH1D *) ftpl->Get(Form("hTemplate_%d", k)); if (!ht) continue; tpl[k] = Template(); for (int b = 0; b < 500; b++) if (ht->GetBinContent(b + 1) != 0) { tpl[k].n[b] = f.minTemplatePulses; tpl[k].sum[b] = f.minTemplatePulses * ht->GetBinContent(b + 1); } fromFile[k] = tpl[k].count() > 0; }
          ftpl->Close(); fout->cd();
        }
      for (int k = 0; k < 32; k++) if (tpl[k].count() < f.minTemplatePulses) { const int kb = (k / 4) * 4 + 2; if (tpl[kb].count() >= f.minTemplatePulses) { tpl[k] = tpl[kb]; borrowed[k] = true; } }
      for (int k = 0; k < 32; k++) { hasTpl[k] = tpl[k].count() >= f.minTemplatePulses; TH1D *ht = new TH1D(Form("hTemplate_%d", k), Form("ch %d template (%.0f pulses);t - t_{cross} [ns];normalised pulse", chan[k], tpl[k].count()), 500, -100, 400); for (int b = 0; b < 500; b++) if (tpl[k].n[b]) ht->SetBinContent(b + 1, tpl[k].sum[b] / tpl[k].n[b]); }
      printf("templates (b = borrowed from the beam-row cell, f = from %s): ", f.templateFile); for (int k = 0; k < 32; k++) printf("%d:%.0f%s ", chan[k], tpl[k].count(), fromFile[k] ? "f" : borrowed[k] ? "b" : ""); printf("\n");
      hRej->Reset();
    }
  const Long64_t pstop = pass == 0 ? std::min(stop, start + f.nTemplate) : stop;
  for (Long64_t i = start; i < pstop; i++)
    {
      h->GetEntry(i);
      if (f.requireComplete && !complete) continue;
      bool bad = false;
      for (int s = 0; s < 20 && !bad; s++) for (int l = 0; l < 8; l++) { if (!present[s][l]) continue; if (f.sfilLinksOnly && (l == 2 || l == 3 || l == 6 || l == 7)) continue; if ((f.rejectCRC && !crc_ok[s][l]) || (f.rejectHamming && ham[s][l])) { bad = true; break; } }
      if (bad || (f.rejectBackward && backward) || (f.rejectOverflow && overflow)) continue;
      pad->GetEntry(i); int npad = 0; padid = 0;
      for (int k = 0; k < 4; k++) { double ped = 0; for (int s = f.pedFirst; s <= f.pedLast; s++) ped += P[k][s]; ped /= (f.pedLast - f.pedFirst + 1); int mn = 1023; for (int s = 8; s <= 16; s++) if (P[k][s] >= 0 && P[k][s] < mn) mn = P[k][s]; if (ped - mn > f.paddleDip) { npad++; padid = 3 + k; } }
      if (f.requirePaddle == -1 && npad == 0) continue; if (f.requirePaddle >= 3 && padid != f.requirePaddle) continue; if (f.rejectMultiPaddle && npad != 1) continue;
      for (int L = 0; L < 4; L++) sf[L]->GetEntry(i);
      trig = i; if (pass == 1) kept++;
      double sq[4] = {0, 0, 0, 0}; bool sqok[4] = {false, false, false, false};
      for (int k = 0; k < 32; k++)
        {
          const int L = k / 8, e = (k / 4) % 2, c = k % 4; const Short_t *w = W[L][e][c];
          ok[k] = false; amp[k] = afit[k] = t0[k] = chi2[k] = pedv[k] = taur[k] = tauf[k] = pw[k] = -1;
          if (masked[k] || w[0] < 0) continue;
          if (pass == 1) hRej->Fill(0);
          double ped = 0; for (int s = f.pedFirst; s <= f.pedLast; s++) ped += w[s]; ped /= (f.pedLast - f.pedFirst + 1); pedv[k] = ped;
          double a = -1e9; int ipk = -1; for (int s = 4; s <= 12; s++) if (w[s] - ped > a) { a = w[s] - ped; ipk = s; } amp[k] = a;
          if (a <= f.minAmp) { if (pass == 1) hRej->Fill(1); continue; }
          // peak sample: the first sample at the maximum (a saturated plateau starts there)
          if (ipk < f.peakLo || ipk > f.peakHi) { if (pass == 1) hRej->Fill(2); continue; }
          double y[20]; int nsat = 0; for (int s = 0; s < 20; s++) { y[s] = w[s] - ped; if (w[s] >= 1023) nsat++; }
          if (nsat && !f.fitSaturated) { hRej->Fill(3); continue; }
          if (pass == 0)   // build the template from clean, unsaturated pulses
            {
              if (nsat || a < 100 || a > 1000) continue; const double tx = crossing(y, a, ipk); if (tx < 0) continue;
              for (int s = 0; s < 20; s++) { if (w[s] < 0) continue; const double u = 25. * s - tx; const int b = int(std::floor(u + 100.)); if (b >= 0 && b < 500) { tpl[k].sum[b] += y[s] / a; tpl[k].n[b] += 1; } }
              continue;
            }
          if (!hasTpl[k]) { hRej->Fill(3); continue; }
          const double tx0 = crossing(y, a, ipk); if (tx0 < 0) { hRej->Fill(3); continue; }
          double best = 1e30, bestA = 0, bestT = tx0; int ndf = 0;
          for (double tx = tx0 - 12.; tx <= tx0 + 12.; tx += 0.25)
            {
              double sym = 0, smm = 0; int n = 0;
              for (int s = 0; s < 20; s++) { if (w[s] < 0 || w[s] >= 1023) continue; const double m = tpl[k].eval(25. * s - tx), sg2 = f.noise * f.noise + std::pow(f.shapeErr * a * m, 2); sym += y[s] * m / sg2; smm += m * m / sg2; n++; }
              if (smm <= 0) continue; const double A = sym / smm; double c2 = 0;
              for (int s = 0; s < 20; s++) { if (w[s] < 0 || w[s] >= 1023) continue; const double m = tpl[k].eval(25. * s - tx), sg2 = f.noise * f.noise + std::pow(f.shapeErr * a * m, 2); c2 += std::pow(y[s] - A * m, 2) / sg2; }
              if (c2 < best) { best = c2; bestA = A; bestT = tx; ndf = n - 2; }
            }
          if (ndf <= 0 || bestA <= 0) { hRej->Fill(3); continue; }
          const double c2 = best / ndf; chi2[k] = c2; afit[k] = bestA; t0[k] = bestT;
          hChi[k]->Fill(c2);
          if (c2 > f.maxChi2) { hRej->Fill(4); continue; }
          ok[k] = true; hRej->Fill(5); hA[k]->Fill(afit[k]); hT0[k]->Fill(t0[k]); hAA[k]->Fill(a, afit[k]);
        }
      if (pass == 1) for (int L = 0; L < 4; L++) { const int ks = 8 * L + 2, kn = 8 * L + 4 + 2; if (ok[ks] && ok[kn]) { hSN[L]->Fill(afit[ks], afit[kn]); hAsym->Fill(L + 1, (afit[kn] - afit[ks]) / (afit[kn] + afit[ks])); sq[L] = std::sqrt(afit[ks] * afit[kn]); sqok[L] = true; } }
      if (pass == 1) for (int L = 0; L < 3; L++) if (sqok[L] && sqok[L + 1]) hLL[L]->Fill(sq[L], sq[L + 1]);
      if (sqok[1] && sqok[2] && sqok[3]) hSum->Fill(sq[0] + sq[1] + sq[2] + sq[3]);
      if (pass == 1) T->Fill();
    }
  }
  printf("triggers %lld-%lld: %lld kept\n", start, stop - 1, kept); for (int i = 0; i < 6; i++) printf("  %-20s %10.0f\n", rn[i], hRej->GetBinContent(i + 1));
  // figures
  gStyle->SetOptStat(0); gStyle->SetOptTitle(0); gStyle->SetPalette(57); gStyle->SetPadLeftMargin(0.16); gStyle->SetPadBottomMargin(0.15); gStyle->SetPadRightMargin(0.04); gStyle->SetPadTopMargin(0.1); gStyle->SetTitleSize(0.06, "xy"); gStyle->SetLabelSize(0.05, "xy"); gStyle->SetTitleOffset(1.2, "y"); gStyle->SetNdivisions(505, "xy");
  auto grid8x4 = [&](const char *name, auto draw) {
    TCanvas *c = new TCanvas(name, "", 3000, 1400); TPad *band = new TPad("band", "", 0, 0.95, 1, 1); band->Draw(); TPad *grid = new TPad("grid", "", 0, 0, 1, 0.95); grid->Draw(); grid->Divide(8, 4, 0.002, 0.004);
    for (int L = 0; L < 4; L++) for (int e = 0; e < 2; e++) for (int cell = 0; cell < 4; cell++) { const int k = 8 * L + 4 * e + cell; grid->cd(cell * 8 + 2 * L + (e ? 0 : 1) + 1); if (masked[k]) { TLatex t; t.SetNDC(); t.SetTextSize(0.1); t.SetTextColor(kGray + 1); t.DrawLatex(0.3, 0.5, Form("ch %d masked", chan[k])); continue; } draw(k); TLatex t; t.SetNDC(); t.SetTextSize(0.085); t.DrawLatex(0.18, 0.92, Form("SFIL%d %s cell %d  (ch %d)", L + 1, e ? "N" : "S", cell + 1, chan[k])); }
    return std::make_pair(c, band); };
  TString runTag = gSystem->BaseName(file); runTag.ReplaceAll("kcu0_", ""); runTag.ReplaceAll(".root", "");
  auto caption = [&](TPad *band, const char *what) { band->cd(); TLatex t; t.SetNDC(); t.SetTextSize(0.3); t.SetTextAlign(12); t.DrawLatex(0.003, 0.7, Form("run %s, triggers %lld-%lld, %lld kept; %s", runTag.Data(), start, stop - 1, kept, what)); t.DrawLatex(0.003, 0.25, Form("cuts: amplitude > %.0f ADC, peak sample %d-%d, template fit chi2/ndf < %.1f (noise %.0f ADC, shape error %.0f%%), saturated samples masked; columns per layer: North | South; rows: cells 1-4", f.minAmp, f.peakLo, f.peakHi, f.maxChi2, f.noise, 100 * f.shapeErr)); };
  { auto cb = grid8x4("cTpl", [&](int k) { TH1D *ht = (TH1D *) fout->Get(Form("hTemplate_%d", k)); if (ht) { ht->GetXaxis()->SetRangeUser(-60, 400); ht->SetLineWidth(2); ht->Draw("hist"); } }); caption(cb.second, "pulse templates: mean normalised pulse vs time from the 20% crossing (clean unsaturated pulses, 100-1000 ADC)"); cb.first->SaveAs("kcu0_fit_templates.png"); }
  { auto cb = grid8x4("cA", [&](int k) { gPad->SetLogy(); hA[k]->Draw(); }); caption(cb.second, "fitted amplitude A_fit"); cb.first->SaveAs("kcu0_fit_amplitude.png"); }
  { auto cb = grid8x4("cAA", [&](int k) { gPad->SetLogz(); hAA[k]->Draw("col"); }); caption(cb.second, "fitted amplitude vs max - pedestal (the fit extrapolates through saturation)"); cb.first->SaveAs("kcu0_fit_afit_vs_amp.png"); }
  { auto cb = grid8x4("cChi", [&](int k) { gPad->SetLogy(); hChi[k]->Draw(); }); caption(cb.second, "chi2/ndf of the pulse fit (before the chi2 cut)"); cb.first->SaveAs("kcu0_fit_chi2.png"); }
  { auto cb = grid8x4("cT0", [&](int k) { hT0[k]->Draw(); }); caption(cb.second, "fitted start time t0"); cb.first->SaveAs("kcu0_fit_t0.png"); }
  {
    TCanvas *c = new TCanvas("cCorr", "", 2400, 1200); c->Divide(4, 2, 0.003, 0.003);
    for (int L = 0; L < 4; L++) { c->cd(L + 1); gPad->SetLogz(); hSN[L]->Draw("col"); TLatex t; t.SetNDC(); t.SetTextSize(0.055); t.DrawLatex(0.2, 0.92, Form("SFIL%d beam row: N vs S", L + 1)); }
    for (int L = 0; L < 3; L++) { c->cd(5 + L); gPad->SetLogz(); hLL[L]->Draw("col"); TLatex t; t.SetNDC(); t.SetTextSize(0.055); t.DrawLatex(0.2, 0.92, Form("#sqrt{S#timesN}: SFIL%d vs SFIL%d", L + 2, L + 1)); }
    c->cd(8); gPad->SetLogz(); hAsym->Draw("col"); TLatex t; t.SetNDC(); t.SetTextSize(0.055); t.DrawLatex(0.2, 0.92, "(N - S)/(N + S) per layer");
    c->SaveAs("kcu0_fit_correlations.png");
  }
  fout->cd(); fout->Write(); printf("wrote %s and kcu0_fit_{amplitude,afit_vs_amp,chi2,t0,correlations}.png\n", out);
}
