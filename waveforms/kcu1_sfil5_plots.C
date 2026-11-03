// kcu1_sfil5_plots.C: the figures of kcu0_sfil_fit.C for KCU1 / SFIL5 South, drawn from the output of kcu1_sfil5_fit.C (tree 'fits' and the
// per-channel histograms) plus a waveform gallery from the raw tree.
//
//   root -l -b -q 'kcu1_sfil5_plots.C+("kcu1_sfil5_fit.root", "../data/kcu1_3062.root", 100000)'
//
// Panels are 4 rows (arrays A B C D, bottom to top in the stack: A < B < D < C in height) x 16 columns (SiPM index in the array):
//   kcu1_fit_templates.png     mean normalised pulse per channel        kcu1_fit_amplitude.png   fitted amplitude
//   kcu1_fit_afit_vs_amp.png   fitted amplitude vs raw peak             kcu1_fit_chi2.png        chi2/ndf
//   kcu1_fit_tot_vs_afit.png   decoded TOT vs fitted amplitude (the calibration overlap for saturated pulses)
//   kcu1_fit_sums.png          array sums B, C, D and the B+C+D total, fitted-where-ok vs raw peaks; fraction of saturated channels per trigger
//   kcu1_wf_pedsub.png         pedestal-subtracted waveforms (sample x ADC, log z) of the first nraw clean triggers, same grid
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

void kcu1_sfil5_plots(const char *fitfile = "kcu1_sfil5_fit.root", const char *rawfile = "../data/kcu1_3062.root", Long64_t nraw = 100000)
{
  gStyle->SetOptStat(0); gStyle->SetOptTitle(0); gStyle->SetPalette(57); gStyle->SetPadLeftMargin(0.16); gStyle->SetPadBottomMargin(0.15); gStyle->SetPadRightMargin(0.04); gStyle->SetPadTopMargin(0.1); gStyle->SetTitleSize(0.065, "xy"); gStyle->SetLabelSize(0.055, "xy"); gStyle->SetNdivisions(505, "xy"); gStyle->SetTitleOffset(1.1, "y");
  TFile *ff = TFile::Open(fitfile); if (!ff || ff->IsZombie()) { printf("cannot open %s\n", fitfile); return; }
  TTree *T = (TTree *) ff->Get("fits"); Int_t chan[64], nsat[64], tot[64]; Bool_t ok[64]; Float_t amp[64], afit[64], chi2[64];
  T->SetBranchAddress("chan", chan); T->SetBranchAddress("nsat", nsat); T->SetBranchAddress("tot", tot); T->SetBranchAddress("ok", ok); T->SetBranchAddress("amp", amp); T->SetBranchAddress("afit", afit); T->SetBranchAddress("chi2", chi2);
  const char *an[4] = {"A", "B", "C", "D"}; const int rowOf[4] = {3, 2, 0, 1};   // draw C on top, then D, B, A (the vertical order in the stack)
  T->GetEntry(0);
  TH2F *hTA[64]; TH1D *hSum[4], *hSumRaw[4], *hNsat = new TH1D("hNsatTrig", ";saturated SFIL5 South channels per trigger;triggers", 20, -0.5, 19.5);
  for (int k = 0; k < 64; k++) hTA[k] = new TH2F(Form("hTotAfit_%d", k), Form("%s%d ch %d;A_{fit} [ADC];decoded TOT", an[k / 16], k % 16, chan[k]), 60, 0, 3000, 60, 0, 3000);
  const char *sn[4] = {"B", "C", "D", "B+C+D"}; for (int i = 0; i < 4; i++) { hSum[i] = new TH1D(Form("hSum_%d", i), Form("%s;sum of fitted amplitudes [ADC];triggers", sn[i]), 120, 0, 24000); hSumRaw[i] = new TH1D(Form("hSumRaw_%d", i), Form("%s raw", sn[i]), 120, 0, 24000); }
  const Long64_t n = T->GetEntries();
  for (Long64_t i = 0; i < n; i++)
    {
      T->GetEntry(i); double s[3] = {0, 0, 0}, r[3] = {0, 0, 0}; int ns = 0;
      for (int k = 0; k < 64; k++)
        {
          if (amp[k] < 0) continue;
          if (ok[k] && tot[k] > 0) hTA[k]->Fill(afit[k], tot[k]);
          const int a = k / 16; if (a >= 1) { const double v = ok[k] ? afit[k] : std::max(0.f, amp[k]); s[a - 1] += v; r[a - 1] += std::max(0.f, amp[k]); }
          if (nsat[k] > 0) ns++;
        }
      for (int a = 0; a < 3; a++) { hSum[a]->Fill(s[a]); hSumRaw[a]->Fill(r[a]); } hSum[3]->Fill(s[0] + s[1] + s[2]); hSumRaw[3]->Fill(r[0] + r[1] + r[2]); hNsat->Fill(ns);
    }
  auto grid = [&](const char *name, auto draw) {
    TCanvas *c = new TCanvas(name, "", 4000, 1200); TPad *band = new TPad("band", "", 0, 0.95, 1, 1); band->Draw(); TPad *g = new TPad("grid", "", 0, 0, 1, 0.95); g->Draw(); g->Divide(16, 4, 0.0015, 0.004);
    for (int a = 0; a < 4; a++) for (int i = 0; i < 16; i++) { const int k = 16 * a + i; g->cd(rowOf[a] * 16 + i + 1); draw(k); TLatex t; t.SetNDC(); t.SetTextSize(0.09); t.DrawLatex(0.2, 0.91, Form("%s%d ch %d", an[a], i, chan[k])); }
    return std::make_pair(c, band); };
  auto caption = [&](TPad *band, const char *what) { band->cd(); TLatex t; t.SetNDC(); t.SetTextSize(0.32); t.SetTextAlign(12); t.DrawLatex(0.003, 0.65, Form("run 3062, KCU1 SFIL5 South single SiPMs, %lld triggers kept by the fit cuts; %s", n, what)); t.DrawLatex(0.003, 0.2, "rows top to bottom: arrays C, D, B, A (their order in height); columns: SiPM index in the array; cuts: complete, no CRC/Hamming on ASIC2, raw peak > 60 ADC, peak sample 5-8, chi2/ndf < 8, saturated samples masked"); };
  { auto cb = grid("cTpl", [&](int k) { TH1D *h = (TH1D *) ff->Get(Form("hTemplate_%d", k)); if (h && h->GetEntries()) { h->GetXaxis()->SetRangeUser(-60, 400); h->SetLineWidth(2); h->Draw("hist"); } }); caption(cb.second, "pulse templates: mean normalised pulse vs time from the 20 % crossing"); cb.first->SaveAs("kcu1_fit_templates.png"); }
  { auto cb = grid("cA", [&](int k) { TH1D *h = (TH1D *) ff->Get(Form("hAfit_%d", k)); if (h) { gPad->SetLogy(); h->Draw(); } }); caption(cb.second, "fitted amplitude A_fit"); cb.first->SaveAs("kcu1_fit_amplitude.png"); }
  { auto cb = grid("cAA", [&](int k) { TH2F *h = (TH2F *) ff->Get(Form("hAfitVsAmp_%d", k)); if (h) { gPad->SetLogz(); h->Draw("col"); } }); caption(cb.second, "fitted amplitude vs raw peak (max - pedestal): the saturated pulses sit at the right edge; the template fit does NOT recover them on KCU1 (shape changes)"); cb.first->SaveAs("kcu1_fit_afit_vs_amp.png"); }
  { auto cb = grid("cChi", [&](int k) { TH1D *h = (TH1D *) ff->Get(Form("hChi2_%d", k)); if (h) { gPad->SetLogy(); h->Draw(); } }); caption(cb.second, "chi2/ndf of the pulse fit (before the chi2 cut)"); cb.first->SaveAs("kcu1_fit_chi2.png"); }
  { auto cb = grid("cTA", [&](int k) { gPad->SetLogz(); hTA[k]->Draw("col"); }); caption(cb.second, "decoded TOT vs fitted amplitude for pulses that carry a TOT word: the overlap below the ceiling calibrates TOT as the amplitude of saturated pulses"); cb.first->SaveAs("kcu1_fit_tot_vs_afit.png"); }
  {
    TCanvas *c = new TCanvas("cSums", "", 2400, 1200); c->Divide(3, 2, 0.003, 0.003);
    for (int i = 0; i < 4; i++) { c->cd(i + 1); gPad->SetLogy(); hSum[i]->SetLineWidth(2); hSum[i]->Draw(); hSumRaw[i]->SetLineColor(2); hSumRaw[i]->SetLineWidth(2); hSumRaw[i]->Draw("same"); TLatex t; t.SetNDC(); t.SetTextSize(0.05); t.DrawLatex(0.4, 0.85, Form("array %s: fit-where-ok (black), raw peaks (red)", sn[i])); t.DrawLatex(0.4, 0.78, Form("mean %.0f / %.0f", hSum[i]->GetMean(), hSumRaw[i]->GetMean())); }
    c->cd(5); gPad->SetLogy(); hNsat->SetLineWidth(2); hNsat->Draw(); TLatex t; t.SetNDC(); t.SetTextSize(0.05); t.DrawLatex(0.3, 0.85, Form("saturated channels per trigger: mean %.2f, zero in %.0f %%", hNsat->GetMean(), 100. * hNsat->GetBinContent(1) / std::max(1., hNsat->GetEntries())));
    c->SaveAs("kcu1_fit_sums.png");
  }
  // waveform gallery from the raw tree: pedestal-subtracted sample x ADC per channel, clean triggers only
  TFile *fr = TFile::Open(rawfile); if (!fr || fr->IsZombie()) { printf("no raw file, skipping the gallery\n"); return; }
  TTree *h = (TTree *) fr->Get("header"), *t5 = (TTree *) fr->Get("sfil5"); Bool_t complete, backward, overflow, crc_ok[20][8], present[20][8]; UChar_t ham[20][8];
  h->SetBranchStatus("*", 0); for (const char *b : {"complete", "backward", "nsamp_overflow", "crc_ok", "ham", "present"}) h->SetBranchStatus(b, 1);
  h->SetBranchAddress("complete", &complete); h->SetBranchAddress("backward", &backward); h->SetBranchAddress("nsamp_overflow", &overflow); h->SetBranchAddress("crc_ok", crc_ok); h->SetBranchAddress("ham", ham); h->SetBranchAddress("present", present);
  Short_t W[4][16][20]; for (int a = 0; a < 4; a++) t5->SetBranchAddress(Form("%s_adc", an[a]), W[a]);
  TH2F *hW[64]; for (int k = 0; k < 64; k++) hW[k] = new TH2F(Form("hWf_%d", k), ";sample;ADC - pedestal", 20, -0.5, 19.5, 120, -100, 1100);
  Long64_t used = 0; const Long64_t stop = std::min(nraw, h->GetEntries());
  for (Long64_t i = 0; i < stop; i++)
    {
      h->GetEntry(i); if (!complete || backward || overflow) continue; bool bad = false;
      for (int s = 0; s < 20 && !bad; s++) for (int l : {4, 5}) if (present[s][l] && (!crc_ok[s][l] || ham[s][l])) bad = true;
      if (bad) continue; t5->GetEntry(i); used++;
      for (int k = 0; k < 64; k++) { const Short_t *w = W[k / 16][k % 16]; if (w[0] < 0) continue; const double ped = (w[1] + w[2] + w[3]) / 3.; for (int s = 0; s < 20; s++) hW[k]->Fill(s, w[s] - ped); }
    }
  { auto cb = grid("cWf", [&](int k) { gPad->SetLogz(); hW[k]->Draw("col"); }); cb.second->cd(); TLatex t; t.SetNDC(); t.SetTextSize(0.32); t.SetTextAlign(12); t.DrawLatex(0.003, 0.65, Form("run 3062, KCU1 SFIL5 South: pedestal-subtracted waveforms of %lld clean triggers (complete, no CRC/Hamming on ASIC2, no backward/overflow), log z", used)); t.DrawLatex(0.003, 0.2, "pulses peak at samples 6-7 (one earlier than KCU0); the ceiling is at 1023 - pedestal ~ 943; rows top to bottom: arrays C, D, B, A"); cb.first->SaveAs("kcu1_wf_pedsub.png"); }
  printf("wrote kcu1_fit_{templates,amplitude,afit_vs_amp,chi2,tot_vs_afit,sums}.png and kcu1_wf_pedsub.png\n");
}
