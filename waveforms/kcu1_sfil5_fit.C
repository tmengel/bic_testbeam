// kcu1_sfil5_fit.C: template pulse fits of the 64 SFIL5 South single-SiPM channels of KCU1 (run 3062), the KCU1 twin of kcu0_sfil_fit.C.
//
//   root -l -b -q 'kcu1_sfil5_fit.C+(20000, 0)'
//
// Differences from the KCU0 macro: no paddle tree on KCU1 (paddle tags come from the KCU0 join through ts0); the pulses peak one
// sample earlier (samples 6-7), so the peak window is 5-8; the single-SiPM channels saturate at the 10-bit ceiling like the summed ones
// (ceiling - pedestal ~ 943 ADC), so saturated samples are masked and the template extrapolates; CRC / Hamming checks only on the links
// of the live ASIC2 (links 4, 5); noise 2 ADC. Templates: per channel from clean pulses of 100-850 ADC aligned on the 20 % crossing.
// Output tree "fits": per trigger trig, ts0, and per channel index k = 16 array + i (arrays A B C D): chan, ok, ped, amp (raw peak), afit, t0, chi2, nsat.
#include "TFile.h"
#include "TTree.h"
#include "TH1D.h"
#include "TH2F.h"
#include "TString.h"
#include <set>
#include <vector>
#include <algorithm>
#include <cmath>
#include <cstdio>

struct Fit1Flags
{
  bool requireComplete = true, rejectCRC = true, rejectHamming = true, rejectBackward = true, rejectOverflow = true;
  std::set<int> links = {4, 5};   // ASIC2 = SFIL5 South
  int pedFirst = 1, pedLast = 3; double minAmp = 60.; int peakLo = 5, peakHi = 8; double maxChi2 = 8., noise = 2., shapeErr = 0.05;
  double tplLo = 100., tplHi = 850.; Long64_t nTemplate = 150000; int minTemplatePulses = 100;
  const char *templateFile = "";   // fit output of another run: channels with too few clean pulses take its hTemplate_k (low-rate run 3065)
};
static double crossing1(const double *y, double a, int ipk) { for (int s = 1; s <= ipk; s++) if (y[s] >= 0.2 * a) return 25. * (s - 1 + (0.2 * a - y[s - 1]) / std::max(1e-6, y[s] - y[s - 1])); return -1; }
struct Template1 { std::vector<double> sum, n; Template1() : sum(500, 0), n(500, 0) {} double eval(double u) const { const double b = u + 100.; if (b < 0 || b >= 499) return 0; const int i = int(b); const double f = b - i; const double v0 = n[i] ? sum[i] / n[i] : 0, v1 = n[i + 1] ? sum[i + 1] / n[i + 1] : 0; return v0 * (1 - f) + v1 * f; } double count() const { double m = 0; for (int i = 150; i < 250; i++) m = std::max(m, n[i]); return m; } };

void kcu1_sfil5_fit(Long64_t nmax = 20000, Long64_t start = 0, Fit1Flags f = Fit1Flags(), const char *file = "../data/kcu1_3062.root", const char *out = "kcu1_sfil5_fit.root")
{
  TFile *fin = TFile::Open(file); if (!fin || fin->IsZombie()) { printf("cannot open %s\n", file); return; }
  TTree *h = (TTree *) fin->Get("header"), *t5 = (TTree *) fin->Get("sfil5");
  Bool_t complete, backward, overflow, crc_ok[20][8], present[20][8]; UChar_t ham[20][8]; Long64_t ts0;
  h->SetBranchStatus("*", 0); for (const char *b : {"complete", "backward", "nsamp_overflow", "crc_ok", "ham", "present", "ts0"}) h->SetBranchStatus(b, 1);
  h->SetBranchAddress("complete", &complete); h->SetBranchAddress("backward", &backward); h->SetBranchAddress("nsamp_overflow", &overflow); h->SetBranchAddress("crc_ok", crc_ok); h->SetBranchAddress("ham", ham); h->SetBranchAddress("present", present); h->SetBranchAddress("ts0", &ts0);
  Int_t CH[4][16]; Short_t W[4][16][20], TT[4][16][20], TA[4][16][20]; const char *an[4] = {"A", "B", "C", "D"};
  for (int a = 0; a < 4; a++) { t5->SetBranchAddress(Form("%s_ch", an[a]), CH[a]); t5->SetBranchAddress(Form("%s_adc", an[a]), W[a]); t5->SetBranchAddress(Form("%s_tot", an[a]), TT[a]); t5->SetBranchAddress(Form("%s_toa", an[a]), TA[a]); }
  TFile *fout = new TFile(out, "RECREATE");
  Long64_t trig; Int_t chan[64], nsat[64], tot[64], tots[64], toa[64], toas[64]; Bool_t ok[64]; Float_t amp[64], afit[64], t0[64], chi2[64], pedv[64];
  TTree *T = new TTree("fits", "KCU1 SFIL5 South pulse fits per trigger; k = 16 array + i (A B C D)");
  T->Branch("trig", &trig, "trig/L"); T->Branch("ts0", &ts0, "ts0/L"); T->Branch("chan", chan, "chan[64]/I"); T->Branch("ok", ok, "ok[64]/O"); T->Branch("ped", pedv, "ped[64]/F"); T->Branch("amp", amp, "amp[64]/F"); T->Branch("afit", afit, "afit[64]/F"); T->Branch("t0", t0, "t0[64]/F"); T->Branch("chi2", chi2, "chi2[64]/F"); T->Branch("nsat", nsat, "nsat[64]/I"); T->Branch("tot", tot, "tot[64]/I"); T->Branch("tots", tots, "tots[64]/I"); T->Branch("toa", toa, "toa[64]/I"); T->Branch("toas", toas, "toas[64]/I");   // decoded TOT (max over samples) and its sample, first TOA code and its sample; -1 = none
  TH1D *hRej = new TH1D("hReasons", "waveforms per outcome;;waveforms", 6, 0, 6); const char *rn[6] = {"looked at", "amp < minAmp", "peak outside window", "fit failed", "chi2 too large", "kept"}; for (int i = 0; i < 6; i++) hRej->GetXaxis()->SetBinLabel(i + 1, rn[i]);
  TH1D *hPk = new TH1D("hPeakSample", "peak sample of pulses > minAmp;sample;pulses", 20, -0.5, 19.5);
  TH2F *hAA[64]; TH1D *hA[64], *hChi[64];
  t5->GetEntry(0); for (int a = 0; a < 4; a++) for (int i = 0; i < 16; i++) { const int k = 16 * a + i; chan[k] = CH[a][i]; hAA[k] = new TH2F(Form("hAfitVsAmp_%d", k), Form("%s%d ch %d;max - pedestal [ADC];A_{fit} [ADC]", an[a], i, chan[k]), 100, 0, 1100, 100, 0, 3000); hA[k] = new TH1D(Form("hAfit_%d", k), Form("%s%d ch %d;A_{fit} [ADC];pulses", an[a], i, chan[k]), 150, 0, 3000); hChi[k] = new TH1D(Form("hChi2_%d", k), Form("%s%d ch %d;#chi^{2}/ndf;pulses", an[a], i, chan[k]), 100, 0, 20); }
  const Long64_t total = h->GetEntries(); if (start < 0) start = 0; const Long64_t stop = nmax < 0 ? total : std::min(start + nmax, total);
  Long64_t kept = 0; Template1 tpl[64]; bool hasTpl[64] = {false};
  for (int pass = 0; pass < 2; pass++)
  {
  if (pass == 1)
    {
      bool fromFile[64] = {false};
      if (f.templateFile && *f.templateFile)
        {
          TFile *ftpl = TFile::Open(f.templateFile); if (!ftpl || ftpl->IsZombie()) { printf("cannot open template file %s\n", f.templateFile); return; }
          for (int k = 0; k < 64; k++) if (tpl[k].count() < f.minTemplatePulses) { TH1D *ht = (TH1D *) ftpl->Get(Form("hTemplate_%d", k)); if (!ht) continue; tpl[k] = Template1(); for (int b = 0; b < 500; b++) if (ht->GetBinContent(b + 1) != 0) { tpl[k].n[b] = f.minTemplatePulses; tpl[k].sum[b] = f.minTemplatePulses * ht->GetBinContent(b + 1); } fromFile[k] = tpl[k].count() > 0; }
          ftpl->Close(); fout->cd();
        }
      printf("templates (f = from %s): ", f.templateFile); for (int k = 0; k < 64; k++) { hasTpl[k] = tpl[k].count() >= f.minTemplatePulses; TH1D *ht = new TH1D(Form("hTemplate_%d", k), Form("ch %d template (%.0f pulses);t - t_{cross} [ns];normalised pulse", chan[k], tpl[k].count()), 500, -100, 400); for (int b = 0; b < 500; b++) if (tpl[k].n[b]) ht->SetBinContent(b + 1, tpl[k].sum[b] / tpl[k].n[b]); if (hasTpl[k]) printf("%d:%.0f%s ", chan[k], tpl[k].count(), fromFile[k] ? "f" : ""); } printf("\n");
      hRej->Reset();
    }
  const Long64_t pstop = pass == 0 ? std::min(stop, start + f.nTemplate) : stop;
  for (Long64_t i = start; i < pstop; i++)
    {
      h->GetEntry(i);
      if (f.requireComplete && !complete) continue;
      bool bad = false;
      for (int s = 0; s < 20 && !bad; s++) for (int l : f.links) { if (!present[s][l]) continue; if ((f.rejectCRC && !crc_ok[s][l]) || (f.rejectHamming && ham[s][l])) bad = true; }
      if (bad || (f.rejectBackward && backward) || (f.rejectOverflow && overflow)) continue;
      t5->GetEntry(i); trig = i; if (pass == 1) kept++;
      for (int k = 0; k < 64; k++)
        {
          const Short_t *w = W[k / 16][k % 16];
          ok[k] = false; amp[k] = afit[k] = t0[k] = chi2[k] = pedv[k] = -1; nsat[k] = 0; tot[k] = tots[k] = toa[k] = toas[k] = -1;
          if (w[0] < 0) continue;
          { const Short_t *tw = TT[k / 16][k % 16], *ta = TA[k / 16][k % 16]; for (int s = 0; s < 20; s++) { if (tw[s] > 0) { const int v = tw[s] < 512 ? tw[s] : 8 * (tw[s] - 512) + 4; if (v > tot[k]) { tot[k] = v; tots[k] = s; } } if (ta[s] > 0 && toa[k] < 0) { toa[k] = ta[s]; toas[k] = s; } } }
          if (pass == 1) hRej->Fill(0);
          double ped = 0; for (int s = f.pedFirst; s <= f.pedLast; s++) ped += w[s]; ped /= (f.pedLast - f.pedFirst + 1); pedv[k] = ped;
          double a = -1e9; int ipk = -1; for (int s = 4; s <= 12; s++) if (w[s] - ped > a) { a = w[s] - ped; ipk = s; } amp[k] = a;
          double y[20]; for (int s = 0; s < 20; s++) { y[s] = w[s] - ped; if (w[s] >= 1023) nsat[k]++; }
          if (a <= f.minAmp) { if (pass == 1) hRej->Fill(1); continue; }
          if (pass == 1) hPk->Fill(ipk);
          if (ipk < f.peakLo || ipk > f.peakHi) { if (pass == 1) hRej->Fill(2); continue; }
          if (pass == 0)
            {
              if (nsat[k] || a < f.tplLo || a > f.tplHi) continue; const double tx = crossing1(y, a, ipk); if (tx < 0) continue;
              for (int s = 0; s < 20; s++) { if (w[s] < 0) continue; const double u = 25. * s - tx; const int b = int(std::floor(u + 100.)); if (b >= 0 && b < 500) { tpl[k].sum[b] += y[s] / a; tpl[k].n[b] += 1; } }
              continue;
            }
          if (!hasTpl[k]) { hRej->Fill(3); continue; }
          const double tx0 = crossing1(y, a, ipk); if (tx0 < 0) { hRej->Fill(3); continue; }
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
          const double c2 = best / ndf; chi2[k] = c2; afit[k] = bestA; t0[k] = bestT; hChi[k]->Fill(c2);
          if (c2 > f.maxChi2) { hRej->Fill(4); continue; }
          ok[k] = true; hRej->Fill(5); hA[k]->Fill(afit[k]); hAA[k]->Fill(a, afit[k]);
        }
      if (pass == 1) T->Fill();
    }
  }
  printf("triggers %lld-%lld: %lld kept\n", start, stop - 1, kept); for (int i = 0; i < 6; i++) printf("  %-20s %10.0f\n", rn[i], hRej->GetBinContent(i + 1));
  printf("peak sample of pulses > minAmp:"); for (int s = 3; s <= 13; s++) printf(" %d:%.0f", s, hPk->GetBinContent(s + 1)); printf("\n");
  fout->cd(); fout->Write(); printf("wrote %s\n", out);
}
