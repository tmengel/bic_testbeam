// cosmic_mip.C: MIP response of the SFIL cells in the cosmic runs (KCU0 3088 / 3091): pulse finding after a per-sample common-mode
// subtraction, per-channel spectra on offline-selected muons, and a compact per-trigger tree for the Landau fits (cosmic_mip_fit.py).
//
//   root -l -b -q 'cosmic_mip.C+(file, out, start, nmax)'       nmax < 0 = to the end
//
// What the 3088 waveforms look like (30 k-trigger study, 2026-10-30): the real pulses sit at the BEAM timing (rise at sample 6,
// peak 7-9, the usual shape, 40-70 ADC for a muon through a cell), not at samples 12-15 as the first 2.6 k-trigger look had said:
// that "late" population is a slow baseline ramp whose rms grows from 7 ADC (sample 2) to 25 ADC (sample 19). The ramp is shared
// by the CONNECTED channels of an ASIC half (correlation 0.66-0.98, strongest within one summing board; the 'unused' channels of the
// half carry none of it, corr -0.1), so per sample the median over the connected channels of the half (the noisy SFIL1 North
// column left out) is subtracted from every channel of that half: rms 5-9 -> 1-4 ADC. Its bias: the median of a half of 6 to 10
// channels is pulled only when at least half of them carry signal (a vertical muon through one board lights 3-4 of them, so the
// 6-channel halves ASIC0/1 (SFIL3 S + SFIL4 S c1-2) and ASIC2/1 (SFIL3 N + SFIL4 N c1-2) are biased for vertical tracks; Ax, the
// same estimate with the channel's own board left out of the reference, is unbiased but noisier (the within-board correlation is
// lost). The clean muon tag is the coincidence of the two ENDS of one cell (different ASICs, independent noise): e.g. S3 c3 69 +
// N3 c3 68 ADC; pairs of neighbouring cells of one board at a constant ~27 ADC are within-board noise, not tracks.
// Channels k = 8 layer + 4 end + cell (SFIL1-4, end 0 South 1 North), 32-35 = SFIL5 North cells 1-4 (on KCU0 ch 153 / 156 / 163 /
// 167 since Aug 14 21:55), 36-39 = paddles P3-P6 (negative pulses: pdip = median of samples 1-5 minus the minimum of 6-11; a
// ~90 ADC dip in exactly one paddle accompanies half of the cosmic triggers, so the paddles are part of that trigger).
// Per channel and trigger: pedestal = mean of samples 1-4; A = max of the corrected samples 5-11, pk its sample, Q = sum of
// 6-12, Al = max of 13-19 (the ramp window, control), Ax / Au = A with the own-board-excluded / unused-channel reference, noise =
// rms of the corrected samples 1-4 (raw = before the subtraction).
// Spectra per channel: hA all triggers; hB = the other end of the same cell has A > 15 (both-ends tag); hV = the other two cells
// of 1-3 of the same board have A > 12 (vertical muon; cell 4 requires cells 2-3; uses Ax); hGM_L_c = sqrt(A_S A_N) of the cell when
// both ends > 10 (attenuation-free); hAl = the late window (control); hPk = peak sample for A > 15.
#include "TFile.h"
#include "TTree.h"
#include "TH1D.h"
#include "TString.h"
#include <algorithm>
#include <vector>
#include <cmath>
#include <cstdio>

void cosmic_mip(const char *file = "../../data/kcu0_3088.root", const char *out = "cosmic_mip.root", Long64_t start = 0, Long64_t nmax = -1)
{
  const int P0 = 1, P1 = 4;                 // pedestal samples
  const int A0 = 5, A1 = 11, Q0 = 6, Q1 = 12, L0 = 13, L1 = 19;   // pulse, charge and late (control) windows
  const double THR_B = 15., THR_V = 12., THR_GM = 10.;
  TFile *fin = TFile::Open(file); if (!fin || fin->IsZombie()) { printf("cannot open %s\n", file); return; }
  TTree *h = (TTree *) fin->Get("header"), *un = (TTree *) fin->Get("unused"), *pad = (TTree *) fin->Get("paddles"), *s5 = (TTree *) fin->Get("sfil5n"), *sf[4], *mp = (TTree *) fin->Get("map");
  for (int L = 0; L < 4; L++) sf[L] = (TTree *) fin->Get(Form("sfil%d", L + 1));
  const bool has5 = s5 != nullptr;
  int asic[288], half[288]; { Int_t ch, a, hf; mp->SetBranchAddress("ch", &ch); mp->SetBranchAddress("asic", &a); mp->SetBranchAddress("half", &hf); for (Long64_t i = 0; i < mp->GetEntries(); i++) { mp->GetEntry(i); asic[ch] = a; half[ch] = hf; } }
  Bool_t complete; Long64_t ts0; h->SetBranchStatus("*", 0); h->SetBranchStatus("complete", 1); h->SetBranchStatus("ts0", 1); h->SetBranchAddress("complete", &complete); h->SetBranchAddress("ts0", &ts0);
  const int NU = 176; Int_t uch[NU]; Short_t uadc[NU][20]; un->SetBranchAddress("ch", uch); un->SetBranchAddress("adc", uadc);
  const int NK = 40; Short_t W[NK][20]; Int_t CH[NK]; const char *pn[4] = {"P3", "P4", "P5", "P6"};
  for (int L = 0; L < 4; L++) for (int c = 0; c < 4; c++) { sf[L]->SetBranchAddress(Form("Sc%d_adc", c + 1), W[8 * L + c]); sf[L]->SetBranchAddress(Form("Nc%d_adc", c + 1), W[8 * L + 4 + c]); sf[L]->SetBranchAddress(Form("Sc%d_ch", c + 1), &CH[8 * L + c]); sf[L]->SetBranchAddress(Form("Nc%d_ch", c + 1), &CH[8 * L + 4 + c]); }
  if (has5) for (int c = 0; c < 4; c++) { s5->SetBranchAddress(Form("Nc%d_adc", c + 1), W[32 + c]); s5->SetBranchAddress(Form("Nc%d_ch", c + 1), &CH[32 + c]); }
  for (int k = 0; k < 4; k++) { pad->SetBranchAddress(Form("%s_adc", pn[k]), W[36 + k]); pad->SetBranchAddress(Form("%s_ch", pn[k]), &CH[36 + k]); }

  TFile *fout = new TFile(out, "RECREATE");
  Long64_t trig; Int_t padbits; Float_t A[NK], Q[NK], Al[NK], Ax[NK], Au[NK], noise[NK], pdip[4]; Short_t pk[NK]; Int_t chan[NK];
  TTree *T = new TTree("cos", "cosmic pulses per trigger; k = 8 layer + 4 end + cell (SFIL1-4), 32-35 SFIL5 North cells 1-4, 36-39 paddles P3-P6 (A = dip)");
  T->Branch("trig", &trig, "trig/L"); T->Branch("ts0", &ts0, "ts0/L"); T->Branch("pad", &padbits, "pad/I"); T->Branch("pdip", pdip, "pdip[4]/F"); T->Branch("chan", chan, Form("chan[%d]/I", NK));
  T->Branch("A", A, Form("A[%d]/F", NK)); T->Branch("pk", pk, Form("pk[%d]/S", NK)); T->Branch("Q", Q, Form("Q[%d]/F", NK)); T->Branch("Al", Al, Form("Al[%d]/F", NK)); T->Branch("Ax", Ax, Form("Ax[%d]/F", NK)); T->Branch("Au", Au, Form("Au[%d]/F", NK)); T->Branch("noise", noise, Form("noise[%d]/F", NK));
  // corrected waveforms (x 10, 0.1 ADC) of the 36 SFIL channels, stored for triggers with any A > 12 outside the SFIL1 North column (shape cuts and template fits offline)
  Int_t nwf = 0; Short_t wf[36 * 20]; T->Branch("nwf", &nwf, "nwf/I"); T->Branch("wf", wf, "wf[nwf]/S");
  TH1D *hA[NK], *hB[NK], *hV[NK], *hAl[NK], *hPk[NK], *hNoise[NK], *hRaw[NK], *hGM[4][4], *hCM[8], *hCMrms[8];
  for (int k = 0; k < NK; k++)
    {
      hA[k] = new TH1D(Form("hA_%d", k), Form("k %d all triggers;A [ADC];triggers", k), 250, -50, 200); hB[k] = new TH1D(Form("hB_%d", k), Form("k %d, other end of the cell > %.0f;A [ADC];triggers", k, THR_B), 250, -50, 200);
      hV[k] = new TH1D(Form("hV_%d", k), Form("k %d, vertical muon (other cells of the board > %.0f);A_{x} [ADC];triggers", k, THR_V), 250, -50, 200); hAl[k] = new TH1D(Form("hAl_%d", k), Form("k %d late window (control);A_{l} [ADC];triggers", k), 250, -50, 200);
      hPk[k] = new TH1D(Form("hPk_%d", k), Form("k %d peak sample (A > 15);sample;pulses", k), 20, -0.5, 19.5); hNoise[k] = new TH1D(Form("hNoise_%d", k), Form("k %d noise after the subtraction;rms of samples 1-4 [ADC];triggers", k), 100, 0, 20); hRaw[k] = new TH1D(Form("hRaw_%d", k), Form("k %d noise before;rms of samples 1-4 [ADC];triggers", k), 100, 0, 20);
    }
  for (int L = 0; L < 4; L++) for (int c = 0; c < 4; c++) hGM[L][c] = new TH1D(Form("hGM_%d_%d", L + 1, c + 1), Form("SFIL%d cell %d #sqrt{A_{S}A_{N}}, both ends > %.0f;#sqrt{A_{S}A_{N}} [ADC];triggers", L + 1, c + 1, THR_GM), 200, 0, 200);
  for (int g = 0; g < 8; g++) { hCM[g] = new TH1D(Form("hCM_%d", g), Form("ASIC %d half %d common mode;median(adc - ped) of the connected channels [ADC];samples", g / 2, g % 2), 200, -50, 50); hCMrms[g] = new TH1D(Form("hCMrms_%d", g), Form("ASIC %d half %d;rms over samples of the common mode [ADC];triggers", g / 2, g % 2), 100, 0, 30); }
  TH1D *hReasons = new TH1D("hReasons", "triggers;;triggers", 3, 0, 3); hReasons->GetXaxis()->SetBinLabel(1, "looked at"); hReasons->GetXaxis()->SetBinLabel(2, "complete"); hReasons->GetXaxis()->SetBinLabel(3, "written");

  const Long64_t total = h->GetEntries(); const Long64_t stop = nmax < 0 ? total : std::min(total, start + nmax);
  std::vector<double> buf; double cm[8][20], cmu[8][20], ped[NK]; int grp[NK], brd[NK]; bool valid[NK];
  auto isRef = [&](int j) { return j < 36 && valid[j] && !(CH[j] == 146 || CH[j] == 150 || CH[j] == 154 || CH[j] == 158); };
  auto median = [&](int g, int exclBoard, int s) { buf.clear(); for (int j = 0; j < 36; j++) if (isRef(j) && grp[j] == g && brd[j] != exclBoard) buf.push_back(W[j][s] - ped[j]); if (buf.size() < 2) return 0.; std::nth_element(buf.begin(), buf.begin() + buf.size() / 2, buf.end()); double m = buf[buf.size() / 2]; if (buf.size() % 2 == 0) { double lo = *std::max_element(buf.begin(), buf.begin() + buf.size() / 2); m = 0.5 * (m + lo); } return m; };
  auto wmax = [](const double *y, int a, int b, int *at = nullptr) { double m = -1e9; for (int s = a; s <= b; s++) if (y[s] > m) { m = y[s]; if (at) *at = s; } return m; };
  for (Long64_t i = start; i < stop; i++)
    {
      h->GetEntry(i); hReasons->Fill(0.5); if (!complete) continue; hReasons->Fill(1.5);
      un->GetEntry(i); pad->GetEntry(i); for (int L = 0; L < 4; L++) sf[L]->GetEntry(i); if (has5) s5->GetEntry(i);
      for (int k = 0; k < NK; k++)
        {
          chan[k] = CH[k]; A[k] = Q[k] = Al[k] = Ax[k] = Au[k] = noise[k] = -999; pk[k] = -1; valid[k] = false; ped[k] = 0; grp[k] = -1; brd[k] = -1;
          if (k >= 32 && k < 36 && !has5) continue;
          bool good = true; for (int s = 1; s <= 19; s++) if (W[k][s] < 0) good = false; if (!good) continue;
          for (int s = P0; s <= P1; s++) ped[k] += W[k][s]; ped[k] /= (P1 - P0 + 1); valid[k] = true;
          grp[k] = 2 * asic[CH[k]] + half[CH[k]]; brd[k] = k < 32 ? 2 * (k / 8) + (k / 4) % 2 : (k < 36 ? 8 : 9 + k - 36); if (k < 32 && k / 8 == 3) brd[k] = 20 + 2 * ((k / 4) % 2) + half[CH[k]];   // SFIL4 boards are split over two halves
        }
      // reference from the unused channels (control only)
      for (int g = 0; g < 8; g++) for (int s = 0; s < 20; s++)
        {
          buf.clear();
          for (int u = 0; u < NU; u++)
            {
              const int c = uch[u]; if (c < 0 || c >= 288 || 2 * asic[c] + half[c] != g || uadc[u][s] < 0) continue;
              double p = 0; int np = 0; for (int q = P0; q <= P1; q++) if (uadc[u][q] >= 0) { p += uadc[u][q]; np++; }
              if (np < 3) continue; buf.push_back(uadc[u][s] - p / np);
            }
          if (buf.size() < 5) { cmu[g][s] = 0; continue; }
          std::nth_element(buf.begin(), buf.begin() + buf.size() / 2, buf.end()); cmu[g][s] = buf[buf.size() / 2];
        }
      // reference from the connected channels of the half
      for (int g = 0; g < 8; g++) { for (int s = 0; s < 20; s++) { cm[g][s] = median(g, -1, s); if (g == 0 || g == 1 || g == 4 || g == 5) hCM[g]->Fill(cm[g][s]); } double m = 0, m2 = 0; for (int s = 1; s <= 19; s++) { m += cm[g][s]; m2 += cm[g][s] * cm[g][s]; } m /= 19; hCMrms[g]->Fill(std::sqrt(std::max(0., m2 / 19 - m * m))); }
      padbits = 0; for (int j = 0; j < 4; j++) pdip[j] = -999;
      for (int k = 0; k < NK; k++)
        {
          if (!valid[k]) continue;
          const int g = grp[k]; const double p = ped[k];
          if (k >= 36) { double med[5]; for (int s = 1; s <= 5; s++) med[s - 1] = W[k][s]; std::sort(med, med + 5); int mn = 1023; for (int s = 6; s <= 11; s++) if (W[k][s] < mn) mn = W[k][s]; pdip[k - 36] = med[2] - mn; A[k] = pdip[k - 36]; if (A[k] > 40) padbits |= 1 << (k - 36); continue; }
          double y[20], yx[20], yu[20], r0 = 0, r1 = 0, c0 = 0, c1 = 0;
          for (int s = 0; s < 20; s++) { y[s] = W[k][s] - p - cm[g][s]; yx[s] = W[k][s] - p - median(g, brd[k], s); yu[s] = W[k][s] - p - cmu[g][s]; }
          for (int s = P0; s <= P1; s++) { const double r = W[k][s] - p; r0 += r; r1 += r * r; c0 += y[s]; c1 += y[s] * y[s]; }
          const int np = P1 - P0 + 1; noise[k] = std::sqrt(std::max(0., c1 / np - c0 * c0 / np / np)); const double raw = std::sqrt(std::max(0., r1 / np - r0 * r0 / np / np));
          int at = -1; A[k] = wmax(y, A0, A1, &at); pk[k] = at; Ax[k] = wmax(yx, A0, A1); Au[k] = wmax(yu, A0, A1); Al[k] = wmax(y, L0, L1);
          Q[k] = 0; for (int s = Q0; s <= Q1; s++) Q[k] += y[s];
          hA[k]->Fill(A[k]); hAl[k]->Fill(Al[k]); hNoise[k]->Fill(noise[k]); hRaw[k]->Fill(raw); if (A[k] > 15) hPk[k]->Fill(pk[k]);
          for (int s = 0; s < 20; s++) wf[20 * k + s] = (Short_t) std::lround(10 * std::max(-3000., std::min(3000., y[s])));
        }
      nwf = 0; for (int k = 0; k < 36; k++) if (valid[k] && isRef(k) && A[k] > 12) nwf = 36 * 20;
      if (nwf) for (int k = 0; k < 36; k++) if (!valid[k]) for (int s = 0; s < 20; s++) wf[20 * k + s] = -32000;
      // muon tags
      for (int k = 0; k < 36; k++)
        {
          if (!valid[k]) continue;
          const int L = k < 32 ? k / 8 : 4, e = k < 32 ? (k / 4) % 2 : 1, c = k % 4;
          auto kk = [&](int LL, int ee, int cc) { return LL < 4 ? 8 * LL + 4 * ee + cc : 32 + cc; };
          if (k < 32) { const int o = kk(L, 1 - e, c); if (valid[o] && A[o] > THR_B) hB[k]->Fill(A[k]); if (e == 0 && valid[o] && A[k] > THR_GM && A[o] > THR_GM) hGM[L][c]->Fill(std::sqrt(A[k] * A[o])); }
          bool V = true; if (c < 3) { for (int cc = 0; cc < 3; cc++) if (cc != c) { const int j = kk(L, e, cc); if (!valid[j] || Ax[j] < THR_V) V = false; } } else { for (int cc = 1; cc <= 2; cc++) { const int j = kk(L, e, cc); if (!valid[j] || Ax[j] < THR_V) V = false; } }
          if (V) hV[k]->Fill(Ax[k]);
        }
      trig = i; T->Fill(); hReasons->Fill(2.5);
    }
  printf("triggers %lld-%lld: looked at %.0f, complete %.0f, written %.0f\n", start, stop - 1, hReasons->GetBinContent(1), hReasons->GetBinContent(2), hReasons->GetBinContent(3));
  for (int k = 0; k < 36; k++) if (hA[k]->GetEntries()) printf("  k %2d ch %3d: noise raw %.1f -> %.1f ADC; A > 15: %.2f %% (late window %.2f %%); both-ends tag %.0f, vertical %.0f\n", k, chan[k], hRaw[k]->GetMean(), hNoise[k]->GetMean(), 100. * hA[k]->Integral(hA[k]->FindBin(15), 250) / hA[k]->GetEntries(), 100. * hAl[k]->Integral(hAl[k]->FindBin(15), 250) / hAl[k]->GetEntries(), hB[k]->GetEntries(), hV[k]->GetEntries());
  for (int L = 0; L < 4; L++) printf("  SFIL%d both-ends > %.0f per cell: %.0f %.0f %.0f %.0f\n", L + 1, THR_GM, hGM[L][0]->GetEntries(), hGM[L][1]->GetEntries(), hGM[L][2]->GetEntries(), hGM[L][3]->GetEntries());
  fout->Write(); fout->Close(); printf("wrote %s\n", out);
}
