// gain_stability.C: response versus time in the run from the template-fit trees (kcu0_sfil_fit.C / kcu1_sfil5_fit.C output), one
// file (= one part of run 3064) per call; the histograms of all parts are hadd-ed and plotted by gain_stability_figs.py.
//
//   root -l -b -q 'gain_stability.C+("kcu0", fitfile, datafile, out, t0, tmax)'
//   root -l -b -q 'gain_stability.C+("kcu1", fitfile, "", out, t0, tmax)'
//
// Time axis: the board counter ts0 x 6.25 ns minus t0 (seconds; the counter does not restart at the run start), bins of 600 s up to
// tmax. KCU0: ts0 comes from the data file's header tree (entry = trig of the fit tree); KCU1: ts0 is in the fit tree.
// KCU0 output, per paddle p (0 = any single paddle, 3-6) : hN_p (triggers with a fitted beam-cell pulse), h2_k_p (A_fit of channel k
// vs time, k = 8 layer + 4 end + cell), hS_L_p (South layer sum of the ok cells, beam cell required), hT_p (SFIL1-4 South total),
// hNS_L_p (ln A_N / A_S of the beam cell, both ends ok, neither saturated-recovered above 900 ADC); per channel k: pPed_k (pedestal vs
// time, all fitted waveforms), pSat_k (fraction of pulses whose raw maximum is at the 1023 ceiling, amp + ped >= 1022).
// KCU1 output (no paddle on KCU1): hN, h2_k (A_fit of the 64 channels, unclipped pulses only), pPed_k, pClip_k (fraction with nsat > 0).
#include "TFile.h"
#include "TTree.h"
#include "TH1D.h"
#include "TH2F.h"
#include "TProfile.h"
#include "TString.h"
#include <cmath>
#include <cstdio>

void gain_stability(const char *board, const char *fitfile, const char *datafile, const char *out, double t0 = 0., double tmax = 45000.)
{
  const double TICK = 6.25e-9; const int NT = int(tmax / 600.) + 1;
  TFile *ff = TFile::Open(fitfile); if (!ff || ff->IsZombie()) { printf("cannot open %s\n", fitfile); return; }
  TTree *T = (TTree *) ff->Get("fits"); const Long64_t n = T->GetEntries();
  TFile *fo = new TFile(out, "RECREATE");
  if (TString(board) == "kcu0")
    {
      TFile *fd = TFile::Open(datafile); if (!fd || fd->IsZombie()) { printf("cannot open %s\n", datafile); return; }
      TTree *H = (TTree *) fd->Get("header"); Long64_t ts0; H->SetBranchStatus("*", 0); H->SetBranchStatus("ts0", 1); H->SetBranchAddress("ts0", &ts0);
      Long64_t trig; Int_t pad; Bool_t ok[32]; Float_t amp[32], afit[32], ped[32];
      T->SetBranchStatus("*", 0); for (const char *b : {"trig", "pad", "ok", "amp", "afit", "ped"}) T->SetBranchStatus(b, 1);
      T->SetBranchAddress("trig", &trig); T->SetBranchAddress("pad", &pad); T->SetBranchAddress("ok", ok); T->SetBranchAddress("amp", amp); T->SetBranchAddress("afit", afit); T->SetBranchAddress("ped", ped);
      fo->cd();
      const int P[5] = {0, 3, 4, 5, 6}; TH1D *hN[5]; TH2F *h2[32][5], *hS[4][5], *hT[5], *hNS[4][5]; TProfile *pPed[32], *pSat[32];
      for (int p = 0; p < 5; p++)
        {
          hN[p] = new TH1D(Form("hN_%d", P[p]), Form("paddle %d: triggers with a fitted beam-cell pulse;time in run [s];triggers / 600 s", P[p]), NT, 0, NT * 600.);
          hT[p] = new TH2F(Form("hT_%d", P[p]), Form("paddle %d: SFIL1-4 South total;time in run [s];#Sigma A_{fit} [ADC]", P[p]), NT, 0, NT * 600., 150, 0, 9000);
          for (int k = 0; k < 32; k++) h2[k][p] = new TH2F(Form("h2_%d_%d", k, P[p]), Form("k %d, paddle %d;time in run [s];A_{fit} [ADC]", k, P[p]), NT, 0, NT * 600., 150, 0, 3000);
          for (int L = 0; L < 4; L++)
            {
              hS[L][p] = new TH2F(Form("hS_%d_%d", L + 1, P[p]), Form("SFIL%d South layer sum, paddle %d;time in run [s];A_{fit} sum [ADC]", L + 1, P[p]), NT, 0, NT * 600., 150, 0, 4500);
              hNS[L][p] = new TH2F(Form("hNS_%d_%d", L + 1, P[p]), Form("SFIL%d beam cell ln(A_{N}/A_{S}), paddle %d;time in run [s];ln(A_{N}/A_{S})", L + 1, P[p]), NT, 0, NT * 600., 100, -2, 2);
            }
        }
      for (int k = 0; k < 32; k++)
        {
          pPed[k] = new TProfile(Form("pPed_%d", k), Form("k %d pedestal;time in run [s];pedestal [ADC]", k), NT, 0, NT * 600.);
          pSat[k] = new TProfile(Form("pSat_%d", k), Form("k %d clipped fraction;time in run [s];fraction at the ceiling", k), NT, 0, NT * 600.);
        }
      Long64_t used = 0;
      for (Long64_t i = 0; i < n; i++)
        {
          T->GetEntry(i); H->GetEntry(trig); const double t = ts0 * TICK - t0; if (t < 0 || t >= NT * 600.) continue;
          const int ip = (pad >= 3 && pad <= 6) ? pad - 2 : -1; used++;
          for (int k = 0; k < 32; k++) if (ok[k]) { pPed[k]->Fill(t, ped[k]); pSat[k]->Fill(t, amp[k] + ped[k] >= 1022 ? 1. : 0.); }
          double tot = 0; bool anyBeam = false;
          for (int L = 0; L < 4; L++)
            {
              const int kb = 8 * L + 2, kn = 8 * L + 4 + 2;
              if (!ok[kb]) continue; anyBeam = true;
              double s = 0; for (int c = 0; c < 4; c++) if (ok[8 * L + c]) s += afit[8 * L + c]; tot += s;
              for (int p = 0; p < 5; p++) { if (p && p != ip) continue; hS[L][p]->Fill(t, s); if (L && ok[kn] && afit[kb] > 100 && afit[kn] > 20 && amp[kb] < 900 && amp[kn] < 900) hNS[L][p]->Fill(t, std::log(afit[kn] / afit[kb])); }
            }
          if (!anyBeam) continue;
          for (int p = 0; p < 5; p++) { if (p && p != ip) continue; hN[p]->Fill(t); hT[p]->Fill(t, tot); for (int k = 0; k < 32; k++) if (ok[k]) h2[k][p]->Fill(t, afit[k]); }
        }
      printf("%s: %lld fit entries, %lld in the time window\n", fitfile, n, used);
    }
  else
    {
      Long64_t ts0; Bool_t ok[64]; Float_t afit[64], ped[64]; Int_t nsat[64];
      T->SetBranchStatus("*", 0); for (const char *b : {"ts0", "ok", "afit", "ped", "nsat"}) T->SetBranchStatus(b, 1);
      T->SetBranchAddress("ts0", &ts0); T->SetBranchAddress("ok", ok); T->SetBranchAddress("afit", afit); T->SetBranchAddress("ped", ped); T->SetBranchAddress("nsat", nsat);
      fo->cd();
      TH1D *hN = new TH1D("hN", "triggers with any fitted pulse;time in run [s];triggers / 600 s", NT, 0, NT * 600.); TH2F *h2[64]; TProfile *pPed[64], *pClip[64];
      for (int k = 0; k < 64; k++)
        {
          h2[k] = new TH2F(Form("h2_%d", k), Form("k %d unclipped;time in run [s];A_{fit} [ADC]", k), NT, 0, NT * 600., 150, 0, 1500);
          pPed[k] = new TProfile(Form("pPed_%d", k), Form("k %d pedestal;time in run [s];pedestal [ADC]", k), NT, 0, NT * 600.);
          pClip[k] = new TProfile(Form("pClip_%d", k), Form("k %d clipped fraction;time in run [s];fraction with nsat > 0", k), NT, 0, NT * 600.);
        }
      Long64_t used = 0;
      for (Long64_t i = 0; i < n; i++)
        {
          T->GetEntry(i); const double t = ts0 * TICK - t0; if (t < 0 || t >= NT * 600.) continue; used++;
          bool any = false;
          for (int k = 0; k < 64; k++) if (ok[k]) { any = true; pPed[k]->Fill(t, ped[k]); pClip[k]->Fill(t, nsat[k] > 0 ? 1. : 0.); if (!nsat[k]) h2[k]->Fill(t, afit[k]); }
          if (any) hN->Fill(t);
        }
      printf("%s: %lld fit entries, %lld in the time window\n", fitfile, n, used);
    }
  fo->Write(); fo->Close(); printf("wrote %s\n", out);
}
