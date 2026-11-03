// bic_streams.C: ROOT macro version of bic_streams.py. Reads the run-3062 stream trees (KCU0, KCU1, CAEN,
// AstroPix) together, matches triggers between streams and offers waveform tools.
//
//   root -l
//   .L bic_streams.C+                      // compile with ACLiC
//   BicStreams S("../data"); S.Index();    // per-trigger clocks of the whole run (cached in cache/index.root)
//   auto j1 = S.MatchKCU1();               // KCU0 entry -> KCU1 entry (-1 if none), exact tick offset
//   auto jc = S.MatchCAEN();               // KCU0 entry -> CAEN entry, linear clock model + running median
//   Short_t Sw[4][20], Nw[4][20]; S.SfilWaveforms(123456, 3, Sw, Nw);
//   Pulse p = PulseFeatures(Nw[2], 20);    // pedestal, amplitude, crossing time, integral, saturation
//   TF1 *f = FitPulse(Nw[2], 20);          // rise^p x exp decay model, saturated samples masked
//   bic_streams()                          // the demo: match, print, draw demo_event_root.png
//
// Clock constants (docs/method.md, docs/root_files.md): KCU tick 6.25 ns; ts(KCU1) = ts(KCU0) + KCU1_OFFSET
// (+0 or +4); KCU0 ticks = CAEN_RATIO x gtt + b(t); AstroPix on-chip ts = ((K+4)/16 - ASTRO_PHASE) mod 256.
#include "TFile.h"
#include "TTree.h"
#include "TH1D.h"
#include "TF1.h"
#include "TGraph.h"
#include "TCanvas.h"
#include "TLegend.h"
#include "TLatex.h"
#include "TStyle.h"
#include "TSystem.h"
#include "TMath.h"
#include <vector>
#include <algorithm>
#include <cstdio>
#include <cmath>

static const double TICK_NS = 6.25;
static const Long64_t KCU1_OFFSET = -1143665946;
static const double CAEN_RATIO = 1.365340679;
static const int ASTRO_PHASE = 189;
static const double TOA_WRAP_NS = 6.37;   // the LED macro's "magic number" (its configuration); only used in macro mode
// nominal TOA code -> ns (fine x 25 ps + mid x 0.2 ns + coarse x 6.25 ns), and the time (ns from the start of sample 0)
// of a TOA word found in `sample`: default 25 s + TOA_ns (what run 3062 supports, method.md 4.9);
// macro = true: the LED macro's rule 25 (s + (TOA_ns < wrap ? 0 : 1)) - TOA_ns
double ToaNs(int code) { return (code & 7) * 0.025 + ((code >> 3) & 31) * 0.2 + ((code >> 8) & 3) * 6.25; }
double ToaTime(int sample, int code, bool macro = false, double wrap = TOA_WRAP_NS) { const double t = ToaNs(code); return macro ? (t < wrap ? sample : sample + 1) * 25.0 - t : 25.0 * sample + t; }

// ------------------------------------------------------------------ waveform tools
struct Pulse { double ped, amp, tcross, integral; bool saturated; };

// pedestal from the first nped samples; amplitude = max above pedestal; crossing time (samples, linear
// interpolation) at frac x amplitude; integral above pedestal; negative = true for paddles and the CAEN
Pulse PulseFeatures(const Short_t *wf, int n, int nped = 4, double frac = 0.2, bool negative = false, int sat = 1023)
{
  Pulse p{0, 0, -1, 0, false}; const double sgn = negative ? -1 : 1;
  for (int i = 0; i < nped; i++) p.ped += sgn * wf[i]; p.ped /= nped;
  int imax = 0;
  for (int i = 0; i < n; i++) { const double v = sgn * wf[i] - p.ped; p.integral += v; if (v > p.amp) { p.amp = v; imax = i; } if (wf[i] == (negative ? 0 : sat)) p.saturated = true; }
  const double thr = frac * p.amp;
  for (int k = 1; k <= imax && p.amp > 0; k++)
    {
      const double y0 = sgn * wf[k - 1] - p.ped, y1 = sgn * wf[k] - p.ped;
      if (y1 >= thr) { p.tcross = (y1 != y0) ? k - 1 + (thr - y0) / (y1 - y0) : k; break; }
    }
  p.ped *= sgn;   // back to ADC units
  return p;
}

// A x ((t - t0)/tau_r)^p x exp(-(t - t0 - p tau_r)/tau_f) for t > t0;  p = {A, t0, tau_r, tau_f, p}
double PulseModel(double *x, double *par)
{
  const double t = x[0]; if (t <= par[1]) return 0;
  const double u = (t - par[1]) / par[2];
  return par[0] * std::pow(u, par[4]) * std::exp(-(t - par[1] - par[4] * par[2]) / par[3]);
}

// fit one waveform (dt = sample spacing in ns), ignoring saturated samples; returns the TF1 (caller owns) or nullptr
TF1 *FitPulse(const Short_t *wf, int n, double dt = 25., int nped = 4, int sat = 1023, bool negative = false)
{
  const double sgn = negative ? -1 : 1; double ped = 0; for (int i = 0; i < nped; i++) ped += sgn * wf[i]; ped /= nped;
  TGraph g; int imax = 0; double vmax = -1e9;
  for (int i = 0; i < n; i++) { if (wf[i] == (negative ? 0 : sat)) continue; const double v = sgn * wf[i] - ped; g.SetPoint(g.GetN(), i * dt, v); if (v > vmax) { vmax = v; imax = i; } }
  if (g.GetN() < 6 || vmax < 20) return nullptr;
  TF1 *f = new TF1("pulse", PulseModel, -dt, n * dt, 5);
  f->SetParameters(vmax, std::max(imax * dt - 2 * dt, 0.), dt, 4 * dt, 2.);
  f->SetParLimits(0, 0, 1e5); f->SetParLimits(1, -dt, n * dt); f->SetParLimits(2, 1, 10 * dt); f->SetParLimits(3, 1, 50 * dt); f->SetParLimits(4, 0.5, 8);
  if (g.Fit(f, "QN") != 0) { delete f; return nullptr; }
  return f;
}

// ------------------------------------------------------------------ the streams
class BicStreams
{
public:
  BicStreams(const char *dir = "../data", const char *cache = "cache")
  {
    fDir = dir; fCache = cache; gSystem->mkdir(cache, true);
    fK0 = TFile::Open(Form("%s/kcu0_3062.root", dir)); fK1 = TFile::Open(Form("%s/kcu1_3062.root", dir)); fC = TFile::Open(Form("%s/caen_3062.root", dir)); fA = TFile::Open(Form("%s/astropix_3062.root", dir));
    for (int L = 0; L < 4; L++) { tSfil[L] = fK0 ? (TTree *) fK0->Get(Form("sfil%d", L + 1)) : nullptr; if (tSfil[L]) for (int c = 0; c < 4; c++) { tSfil[L]->SetBranchAddress(Form("Sc%d_adc", c + 1), bS[L][c]); tSfil[L]->SetBranchAddress(Form("Nc%d_adc", c + 1), bN[L][c]); } }
    tPad = fK0 ? (TTree *) fK0->Get("paddles") : nullptr; const char *pn[4] = {"P3", "P4", "P5", "P6"};
    if (tPad) for (int k = 0; k < 4; k++) tPad->SetBranchAddress(Form("%s_adc", pn[k]), bP[k]);
    tC = fC ? (TTree *) fC->Get("bbcal") : nullptr; tCp = fC ? (TTree *) fC->Get("paddles") : nullptr;
    if (tC) tC->SetBranchAddress("adc", bBB); if (tCp) for (int k = 0; k < 4; k++) tCp->SetBranchAddress(Form("%s_adc", pn[k]), bCP[k]);
    tS5 = fK1 ? (TTree *) fK1->Get("sfil5") : nullptr; const char *an[4] = {"A", "B", "C", "D"};
    if (tS5) for (int a = 0; a < 4; a++) tS5->SetBranchAddress(Form("%s_adc", an[a]), bS5[a]);
  }

  // ---- per-trigger clocks of the whole run, cached
  void Index(bool force = false)
  {
    const TString cf = Form("%s/index.root", fCache.Data());
    if (!force && !gSystem->AccessPathName(cf)) { if (ReadCache(cf)) return; }
    Long64_t ts; Bool_t comp; Long64_t rt;
    if (fK0) { TTree *h = (TTree *) fK0->Get("header"); h->SetBranchStatus("*", 0); h->SetBranchStatus("ts0", 1); h->SetBranchStatus("complete", 1); h->SetBranchStatus("rt_first", 1); h->SetBranchAddress("ts0", &ts); h->SetBranchAddress("complete", &comp); h->SetBranchAddress("rt_first", &rt);
      const Long64_t n = h->GetEntries(); k0_ts.resize(n); k0_complete.resize(n); k0_rt.resize(n); for (Long64_t i = 0; i < n; i++) { h->GetEntry(i); k0_ts[i] = ts; k0_complete[i] = comp; k0_rt[i] = rt; } h->SetBranchStatus("*", 1); }
    if (fK1) { TTree *h = (TTree *) fK1->Get("header"); h->SetBranchStatus("*", 0); h->SetBranchStatus("ts0", 1); h->SetBranchStatus("complete", 1); h->SetBranchAddress("ts0", &ts); h->SetBranchAddress("complete", &comp);
      const Long64_t n = h->GetEntries(); k1_ts.resize(n); k1_complete.resize(n); for (Long64_t i = 0; i < n; i++) { h->GetEntry(i); k1_ts[i] = ts; k1_complete[i] = comp; } h->SetBranchStatus("*", 1); }
    if (fC) { TTree *h = (TTree *) fC->Get("header"); Long64_t gtt, crt; h->SetBranchStatus("*", 0); h->SetBranchStatus("gtt", 1); h->SetBranchStatus("rt", 1); h->SetBranchAddress("gtt", &gtt); h->SetBranchAddress("rt", &crt);
      const Long64_t n = h->GetEntries(); c_gtt.resize(n); c_rt.resize(n); for (Long64_t i = 0; i < n; i++) { h->GetEntry(i); c_gtt[i] = gtt; c_rt[i] = crt; } h->SetBranchStatus("*", 1); }
    if (fA) { TTree *h = (TTree *) fA->Get("hits"); UInt_t fp; UChar_t ats, id; Bool_t val; h->SetBranchStatus("*", 0); for (const char *b : {"fpga_ts", "ts", "id", "valid"}) h->SetBranchStatus(b, 1);
      h->SetBranchAddress("fpga_ts", &fp); h->SetBranchAddress("ts", &ats); h->SetBranchAddress("id", &id); h->SetBranchAddress("valid", &val);
      const Long64_t n = h->GetEntries(); a_fpga.resize(n); a_ts.resize(n); a_chip.resize(n); a_valid.resize(n); Long64_t wraps = 0, prev = 0;
      for (Long64_t i = 0; i < n; i++) { h->GetEntry(i); if (i && (Long64_t) fp - prev < -(1LL << 31)) wraps++; prev = fp; a_fpga[i] = (Long64_t) fp + wraps * (1LL << 32); a_ts[i] = ats; a_chip[i] = id; a_valid[i] = val; } h->SetBranchStatus("*", 1); }
    WriteCache(cf);
    printf("index: KCU0 %zu, KCU1 %zu, CAEN %zu triggers, AstroPix %zu hits\n", k0_ts.size(), k1_ts.size(), c_gtt.size(), a_fpga.size());
  }

  // ---- KCU0 <-> KCU1 by the exact tick offset
  std::vector<Long64_t> MatchKCU1(int tol = 8, Long64_t offset = KCU1_OFFSET)
  {
    std::vector<Long64_t> order(k1_ts.size()); for (size_t i = 0; i < order.size(); i++) order[i] = i;
    std::sort(order.begin(), order.end(), [&](Long64_t a, Long64_t b) { return k1_ts[a] < k1_ts[b]; });
    std::vector<Long64_t> ks(order.size()); for (size_t i = 0; i < order.size(); i++) ks[i] = k1_ts[order[i]];
    std::vector<Long64_t> out(k0_ts.size(), -1); Long64_t nok = 0;
    for (size_t i = 0; i < k0_ts.size(); i++)
      {
        const Long64_t pred = k0_ts[i] + offset; auto it = std::lower_bound(ks.begin(), ks.end(), pred); Long64_t j = it - ks.begin(); if (j >= (Long64_t) ks.size()) j = ks.size() - 1; if (j > 0 && std::llabs(ks[j - 1] - pred) < std::llabs(ks[j] - pred)) j--;
        if (std::llabs(ks[j] - pred) <= tol) { out[i] = order[j]; nok++; }
      }
    printf("KCU1 partners: %lld of %zu KCU0 triggers (%.1f%%)\n", nok, k0_ts.size(), 100. * nok / k0_ts.size()); return out;
  }

  // offset b - a (in KCU ticks) maximising coincidences: coarse then fine histogram of pairwise differences
  double FindOffset(const std::vector<double> &a, const std::vector<double> &b, double guess, double window, double fineBin, double coarseBin, size_t n, double *sig = nullptr)
  {
    TH1D hc("hc", "", int(2 * window / coarseBin), guess - window, guess + window);
    const size_t na = std::min(n, a.size());
    for (size_t i = 0; i < na; i++) { auto lo = std::lower_bound(b.begin(), b.end(), a[i] + guess - window), hi = std::upper_bound(b.begin(), b.end(), a[i] + guess + window); for (auto it = lo; it != hi; ++it) hc.Fill(*it - a[i]); }
    const double c0 = hc.GetBinCenter(hc.GetMaximumBin());
    TH1D hf("hf", "", int(4 * coarseBin / fineBin), c0 - 2 * coarseBin, c0 + 2 * coarseBin);
    for (size_t i = 0; i < na; i++) { auto lo = std::lower_bound(b.begin(), b.end(), a[i] + c0 - 2 * coarseBin), hi = std::upper_bound(b.begin(), b.end(), a[i] + c0 + 2 * coarseBin); for (auto it = lo; it != hi; ++it) hf.Fill(*it - a[i]); }
    if (sig) { std::vector<double> v; for (int k = 1; k <= hf.GetNbinsX(); k++) v.push_back(hf.GetBinContent(k)); std::nth_element(v.begin(), v.begin() + v.size() / 2, v.end()); *sig = hf.GetMaximum() / std::max(v[v.size() / 2], 1.); }
    return hf.GetBinCenter(hf.GetMaximumBin());
  }

  // ---- KCU0 <-> CAEN: linear clock model (ratio + offset + residual slope from scans at both ends of the run),
  // then the few-us wander followed by a running median in two passes (20 us, 2 us); final tolerance tol_ns
  std::vector<Long64_t> MatchCAEN(double tol_ns = 100., double ratio = CAEN_RATIO, int block = 201)
  {
    std::vector<double> k0(k0_ts.begin(), k0_ts.end()), g(c_gtt.size()); for (size_t i = 0; i < g.size(); i++) g[i] = c_gtt[i] * ratio;
    double b0, slope; const TString cf = Form("%s/caen_offset.txt", fCache.Data()); FILE *fc = fopen(cf, "r");
    if (fc && fscanf(fc, "%lf %lf", &b0, &slope) == 2) fclose(fc);
    else
      {
        if (fc) fclose(fc);
        const size_t m = std::min<size_t>(20000, std::min(k0.size(), g.size())); std::vector<double> d1(m), d2(m);
        for (size_t i = 0; i < m; i++) { d1[i] = k0_rt[i] - k0[i] * TICK_NS; d2[i] = c_rt[i] - g[i] * TICK_NS; }
        std::nth_element(d1.begin(), d1.begin() + m / 2, d1.end()); std::nth_element(d2.begin(), d2.begin() + m / 100, d2.end());
        const double rtGuess = d1[m / 2] - d2[m / 100];   // CAEN lags up to 2 s: take its lower envelope
        double sig0 = 0, sig1 = 0;
        b0 = FindOffset(g, k0, -rtGuess / TICK_NS, 3e9 / TICK_NS, 4., 160000., 5000, &sig0);
        std::vector<double> gtail(g.end() - 6000, g.end());
        const double b1 = FindOffset(gtail, k0, b0, 2e6 / TICK_NS, 4., 800., 6000, &sig1);
        slope = (b1 - b0) / (g[g.size() - 3000] - g[3000]);
        printf("CAEN offset: %.6f s at the start (peak/background %.0f), %+.1f us later at the end (peak %.0f) -> residual slope %+.3f ppm\n", b0 * TICK_NS * 1e-9, sig0, (b1 - b0) * TICK_NS * 1e-3, sig1, slope * 1e6);
        fc = fopen(cf, "w"); if (fc) { fprintf(fc, "%.6f %.12e\n", b0, slope); fclose(fc); }
      }
    std::vector<double> pred(g.size()); for (size_t i = 0; i < g.size(); i++) pred[i] = g[i] + slope * (g[i] - g[3000]) + b0;
    auto nearest = [&](const std::vector<double> &p, std::vector<Long64_t> &near, std::vector<double> &res)
      { near.resize(p.size()); res.resize(p.size()); for (size_t i = 0; i < p.size(); i++) { auto it = std::lower_bound(k0.begin(), k0.end(), p[i]); Long64_t j = it - k0.begin(); if (j >= (Long64_t) k0.size()) j = k0.size() - 1; if (j > 0 && std::fabs(k0[j - 1] - p[i]) < std::fabs(k0[j] - p[i])) j--; near[i] = j; res[i] = (k0[j] - p[i]) * TICK_NS; } };
    std::vector<Long64_t> near; std::vector<double> res;
    for (double tol : {20000., 2000.})
      {
        nearest(pred, near, res); std::vector<size_t> idx; for (size_t i = 0; i < res.size(); i++) if (std::fabs(res[i]) < tol) idx.push_back(i);
        if ((int) idx.size() < block) break;
        std::vector<double> rm(idx.size()), win;
        for (size_t k = 0; k < idx.size(); k++) { const size_t lo = k >= (size_t) block / 2 ? k - block / 2 : 0, hi = std::min(idx.size(), k + block / 2 + 1); win.assign(hi - lo, 0); for (size_t q = lo; q < hi; q++) win[q - lo] = res[idx[q]]; std::nth_element(win.begin(), win.begin() + win.size() / 2, win.end()); rm[k] = win[win.size() / 2]; }
        size_t k = 0;
        for (size_t i = 0; i < pred.size(); i++) { while (k + 1 < idx.size() && idx[k + 1] <= i) k++; double corr = rm[k]; if (k + 1 < idx.size() && i > idx[k]) corr = rm[k] + (rm[k + 1] - rm[k]) * double(i - idx[k]) / double(idx[k + 1] - idx[k]); pred[i] += corr / TICK_NS; }
      }
    nearest(pred, near, res);
    std::vector<Long64_t> out(k0_ts.size(), -1); Long64_t nok = 0; std::vector<double> good;
    for (size_t i = 0; i < res.size(); i++) if (std::fabs(res[i]) < tol_ns) { out[near[i]] = i; nok++; good.push_back(res[i]); }
    double mad = 0; if (!good.empty()) { std::vector<double> t(good); std::nth_element(t.begin(), t.begin() + t.size() / 2, t.end()); const double med = t[t.size() / 2]; for (auto &v : t) v = std::fabs(v - med); std::nth_element(t.begin(), t.begin() + t.size() / 2, t.end()); mad = t[t.size() / 2]; }
    printf("CAEN partners: %lld of %zu CAEN triggers within %.0f ns (%.1f%%), residual MAD %.1f ns\n", nok, g.size(), tol_ns, 100. * nok / g.size(), mad);
    return out;
  }

  // ---- AstroPix on-chip phase test for a KCU0 trigger at tick K: 0 for a true partner
  static int AstroPhase(Long64_t K, int ts) { return int((((K + 4) / 16 - ts - ASTRO_PHASE + 128) % 256 + 256) % 256) - 128; }

  // ---- waveforms of one trigger (-1 = line missing)
  bool SfilWaveforms(Long64_t i, int layer, Short_t S[4][20], Short_t N[4][20]) { TTree *t = tSfil[layer - 1]; if (!t || i < 0 || i >= t->GetEntries()) return false; t->GetEntry(i); for (int c = 0; c < 4; c++) for (int s = 0; s < 20; s++) { S[c][s] = bS[layer - 1][c][s]; N[c][s] = bN[layer - 1][c][s]; } return true; }
  bool PaddleWaveforms(Long64_t i, Short_t P[4][20]) { if (!tPad || i < 0 || i >= tPad->GetEntries()) return false; tPad->GetEntry(i); for (int k = 0; k < 4; k++) for (int s = 0; s < 20; s++) P[k][s] = bP[k][s]; return true; }
  bool CaenWaveforms(Long64_t j, Short_t bb[12][256], Short_t pad[4][256]) { if (!tC || j < 0 || j >= tC->GetEntries()) return false; tC->GetEntry(j); tCp->GetEntry(j); for (int c = 0; c < 12; c++) for (int s = 0; s < 256; s++) bb[c][s] = bBB[c][s]; for (int k = 0; k < 4; k++) for (int s = 0; s < 256; s++) pad[k][s] = bCP[k][s]; return true; }
  bool Sfil5Waveforms(Long64_t j, Short_t arr[4][16][20]) { if (!tS5 || j < 0 || j >= tS5->GetEntries()) return false; tS5->GetEntry(j); for (int a = 0; a < 4; a++) for (int i = 0; i < 16; i++) for (int s = 0; s < 20; s++) arr[a][i][s] = bS5[a][i][s]; return true; }

  std::vector<Long64_t> k0_ts, k1_ts, c_gtt, k0_rt, c_rt, a_fpga; std::vector<Bool_t> k0_complete, k1_complete, a_valid; std::vector<UChar_t> a_ts, a_chip;

private:
  void WriteCache(const TString &cf)
  {
    TFile f(cf, "RECREATE"); TTree t("index", "per-trigger clocks of run 3062");
    auto *p0 = &k0_ts; auto *p1 = &k0_complete; auto *p2 = &k0_rt; auto *p3 = &k1_ts; auto *p4 = &k1_complete; auto *p5 = &c_gtt; auto *p6 = &c_rt; auto *p7 = &a_fpga; auto *p8 = &a_ts; auto *p9 = &a_chip; auto *p10 = &a_valid;
    t.Branch("k0_ts", &p0); t.Branch("k0_complete", &p1); t.Branch("k0_rt", &p2); t.Branch("k1_ts", &p3); t.Branch("k1_complete", &p4); t.Branch("c_gtt", &p5); t.Branch("c_rt", &p6); t.Branch("a_fpga", &p7); t.Branch("a_ts", &p8); t.Branch("a_chip", &p9); t.Branch("a_valid", &p10);
    t.Fill(); t.Write(); f.Close();
  }
  bool ReadCache(const TString &cf)
  {
    TFile f(cf); TTree *t = (TTree *) f.Get("index"); if (!t) return false;
    auto *p0 = &k0_ts; auto *p1 = &k0_complete; auto *p2 = &k0_rt; auto *p3 = &k1_ts; auto *p4 = &k1_complete; auto *p5 = &c_gtt; auto *p6 = &c_rt; auto *p7 = &a_fpga; auto *p8 = &a_ts; auto *p9 = &a_chip; auto *p10 = &a_valid;
    t->SetBranchAddress("k0_ts", &p0); t->SetBranchAddress("k0_complete", &p1); t->SetBranchAddress("k0_rt", &p2); t->SetBranchAddress("k1_ts", &p3); t->SetBranchAddress("k1_complete", &p4); t->SetBranchAddress("c_gtt", &p5); t->SetBranchAddress("c_rt", &p6); t->SetBranchAddress("a_fpga", &p7); t->SetBranchAddress("a_ts", &p8); t->SetBranchAddress("a_chip", &p9); t->SetBranchAddress("a_valid", &p10);
    t->GetEntry(0); printf("index from cache: KCU0 %zu, KCU1 %zu, CAEN %zu triggers, AstroPix %zu hits\n", k0_ts.size(), k1_ts.size(), c_gtt.size(), a_fpga.size()); return !k0_ts.empty();
  }
  TString fDir, fCache; TFile *fK0 = nullptr, *fK1 = nullptr, *fC = nullptr, *fA = nullptr;
  TTree *tSfil[4] = {nullptr, nullptr, nullptr, nullptr}, *tPad = nullptr, *tC = nullptr, *tCp = nullptr, *tS5 = nullptr;
  Short_t bS[4][4][20], bN[4][4][20], bP[4][20], bBB[12][256], bCP[4][256], bS5[4][16][20];
};

// ------------------------------------------------------------------ demo
void bic_streams(const char *dir = "../data")
{
  BicStreams S(dir); S.Index();
  auto j1 = S.MatchKCU1(); auto jc = S.MatchCAEN();
  std::vector<Long64_t> both; for (size_t i = 0; i < j1.size(); i++) if (j1[i] >= 0 && jc[i] >= 0 && S.k0_complete[i]) both.push_back(i);
  printf("KCU0 triggers with KCU1 and CAEN partners and a complete waveform: %zu\n", both.size());
  if (both.empty()) return;
  const Long64_t i = both[both.size() / 2];
  Short_t Sw[4][20], Nw[4][20]; S.SfilWaveforms(i, 3, Sw, Nw);
  printf("trigger %lld: SFIL3 South amp", i); for (int c = 0; c < 4; c++) printf(" %.0f", PulseFeatures(Sw[c], 20).amp);
  printf("  North amp"); for (int c = 0; c < 4; c++) printf(" %.0f", PulseFeatures(Nw[c], 20).amp);
  printf("  crossing (samples) S"); for (int c = 0; c < 4; c++) printf(" %.2f", PulseFeatures(Sw[c], 20).tcross); printf(" N"); for (int c = 0; c < 4; c++) printf(" %.2f", PulseFeatures(Nw[c], 20).tcross); printf("\n");
  TF1 *fit = FitPulse(Nw[2], 20); if (fit) printf("North cell 3 fit: A %.1f t0 %.1f tau_r %.1f tau_f %.1f p %.2f\n", fit->GetParameter(0), fit->GetParameter(1), fit->GetParameter(2), fit->GetParameter(3), fit->GetParameter(4));
  Short_t bb[12][256], cp[4][256]; const bool hasC = S.CaenWaveforms(jc[i], bb, cp); if (hasC) printf("CAEN partner: paddle P4 min %d\n", *std::min_element(cp[1], cp[1] + 256));
  Long64_t nph = 0, nsel = 0; const Long64_t K = S.k0_ts[i]; for (size_t h = 0; h < S.a_ts.size(); h++) if (S.a_valid[h] && S.a_fpga[h] * 4 > K - 1) { nsel++; if (BicStreams::AstroPhase(K, S.a_ts[h]) == 0) nph++; }
  printf("AstroPix hits whose on-chip phase matches this trigger (whole run, no latency cut): %lld of %lld\n", nph, nsel);
  // picture
  gStyle->SetOptStat(0); gStyle->SetOptTitle(0); gStyle->SetPadLeftMargin(0.14); gStyle->SetPadBottomMargin(0.13); gStyle->SetTitleSize(0.05, "xy"); gStyle->SetLabelSize(0.045, "xy");
  TCanvas *c = new TCanvas("c", "", 1600, 500); c->Divide(3, 1);
  auto ped = [](const Short_t *w, int n, int nped, double *out) { double p = 0; for (int k = 0; k < nped; k++) p += w[k]; p /= nped; for (int k = 0; k < n; k++) out[k] = w[k] - p; };
  double x20[20], y[20], x256[256], y256[256]; for (int k = 0; k < 20; k++) x20[k] = k; for (int k = 0; k < 256; k++) x256[k] = k;
  c->cd(1); TH1 *fr = gPad->DrawFrame(-0.5, -50, 19.5, 1100, ";sample (25 ns);ADC - pedestal"); TLegend *l1 = new TLegend(0.2, 0.68, 0.55, 0.88); l1->SetBorderSize(0); l1->SetFillStyle(0);
  ped(Sw[2], 20, 4, y); TGraph *gs = new TGraph(20, x20, y); gs->SetLineColor(2); gs->SetMarkerColor(2); gs->SetMarkerStyle(20); gs->Draw("LP"); l1->AddEntry(gs, "SFIL3 South cell 3", "lp");
  ped(Nw[2], 20, 4, y); TGraph *gn = new TGraph(20, x20, y); gn->SetLineColor(4); gn->SetMarkerColor(4); gn->SetMarkerStyle(21); gn->Draw("LP"); l1->AddEntry(gn, "SFIL3 North cell 3", "lp");
  if (fit) { TF1 *fd = (TF1 *) fit->Clone("fd"); fd->SetRange(0, 475); fd->SetLineColor(4); fd->SetLineStyle(2); TGraph *gf = new TGraph(); for (int k = 0; k < 96; k++) gf->SetPoint(k, k * 5. / 25., fd->Eval(k * 5.)); gf->SetLineColor(4); gf->SetLineStyle(2); gf->Draw("L"); l1->AddEntry(gf, "pulse fit (North)", "l"); }
  l1->Draw(); (void) fr;
  c->cd(2); gPad->DrawFrame(-0.5, -50, 19.5, 1100, ";sample (25 ns);ADC - pedestal"); TLegend *l2 = new TLegend(0.2, 0.6, 0.6, 0.88); l2->SetBorderSize(0); l2->SetFillStyle(0);
  const int col[4] = {1, 2, 4, 6}, mk[4] = {20, 21, 22, 23};
  for (int L = 1; L <= 4; L++) { Short_t s2[4][20], n2[4][20]; S.SfilWaveforms(i, L, s2, n2); ped(n2[2], 20, 4, y); TGraph *g = new TGraph(20, x20, y); g->SetLineColor(col[L - 1]); g->SetMarkerColor(col[L - 1]); g->SetMarkerStyle(mk[L - 1]); g->Draw("LP"); l2->AddEntry(g, Form("SFIL%d North cell 3", L), "lp"); }
  l2->Draw();
  c->cd(3); gPad->DrawFrame(0, -3000, 256, 300, ";sample (1 ns);CAEN ADC - pedestal"); TLegend *l3 = new TLegend(0.2, 0.2, 0.5, 0.45); l3->SetBorderSize(0); l3->SetFillStyle(0);
  if (hasC) { const char *pn[4] = {"P3", "P4", "P5", "P6"}; for (int k = 0; k < 4; k++) { ped(cp[k], 256, 20, y256); TGraph *g = new TGraph(256, x256, y256); g->SetLineColor(col[k]); g->Draw("L"); l3->AddEntry(g, Form("CAEN %s", pn[k]), "l"); } l3->Draw(); }
  c->SaveAs("demo_event_root.png"); printf("wrote demo_event_root.png\n");
}
