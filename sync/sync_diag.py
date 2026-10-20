"""Event-synchronization diagnostics vs run time, run 3062 -> out/sync_diag.png

 1 recorded trigger rates per stream and the union (10 s bins)
 2 pairwise overlap: share of each stream's triggers also seen by another stream
 3 KCU0 - KCU1 tick residual: share at 0, -4 (KCU1 one bunch late) and anything else
 4 CAEN - KCU0 residual [ns] vs time (matched pairs)
 5 residual distributions: CAEN - KCU0, CAEN - KCU1 on triggers without KCU0 (independent of KCU0)
 6 paddle pattern agreement KCU0 vs CAEN on shared triggers, real and neighbour-trigger decoy
 7 AstroPix 10 MHz phase lock vs time: signal-band minus control-band candidates per cluster
 8 AstroPix readout latency of log10LR > 2.5 matches vs time
 9 AstroPix clusters and log10LR > 2.5 matches per 10 s"""
import numpy as np, uproot, glob, ROOT
from plotstyle import setup, finish, legend, headroom, colz
setup()
TT = np.load('out/triggers_3062.npz'); T = TT['T']; t0 = T[0]; sec = lambda x: (x - t0) * 6.25e-9
i0, i1, ic = TT['i_kcu0'], TT['i_kcu1'], TT['i_caen']; h0, h1, hc = i0 >= 0, i1 >= 0, ic >= 0
E = np.arange(0, 1500, 10.0); tc = 0.5 * (E[1:] + E[:-1]); ts = sec(T)
H = lambda m: np.histogram(ts[m], E)[0].astype(float)
keep = []; cols = (ROOT.kBlack, ROOT.kRed + 1, ROOT.kBlue + 1, ROOT.kGreen + 2, ROOT.kMagenta + 1, ROOT.kOrange + 7)
cv = ROOT.TCanvas("cv", "", 2100, 1500); cv.Divide(3, 3, 0.004, 0.006)

def graphs(pad, title, series, ymin, ymax, leg=(0.6, 0.72, 0.9, 0.9), logy=False):
    cv.cd(pad); ROOT.gPad.SetGrid(); ROOT.gPad.SetLogy(logy)
    fr = ROOT.TH1D(f"fr{pad}", title, 1, 0, 1500); fr.SetMinimum(ymin); fr.SetMaximum(ymax); fr.Draw(); keep.append(fr)
    lg = legend(*leg); lg.SetFillColor(0); lg.SetNColumns(2 if len(series) > 3 else 1); keep.append(lg)
    for k, (name, x, y) in enumerate(series):
        ok = np.isfinite(y); g = ROOT.TGraph(int(ok.sum()), x[ok].astype(float), y[ok].astype(float))
        g.SetLineColor(cols[k]); g.SetMarkerColor(cols[k]); g.SetLineWidth(2); g.Draw("L same"); keep.append(g); lg.AddEntry(g, name, "l")
    lg.Draw()

# 1 rates
graphs(1, "Recorded triggers per stream;run time [s];rate [Hz]",
       [("union", tc, H(np.ones(len(T), bool)) / 10), ("KCU0 complete", tc, H(h0) / 10), ("CAEN", tc, H(hc) / 10), ("KCU1", tc, H(h1) / 10),
        ("KCU0 fragments only", tc, H(TT['kcu0_frag']) / 10)], 0, 6800, (0.2, 0.67, 0.93, 0.86))
# 2 overlaps
r = lambda a, b: np.where(H(a) > 0, H(a & b) / np.maximum(H(a), 1), np.nan)
graphs(2, "Share of a stream's triggers also recorded by another;run time [s];share", [("KCU1 also in KCU0", tc, r(h1, h0)), ("CAEN also in KCU0", tc, r(hc, h0)),
       ("KCU1 also in CAEN", tc, r(h1, hc)), ("KCU0 also in CAEN", tc, r(h0, hc))], 0, 1.6, (0.2, 0.7, 0.93, 0.86))
# 3 KCU0-KCU1 residual
A1 = [uproot.open(f)['kcu'].arrays(['ts', 'nsamp'], library='np') for f in sorted(glob.glob('out/beam_ROC_1-00003062-*.root'))]
tt1 = np.concatenate([a['ts'][a['nsamp'] == 20] for a in A1]).astype(np.int64); tt1 = np.unique(tt1[tt1 > 1e8])
b01 = h0 & h1; d10 = T[b01] - (tt1[i1[b01]] + 1143665946); tb = ts[b01]
hb = lambda m: np.histogram(tb[m], E)[0].astype(float); nb = np.maximum(hb(np.ones(len(tb), bool)), 1)
graphs(3, "KCU0 - KCU1 residual after the fixed offset;run time [s];share of shared triggers",
       [("0 ticks", tc, hb(d10 == 0) / nb), ("-4 ticks (KCU1 +25 ns)", tc, hb(d10 == -4) / nb), ("other", tc, hb(~np.isin(d10, [0, -4])) / nb)], 1e-4, 2, (0.55, 0.35, 0.9, 0.6), logy=True)
print("KCU0-KCU1: 0 %.4f  -4 %.4f  other %.5f over %d shared triggers" % ((d10 == 0).mean(), (d10 == -4).mean(), (~np.isin(d10, [0, -4])).mean(), len(d10)))
# 4 CAEN-KCU0 residual vs time
M = np.load('out/caen_kcu0_match.npz'); tp = M['tpred']
b0c = h0 & hc; rc = (tp[ic[b0c]] / 6.25 - T[b0c]) * 6.25; rc -= np.median(rc)
cv.cd(4); colz(); h4 = ROOT.TH2D("h4", "CAEN - KCU0 residual vs time;run time [s];residual [ns]", 150, 0, 1500, 60, -30, 30); keep.append(h4)
for x, yv in zip(ts[b0c][::5], rc[::5]): h4.Fill(x, yv)
h4.Draw("colz")
med = [np.median(rc[(ts[b0c] >= a) & (ts[b0c] < a + 50)]) if ((ts[b0c] >= a) & (ts[b0c] < a + 50)).sum() > 50 else np.nan for a in np.arange(0, 1500, 50)]
g = ROOT.TGraph(30, np.arange(25, 1500, 50.0), np.array(med, float)); g.SetLineColor(ROOT.kRed); g.SetLineWidth(3); g.Draw("L same"); keep.append(g)
# 5 residual distributions
b1c = h1 & hc & ~h0; r1c = (tp[ic[b1c]] / 6.25 - (tt1[i1[b1c]] + 1143665946)) * 6.25; r1c -= np.median(r1c)
cv.cd(5); ROOT.gPad.SetGrid(); ROOT.gPad.SetLogy()
ha = ROOT.TH1D("ha", "Trigger-time residuals;residual [ns];triggers", 120, -60, 60); hb_ = ROOT.TH1D("hb", "", 120, -60, 60); keep += [ha, hb_]
for v in rc: ha.Fill(v)
for v in r1c: hb_.Fill(v)
ha.SetLineColor(ROOT.kBlack); ha.SetLineWidth(2); hb_.SetLineColor(ROOT.kRed + 1); hb_.SetLineWidth(2); ha.Draw("hist"); hb_.Draw("hist same")
mad = lambda x: np.median(np.abs(x - np.median(x)))
ha.SetMaximum(3e7); lg = legend(0.2, 0.72, 0.93, 0.86, 0.04); lg.AddEntry(ha, "CAEN - KCU0, MAD %.1f ns (%d)" % (mad(rc), len(rc)), "l")
lg.AddEntry(hb_, "CAEN - KCU1, no KCU0, MAD %.1f ns (%d)" % (mad(r1c), len(r1c)), "l"); lg.Draw(); keep.append(lg)
# 6 paddle agreement
p0, pc = TT['pad_kcu0'], TT['pad_caen']; sh = h0 & hc; same = (p0 == pc); dec = np.zeros(len(T), bool)
ish = np.nonzero(sh)[0]; dec[ish] = p0[ish] == np.roll(pc[ish], 1)
hs = lambda m: np.histogram(ts[m], E)[0].astype(float); ns_ = np.maximum(hs(sh), 1)
graphs(6, "Paddle pattern, KCU0 vs CAEN, on shared triggers;run time [s];identical pattern", [("same trigger", tc, hs(sh & same) / ns_),
       ("neighbour trigger (decoy)", tc, hs(sh & dec) / ns_)], 0, 1.05, (0.5, 0.4, 0.9, 0.55))
print("paddle agreement %.4f, decoy %.4f" % (same[sh].mean(), dec[sh].mean()))
# 7 AstroPix phase lock vs time
R = np.load('out/astro_matcher_rows.npz'); ka, band = R['ki'], R['band']; tA = ts[ka]
c = np.load('out/astro_clusters.npz'); Fm = c['Fmin']; O = int(c['coarse']) - int(0.43e9 / 6.25)
EA = np.arange(0, 1500, 50.0); tca = 0.5 * (EA[1:] + EA[:-1])
ncl = np.histogram((Fm + O - t0) * 6.25e-9, EA)[0].astype(float)
sig = np.histogram(tA[band == 0], EA)[0]; ctl = np.histogram(tA[band > 0], EA)[0] / 2.0
graphs(7, "AstroPix phase lock: (signal - control) per cluster;run time [s];excess candidates / cluster",
       [("phase 189 minus control bands", tca, (sig - ctl) / np.maximum(ncl, 1))], 0, 0.15, (0.4, 0.15, 0.9, 0.3))
# 8 latency of good matches vs time
Mx = np.load('out/astro_kcu0_matches.npz'); gm = Mx['log10LR'] > 2.5
cv.cd(8); colz()
h8 = ROOT.TH2D("h8", "AstroPix readout latency, log_{10}LR > 2.5;run time [s];trigger - AstroPix stamp [ms]", 75, 0, 1500, 60, -110, 10); keep.append(h8)
for x, yv in zip(sec(Mx['trig_ts'][gm]), -Mx['latency_ms'][gm]): h8.Fill(x, yv)
h8.Draw("colz")
# 9 AstroPix match rate
nm = np.histogram(sec(Mx['trig_ts'][gm]), E)[0].astype(float); nc = np.histogram((Fm + O - t0) * 6.25e-9, E)[0].astype(float)
graphs(9, "AstroPix clusters and matches;run time [s];per 10 s", [("clusters / 100", tc, nc / 100), ("matches, log_{10}LR > 2.5", tc, nm)], 0, 900, (0.2, 0.72, 0.8, 0.86))
finish(cv); cv.SaveAs("out/sync_diag.png")
