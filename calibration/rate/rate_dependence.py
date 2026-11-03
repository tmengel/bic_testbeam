"""Rate and pile-up dependence of the SFIL response in run 3062.  Per shower (fit cache, KCU0 South layer sums), the time since
the previous trigger (dt_prev, from the KCU0 header timestamps, 6.25 ns ticks), the time to the next trigger, and the
trigger rate in a +-50 ms window; the mean layer sums (normalised per paddle to remove the energy ordering) against them.
    python3 rate_dependence.py -> rate_dependence.txt/.png
"""
import numpy as np, uproot, sys, ROOT
sys.path.insert(0, '/sphenix/user/tmengel/eic/bic_testbeam/common'); from plotstyle import setup, legend, label
setup(); ROOT.gROOT.SetBatch(True); out = open('rate_dependence.txt', 'w')
def say(*a): print(*a, flush=True); print(*a, file=out)
D = np.load('/sphenix/user/tmengel/eic/bic_testbeam/sim/analysis/data_fit_layers.npz'); S = D['S'].astype(float); pad = D['pad']; trig = D['trig']
ts_all = uproot.open('/sphenix/user/tmengel/eic/bic_testbeam/data/kcu0_3062.root')['header']['ts0'].array(library='np').astype(np.int64)
t = ts_all * 6.25e-9                                         # s, all KCU0 triggers (issued triggers with a complete waveform in the file)
order = np.argsort(t); ts = t[order]; rank = np.empty_like(order); rank[order] = np.arange(len(t))
dprev = np.full(len(t), np.nan); dnext = np.full(len(t), np.nan); dprev[order[1:]] = np.diff(ts); dnext[order[:-1]] = np.diff(ts)
# rate in +-50 ms: count of triggers in the window / 0.1 s
lo = np.searchsorted(ts, ts - 0.05); hi = np.searchsorted(ts, ts + 0.05); rate_sorted = (hi - lo) / 0.1; rate = np.empty(len(t)); rate[order] = rate_sorted
dp, dn, rt = dprev[trig], dnext[trig], rate[trig]; tot = S.sum(axis=1)
shower = (tot > 600) & np.isfinite(dp)
say('showers %d; dt_prev percentiles 1/10/50/90 %%: %s ms; rate percentiles: %s Hz' % (shower.sum(), np.round(1e3 * np.percentile(dp[shower], [1, 10, 50, 90]), 3), np.round(np.percentile(rt[shower], [10, 50, 90]), 0)))
say('note: dt_prev is measured between recorded triggers; the 100 us busy means no dt below 0.1 ms, and UDP loss removes some neighbours')
# normalise each shower's layer sums by the paddle mean (removes the energy ordering); then profiles
norm = np.ones(len(S))
for p in (3, 4, 5, 6):
    m = shower & (pad == p); norm[m] = tot[m].mean()
r = tot / norm
edges_dt = np.array([0.1, 0.15, 0.2, 0.3, 0.5, 0.7, 1, 1.5, 2, 3, 5, 10, 30]) * 1e-3
say('\n<total SFIL1-4 South> / paddle mean against the time since the previous trigger:')
xs, ys, es = [], [], []
for a, b in zip(edges_dt[:-1], edges_dt[1:]):
    m = shower & (dp >= a) & (dp < b)
    if m.sum() < 200: continue
    xs.append(np.sqrt(a * b) * 1e3); ys.append(r[m].mean()); es.append(r[m].std() / np.sqrt(m.sum())); say('  dt_prev %6.2f-%6.2f ms: n %7d  ratio %.4f +- %.4f   per layer %s' % (a * 1e3, b * 1e3, m.sum(), ys[-1], es[-1], np.round((S[m] / norm[m, None]).mean(axis=0) / (S[shower] / norm[shower, None]).mean(axis=0), 3)))
say('\n... against the time to the NEXT trigger (a check: no causal effect expected):')
xs2, ys2 = [], []
for a, b in zip(edges_dt[:-1], edges_dt[1:]):
    m = shower & (dn >= a) & (dn < b)
    if m.sum() < 200: continue
    xs2.append(np.sqrt(a * b) * 1e3); ys2.append(r[m].mean()); say('  dt_next %6.2f-%6.2f ms: n %7d  ratio %.4f' % (a * 1e3, b * 1e3, m.sum(), ys2[-1]))
say('\n... against the trigger rate in a +-50 ms window:')
edges_r = np.array([0, 1000, 2000, 3000, 4000, 5000, 6000, 8000, 12000])
xr, yr, er = [], [], []
for a, b in zip(edges_r[:-1], edges_r[1:]):
    m = shower & (rt >= a) & (rt < b)
    if m.sum() < 200: continue
    xr.append(0.5 * (a + b)); yr.append(r[m].mean()); er.append(r[m].std() / np.sqrt(m.sum())); say('  rate %5d-%5d Hz: n %7d  ratio %.4f +- %.4f  per layer %s' % (a, b, m.sum(), yr[-1], er[-1], np.round((S[m] / norm[m, None]).mean(axis=0) / (S[shower] / norm[shower, None]).mean(axis=0), 3)))
# amplitude of the shower in the previous trigger: pile-up from a large preceding shower (tails) - need previous trigger's sum: map via rank
prev_idx = np.full(len(t), -1); prev_idx[order[1:]] = order[:-1]
tot_all = np.full(len(t), np.nan); tot_all[trig] = tot; prev_tot = tot_all[prev_idx[trig]]
say('\nratio against the previous trigger\'s own SFIL sum (dt_prev < 0.5 ms):')
for a, b in ((0, 100), (100, 600), (600, 2000), (2000, 4000), (4000, 20000)):
    m = shower & (dp < 0.5e-3) & (prev_tot >= a) & (prev_tot < b)
    if m.sum() > 200: say('  previous sum %5d-%5d: n %6d  ratio %.4f' % (a, b, m.sum(), r[m].mean()))
# figure
c = ROOT.TCanvas('c', '', 2000, 700); c.Divide(3, 1); keep = []
c.cd(1); ROOT.gPad.SetLogx(); fr = ROOT.gPad.DrawFrame(0.08, 0.9, 40, 1.1, ';time since the previous trigger [ms];<SFIL1-4 South> / paddle mean'); keep.append(fr)
g = ROOT.TGraphErrors(len(xs), np.array(xs), np.array(ys), np.zeros(len(xs)), np.array(es)); g.SetMarkerStyle(20); g.Draw('p same'); keep.append(g)
g2 = ROOT.TGraph(len(xs2), np.array(xs2), np.array(ys2)); g2.SetMarkerStyle(24); g2.SetMarkerColor(ROOT.kBlue + 1); g2.Draw('p same'); keep.append(g2)
lg = legend(0.45, 0.7, 0.88, 0.88); lg.AddEntry(g, 'previous trigger', 'p'); lg.AddEntry(g2, 'next trigger (control)', 'p'); lg.Draw(); keep.append(lg)
c.cd(2); fr2 = ROOT.gPad.DrawFrame(0, 0.9, 12000, 1.1, ';trigger rate in #pm50 ms [Hz];<SFIL1-4 South> / paddle mean'); keep.append(fr2)
g3 = ROOT.TGraphErrors(len(xr), np.array(xr), np.array(yr), np.zeros(len(xr)), np.array(er)); g3.SetMarkerStyle(20); g3.Draw('p same'); keep.append(g3)
c.cd(3); ROOT.gPad.SetLogy(); h = ROOT.TH1F('hdt', ';time since the previous trigger [ms];showers', 100, 0, 5); [h.Fill(v * 1e3) for v in dp[shower]]; h.Draw(); keep.append(h)
c.SaveAs('rate_dependence.png'); c.SaveAs('rate_dependence.pdf'); say('wrote rate_dependence.png')
