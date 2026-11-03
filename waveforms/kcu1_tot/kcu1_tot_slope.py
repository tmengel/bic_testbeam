"""TOT slope for the SFIL5 South single SiPMs.  Below the ceiling A_fit is used; a clipped pulse carries only its time
over the ceiling (decoded TOT, ~400 units per 25 ns sample; nsat gives the same at 25 ns resolution).  The HGCROC TOT is
linear in charge above its toe, so A(TOT) = A_c + b (TOT - T_c) with the ceiling crossing T_c and slope b unknown.
b is chosen so that <A5 | A4> is as linear as the model's <L5 | L4> (adjacent layers, cells 1-3), using triggers whose
SFIL4 South channels are unsaturated (unbiased A4 bins).  The absolute SFIL5 scale is not touched here.
    python3 kcu1_tot_slope.py  -> kcu1_tot_slope.txt/.png, and sfil5_recovered.npz rewritten with the chosen b
"""
import numpy as np, uproot, glob, sys, ROOT
sys.path.insert(0, '/sphenix/user/tmengel/eic/bic_testbeam/common'); from plotstyle import setup, legend, label
setup(); out = open('kcu1_tot_slope.txt', 'w')
def say(*a): print(*a, flush=True); print(*a, file=out)
D = np.load('/sphenix/tg/tg01/jets/tmengel/eic/kcu1/fit_full/kcu1_fit_arrays.npz')
afit, amp, nsat, tot, ok, ts0 = D['afit'], D['amp'], D['nsat'], D['tot'], D['ok'], D['ts0']; n = len(ts0)
live = ok.sum(axis=0) > 20000; ks = [k for k in range(16, 64) if live[k]]
unsat = ok & (nsat == 0); sat = nsat > 0; hastot = (tot > 0) & (tot < 4092)
# ---- TOT <-> nsat: median TOT per nsat (for clipped pulses without a TOT word) and the ceiling crossing T_c
say('clipped pulses (live B+C+D channels): fraction with a TOT word by nsat, and median TOT')
Tn = {}
for s in range(1, 12):
    m = sat[:, ks] & (nsat[:, ks] == s); mt = m & hastot[:, ks]
    if m.sum() < 50: continue
    Tn[s] = np.median(tot[:, ks][mt]) if mt.sum() > 20 else np.nan
    say('  nsat %2d: n %7d  with TOT %.2f  median TOT %5.0f  16-84 %% %5.0f-%5.0f' % (s, m.sum(), mt.sum() / m.sum(), Tn[s], np.percentile(tot[:, ks][mt], 16), np.percentile(tot[:, ks][mt], 84)))
T1 = Tn[1]; T2 = Tn[2]; Tc = T1 - 0.5 * (T2 - T1)       # TOT at which the pulse just touches the ceiling: half a sample below the nsat = 1 median
say('ceiling crossing T_c = %.0f (nsat 1 median %.0f, per-sample step %.0f units)' % (Tc, T1, T2 - T1))
CEIL = 943.
def rebuild(b):
    A = np.zeros((n, len(ks)), np.float32)
    for i, k in enumerate(ks):
        f = unsat[:, k]; A[f, i] = afit[f, k]
        r = ~ok[:, k] & (nsat[:, k] == 0) & (amp[:, k] > 60); A[r, i] = amp[r, k]
        s = sat[:, k]; T = np.where(hastot[:, k], tot[:, k], 0).astype(float)
        for ns_, tv in Tn.items(): T[s & ~hastot[:, k] & (nsat[:, k] == ns_)] = tv
        T[s & (T == 0)] = Tn[max(Tn)]
        A[s, i] = CEIL + b * np.maximum(T[s] - Tc, 0)
    return A.sum(axis=1)
# ---- join: KCU1 -> events rows; KCU0 fit layers -> events rows
ev = uproot.open('/sphenix/user/tmengel/eic/bic_testbeam/data/events_3062.root')['events'].arrays(['trig_ts', 'pad'], library='np')
ets = ev['trig_ts']; order = np.argsort(ets); ets_s = ets[order]
def tojoin(t):
    pos = np.clip(np.searchsorted(ets_s, t), 0, len(ets_s) - 1); hit = ets_s[pos] == t; r = np.full(len(t), -1, np.int64); r[hit] = order[pos[hit]]; return r
row1 = tojoin(ts0 + 1143665946)
F = np.load('/sphenix/user/tmengel/eic/bic_testbeam/sim/analysis/data_fit_layers.npz'); S4 = F['S'][:, 3].astype(float); sat4 = F['satS'][:, 3] > 0; pad = F['pad']
ts0k = uproot.open('/sphenix/user/tmengel/eic/bic_testbeam/data/kcu0_3062.root')['header']['ts0'].array(library='np'); row0 = tojoin(ts0k[F['trig']])
A4 = np.full(len(ets), np.nan); A4[row0[row0 >= 0]] = S4[row0 >= 0]; SAT4 = np.zeros(len(ets), bool); SAT4[row0[row0 >= 0]] = sat4[row0 >= 0]; PAD = np.zeros(len(ets), np.int16); PAD[row0[row0 >= 0]] = pad[row0 >= 0]
h1 = row1 >= 0; a4 = A4[row1[h1]]; s4 = SAT4[row1[h1]]; pd = PAD[row1[h1]]; good = np.isfinite(a4) & ~s4 & (a4 > 100)
say('KCU1 triggers joined %d; with KCU0 fit layers %d; SFIL4 unsaturated and A4 > 100: %d' % (h1.sum(), np.isfinite(a4).sum(), good.sum()))
# ---- model <L5 | L4>
sim = {}
for f in sorted(glob.glob('/sphenix/tg/tg01/jets/tmengel/eic/sim/bic_sim_e*_s*.root')):
    g = uproot.open(f)['sim'].arrays(['ebeam', 'sfil_ls'], library='np'); E = int(round(g['ebeam'].mean()))
    if 2500 <= E <= 6500 and '_x80' not in f: sim.setdefault(E, []).append(g['sfil_ls'][:, 3:5, :3].sum(axis=2))
sim = {E: np.concatenate(v) for E, v in sim.items()}
say('\nmodel <L5|L4> in L4 quantile bins (cells 1-3, South light, beam under P6), ratio to a straight line through the origin fitted to the bins:')
eb = np.array([0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 0.97])
mshape = {}
for E in sorted(sim):
    L4, L5 = sim[E][:, 0], sim[E][:, 1]; q = np.quantile(L4, eb); xs, ys = [], []
    for lo, hi in zip(q[:-1], q[1:]):
        m = (L4 >= lo) & (L4 < hi); xs.append(L4[m].mean()); ys.append(L5[m].mean())
    xs, ys = np.array(xs), np.array(ys); p = np.polyfit(xs, ys, 1); mshape[E] = (xs, ys, p)
    say('  %d GeV: slope %.3f intercept %.1f (= %.0f %% of <L5>); <L5>/<L4> %.3f; bin ratios L5/L4: %s' % (E // 1000, p[0], p[1], 100 * p[1] / L5.mean(), L5.mean() / L4.mean(), ' '.join('%.2f' % (y / x) for x, y in zip(xs, ys))))
# ---- data <A5(b) | A4> per b, in A4 bins, unsaturated SFIL4
bins = np.array([100, 200, 300, 400, 500, 600, 800, 1000, 1200, 1500, 2000, 3000])
res = {}
say('\ndata <A5(b)> / <A4> in A4 bins (SFIL4 South unsaturated), normalised to the 300-400 bin; the model (4 GeV) bin ratios normalised the same way are the target shape')
xs4, ys4, _ = mshape[4000]; say('  model 4 GeV L5/L4 normalised: %s (L4 quantile bins)' % ' '.join('%.2f' % ((y / x) / (ys4[2] / xs4[2])) for x, y in zip(xs4, ys4)))
for b in (0.0, 0.25, 0.5, 0.75, 1.0, 1.5, 2.0, 3.0):
    A5 = rebuild(b)[h1]; rat = []; cnt = []
    for lo, hi in zip(bins[:-1], bins[1:]):
        m = good & (a4 >= lo) & (a4 < hi); cnt.append(m.sum()); rat.append(A5[m].mean() / a4[m].mean() if m.sum() > 50 else np.nan)
    rat = np.array(rat); res[b] = (rat, cnt)
    say('  b = %.2f: %s   (ratio 1500-3000 / 300-400 = %.2f)' % (b, ' '.join('%5.2f' % (r / rat[2]) for r in rat), rat[-1] / rat[2]))
say('  A4 bin counts: %s' % cnt)
# choose b: the one whose high/low ratio is closest to the model's (4 GeV, last/third bin)
target = (ys4[-1] / xs4[-1]) / (ys4[2] / xs4[2]); bs = sorted(res); vals = [res[b][0][-1] / res[b][0][2] for b in bs]
bbest = float(np.interp(target, vals, bs)) if np.all(np.diff(vals) > 0) else bs[int(np.argmin(np.abs(np.array(vals) - target)))]
say('\nmodel target for (1500-3000)/(300-400): %.2f  ->  b = %.2f ADC per TOT unit (%.0f ADC per clipped 25 ns sample)' % (target, bbest, bbest * (T2 - T1)))
A5 = rebuild(bbest); A5flat = rebuild(0.0)
np.savez('/sphenix/user/tmengel/eic/bic_testbeam/sim/analysis/sfil5_recovered.npz', row=row1[h1], A5=A5[h1], A5old=A5flat[h1], nsat5=sat[:, ks].sum(axis=1)[h1], b=bbest, Tc=Tc, step=T2 - T1)
say('wrote sfil5_recovered.npz with b = %.2f (A5old field = ceiling-only, b = 0); mean A5: b = 0 %.0f, chosen %.0f, per-trigger ratio for >= 4 clipped channels %.2f' % (bbest, A5flat.mean(), A5.mean(), A5[sat[:, ks].sum(axis=1) >= 4].mean() / A5flat[sat[:, ks].sum(axis=1) >= 4].mean()))
# ---- figure
c = ROOT.TCanvas('c', '', 2000, 800); c.Divide(2, 1); keep = []
c.cd(1); fr = ROOT.gPad.DrawFrame(0, 0, 3000, 2.0, ';A4 (SFIL4 South, fitted sum, unsaturated) [KCU0 ADC];<A5>/<A4> normalised to A4 = 300-400'); keep.append(fr)
cols = [ROOT.kGray + 1, ROOT.kBlue + 1, ROOT.kAzure + 7, ROOT.kGreen + 2, ROOT.kOrange + 7, ROOT.kRed + 1, ROOT.kMagenta + 1, ROOT.kBlack]
lg = legend(0.15, 0.6, 0.5, 0.88); xc = 0.5 * (bins[1:] + bins[:-1])
for b, colr in zip(bs, cols):
    rat = res[b][0] / res[b][2]; f = np.isfinite(rat); g = ROOT.TGraph(int(f.sum()), xc[f].astype(float), rat[f]); g.SetLineColor(colr); g.SetMarkerColor(colr); g.SetMarkerStyle(20); g.SetLineWidth(2); g.Draw('lp same'); keep.append(g); lg.AddEntry(g, 'b = %.2f' % b, 'lp')
lg.Draw(); keep.append(lg); keep.append(label('data, SFIL5 rebuilt with A = 943 + b (TOT - T_{c}) for clipped pulses', x=0.15, y=0.93, size=0.03))
c.cd(2); fr = ROOT.gPad.DrawFrame(0, 0, 1.0, 2.0, ';L4 quantile (model, cells 1-3);<L5>/<L4> normalised to the third bin'); keep.append(fr)
lg2 = legend(0.15, 0.65, 0.45, 0.88)
for E, colr in zip(sorted(mshape), (ROOT.kBlue + 1, ROOT.kGreen + 2, ROOT.kOrange + 7, ROOT.kRed + 1)):
    xs, ys, p = mshape[E]; r = (ys / xs) / (ys[2] / xs[2]); g = ROOT.TGraph(len(r), 0.5 * (eb[1:] + eb[:-1]), r); g.SetLineColor(colr); g.SetMarkerColor(colr); g.SetMarkerStyle(20); g.SetLineWidth(2); g.Draw('lp same'); keep.append(g); lg2.AddEntry(g, '%d GeV' % (E // 1000), 'lp')
lg2.Draw(); keep.append(lg2); keep.append(label('model: adjacent-layer relation', x=0.15, y=0.93, size=0.03))
c.SaveAs('kcu1_tot_slope.png'); c.SaveAs('kcu1_tot_slope.pdf'); say('wrote kcu1_tot_slope.png')
