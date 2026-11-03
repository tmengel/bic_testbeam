"""TOT slope b for clipped SFIL5 South pulses from the WIDTH of the per-paddle SFIL5 distribution: clipping compresses the
upper tail of A5, and A = 943 + b (TOT - T_c) restores it.  The model gives the relative width of L5 (cells 1-3) at
3-6 GeV; b is read off where the data quantile ratios q90/q50 and q97/q50 match.  Model-shape anchored (not scale).
    python3 kcu1_tot_shape.py  -> kcu1_tot_shape.txt/.png, sfil5_recovered.npz (A5 at the central b, A5lo/A5hi at the range)
"""
import numpy as np, uproot, glob, sys, ROOT
sys.path.insert(0, '/sphenix/user/tmengel/eic/bic_testbeam/common'); from plotstyle import setup, legend, label
setup(); out = open('kcu1_tot_shape.txt', 'w')
def say(*a): print(*a, flush=True); print(*a, file=out)
D = np.load('/sphenix/tg/tg01/jets/tmengel/eic/kcu1/fit_full/kcu1_fit_arrays.npz')
afit, amp, nsat, tot, ok, ts0 = D['afit'], D['amp'], D['nsat'], D['tot'], D['ok'], D['ts0']; n = len(ts0)
live = ok.sum(axis=0) > 20000; ks = [k for k in range(16, 64) if live[k]]
unsat = ok & (nsat == 0); sat = nsat > 0; hastot = (tot > 0) & (tot < 4092)
Tn = {}
for s in range(1, 12):
    mt = sat[:, ks] & (nsat[:, ks] == s) & hastot[:, ks]
    if mt.sum() > 20: Tn[s] = np.median(tot[:, ks][mt])
Tc = Tn[1] - 0.5 * (Tn[2] - Tn[1]); CEIL = 943.
base = np.zeros(n, np.float32); T = np.zeros((n, len(ks)), np.float32)
for i, k in enumerate(ks):
    f = unsat[:, k]; base += np.where(f, afit[:, k], 0); r = ~ok[:, k] & (nsat[:, k] == 0) & (amp[:, k] > 60); base += np.where(r, amp[:, k], 0)
    s = sat[:, k]; t = np.where(hastot[:, k], tot[:, k], 0).astype(np.float32)
    for ns_, tv in Tn.items(): t[s & ~hastot[:, k] & (nsat[:, k] == ns_)] = tv
    t[s & (t == 0)] = Tn[max(Tn)]; T[:, i] = np.where(s, np.maximum(t - Tc, 0), 0)
ncl = sat[:, ks].sum(axis=1); excess = T.sum(axis=1)                     # A5(b) = base + 943 ncl + b * excess
def A5(b): return base + CEIL * ncl + b * excess
ev = uproot.open('/sphenix/user/tmengel/eic/bic_testbeam/data/events_3062.root')['events'].arrays(['trig_ts'], library='np'); ets = ev['trig_ts']; order = np.argsort(ets); ets_s = ets[order]
def tojoin(t):
    pos = np.clip(np.searchsorted(ets_s, t), 0, len(ets_s) - 1); hit = ets_s[pos] == t; r = np.full(len(t), -1, np.int64); r[hit] = order[pos[hit]]; return r
row1 = tojoin(ts0 + 1143665946)
F = np.load('/sphenix/user/tmengel/eic/bic_testbeam/sim/analysis/data_fit_layers.npz'); Stot = F['S'].sum(axis=1).astype(float); pad = F['pad']
ts0k = uproot.open('/sphenix/user/tmengel/eic/bic_testbeam/data/kcu0_3062.root')['header']['ts0'].array(library='np'); row0 = tojoin(ts0k[F['trig']])
N = len(ets); ST = np.full(N, np.nan); ST[row0[row0 >= 0]] = Stot[row0 >= 0]; PAD = np.zeros(N, np.int16); PAD[row0[row0 >= 0]] = pad[row0 >= 0]
h1 = row1 >= 0; st = ST[row1[h1]]; pd = PAD[row1[h1]]; shower = np.isfinite(st) & (st > 600)
say('KCU1 triggers joined %d, with KCU0 layers and SFIL1-4 sum > 600: %d; T_c = %.0f, per-sample TOT step %.0f' % (h1.sum(), shower.sum(), Tc, Tn[2] - Tn[1]))
# ---- model L5 quantile ratios
sim = {}
for fn in sorted(glob.glob('/sphenix/tg/tg01/jets/tmengel/eic/sim/bic_sim_e*_s*.root')):
    g = uproot.open(fn)['sim'].arrays(['ebeam', 'sfil_ls'], library='np'); E = int(round(g['ebeam'].mean())); x = 80 if '_x80' in fn else 30
    if 2500 <= E <= 6500: sim.setdefault((x, E), []).append(g['sfil_ls'][:, 4, :3].sum(axis=1))
sim = {k: np.concatenate(v) for k, v in sim.items()}
qs = (50, 90, 97); 
def qr(v): q = np.percentile(v, qs); return q[1] / q[0], q[2] / q[0], q[0]
say('\nmodel L5 (cells 1-3, South light) quantile ratios q90/q50, q97/q50:')
for (x, E) in sorted(sim): r = qr(sim[(x, E)]); say('  x = %d mm, %d GeV: %.3f  %.3f  (n %d)' % (x, E // 1000, r[0], r[1], len(sim[(x, E)])))
# ---- data per b
bs = np.array([0, 0.25, 0.5, 0.75, 1.0, 1.25, 1.5, 2.0, 2.5, 3.0, 4.0])
PB = {4: 4, 3: 3, 6: 6, 5: 5}; res = {}
say('\ndata A5(b) quantile ratios per paddle (showers):')
for p in (4, 3, 6, 5):
    m = shower & (pd == p); rows = []
    for b in bs: rows.append(qr(A5(b)[h1][m]))
    res[p] = np.array(rows); say('  P%d (n %d): ' % (p, m.sum()) + '  '.join('b %.2f: %.3f/%.3f' % (b, r[0], r[1]) for b, r in zip(bs, rows)))
# ---- read off b per paddle and energy hypothesis, with the model L5 smeared by a per-event beam-energy spread sigma_E/E
# (L5 scales roughly as E^0.9 between the simulated energies); the paddle energy order is P4 < P3 < P6 < P5 within 3-6 GeV
rng = np.random.default_rng(1)
def model_q(x, E, spread):
    v = sim[(x, E)]
    if spread > 0: v = v * np.clip(1 + spread * rng.standard_normal(len(v)), 0.2, None) ** 0.9
    return qr(v)
say('\nb (ADC per TOT unit) from q90/q50 | q97/q50 per paddle, for energy hypotheses and model energy spreads:')
hyp = {4: (3000, 4000), 3: (3000, 4000, 5000), 6: (4000, 5000), 5: (5000, 6000)}
table = {}
for spread in (0.0, 0.1, 0.2, 0.3):
    vals = []
    for p in (4, 3, 6, 5):
        x = 80 if p == 5 else 30; line = []
        for E in hyp[p]:
            if (x, E) not in sim: continue
            mr = model_q(x, E, spread)
            b90 = np.interp(mr[0], res[p][:, 0], bs) if mr[0] <= res[p][-1, 0] else np.nan; b97 = np.interp(mr[1], res[p][:, 1], bs) if mr[1] <= res[p][-1, 1] else np.nan
            if mr[0] < res[p][0, 0]: b90 = 0.
            if mr[1] < res[p][0, 1]: b97 = 0.
            line.append('%d GeV %.2f|%.2f' % (E // 1000, b90, b97)); vals += [b90, b97]
        say('  spread %2.0f %%  P%d: ' % (100 * spread, p) + '   '.join(line))
    table[spread] = np.nanmedian(vals); say('    -> median b = %.2f' % table[spread])
bcen = table[0.2]; blo = table[0.0]; bhi = table[0.3]
say('\ncentral b = %.2f (20 %% energy spread), band %.2f (no spread) to %.2f (30 %% spread); the toe slope below the ceiling is 0.33' % (bcen, blo, bhi))
for b in (0, blo, bcen, bhi, 1.0):
    a = A5(b)[h1][shower]; say('  b = %.2f: mean A5 in showers %.0f, mean over triggers with >= 4 clipped channels %.0f' % (b, a.mean(), A5(b)[h1][shower & (ncl[h1] >= 4)].mean()))
np.savez('/sphenix/user/tmengel/eic/bic_testbeam/sim/analysis/sfil5_recovered.npz', row=row1[h1], A5=A5(bcen)[h1], A5lo=A5(blo)[h1], A5hi=A5(bhi)[h1], A5old=A5(0)[h1], nsat5=ncl[h1], b=bcen, blo=blo, bhi=bhi, Tc=Tc, step=Tn[2] - Tn[1])
say('wrote sfil5_recovered.npz (A5 at b = %.2f; A5lo/A5hi at %.2f/%.2f; A5old = ceiling only)' % (bcen, blo, bhi))
# ---- figure: q90/q50 vs b per paddle with model bands; A5 distributions at b = 0 and central for P6 vs model 4-5 GeV (scaled to the median)
c = ROOT.TCanvas('c', '', 2000, 800); c.Divide(2, 1); keep = []
c.cd(1); fr = ROOT.gPad.DrawFrame(0, 1.0, 4.2, 2.4, ';b  [ADC per TOT unit];q_{90} / q_{50} of the SFIL5 South sum'); keep.append(fr)
lg = legend(0.15, 0.62, 0.5, 0.88)
for p, colr in zip((4, 3, 6, 5), (ROOT.kBlue + 1, ROOT.kGreen + 2, ROOT.kOrange + 7, ROOT.kRed + 1)):
    g = ROOT.TGraph(len(bs), bs.astype(float), res[p][:, 0].astype(float)); g.SetLineColor(colr); g.SetMarkerColor(colr); g.SetMarkerStyle(20); g.SetLineWidth(2); g.Draw('lp same'); keep.append(g); lg.AddEntry(g, 'data P%d' % p, 'lp')
for E, sty in zip((3000, 4000, 5000, 6000), (2, 3, 4, 5)):
    r = qr(sim[(30, E)]); ln = ROOT.TLine(0, r[0], 4.2, r[0]); ln.SetLineStyle(sty); ln.SetLineWidth(2); ln.Draw(); keep.append(ln); lg.AddEntry(ln, 'model %d GeV' % (E // 1000), 'l')
lg.Draw(); keep.append(lg)
c.cd(2); m = shower & (pd == 6); fr = ROOT.gPad.DrawFrame(0, 0, 4, 0.12, ';A5 / median;fraction per bin'); keep.append(fr); lg2 = legend(0.5, 0.65, 0.88, 0.88)
for b, colr, nm in ((0, ROOT.kGray + 2, 'data P6, b = 0 (ceiling only)'), (bcen, ROOT.kRed + 1, 'data P6, b = %.2f' % bcen), (bhi, ROOT.kOrange + 7, 'data P6, b = %.2f' % bhi)):
    a = A5(b)[h1][m]; h = ROOT.TH1F('h%.2f' % b, '', 80, 0, 4); [h.Fill(x) for x in a / np.median(a)]; h.Scale(1. / h.GetEntries()); h.SetLineColor(colr); h.SetLineWidth(2); h.Draw('hist same'); keep.append(h); lg2.AddEntry(h, nm, 'l')
for E, colr in ((4000, ROOT.kBlue + 1), (5000, ROOT.kGreen + 2)):
    v = sim[(30, E)]; h = ROOT.TH1F('hm%d' % E, '', 80, 0, 4); [h.Fill(x) for x in v / np.median(v)]; h.Scale(1. / h.GetEntries()); h.SetLineColor(colr); h.SetLineWidth(2); h.SetLineStyle(2); h.Draw('hist same'); keep.append(h); lg2.AddEntry(h, 'model L5, %d GeV' % (E // 1000), 'l')
lg2.Draw(); keep.append(lg2)
c.SaveAs('kcu1_tot_shape.png'); c.SaveAs('kcu1_tot_shape.pdf'); say('wrote kcu1_tot_shape.png')
