"""Recreate the 2017 dual-sided sPHENIX EmCal study (Craig, Sean, Bob, Martin, FTBF, DRS4 at 2 GS/s) with run 3062.
Their module position scan -> our four paddles (x = -80, -30, +30, +80 mm) and the AstroPix x; their two SiPMs -> the North
and South ends of the beam cell (cell 3) of SFIL2 and SFIL3; their cuts -> one paddle, neither end clipped (raw max < 1000).
Positive x is towards the South end, so South plays their "right" SiPM:  r = (S_S - S_N) / (S_S + S_N).
Gain: the South amplitude is scaled so the mean asymmetry is 0 at x = 0 (average of P3 and P6), as they did at position 0.
Timing per end, from the same 25 ns waveforms:
  CFD      c[n] = w[n] - f w[n+d], f = 0.8, d = 1 sample (25 ns, their 50 bins of 0.5 ns), zero crossing interpolated
  ADC 50 % linear interpolation of the 50 % crossing on the rising edge
  template 20 % crossing of the aligned template (kcu0_sfil_fit_full.root)
    python3 dual_sided.py -> dual_sided.txt, fig_asym_dist.png, fig_asym_pos.png, fig_dt_dist.png, fig_dt_pos.png"""
import numpy as np, uproot, sys, ROOT
sys.path.insert(0, '/sphenix/user/tmengel/eic/bic_testbeam/common'); from plotstyle import setup, legend, headroom, finish, label
setup(); out = open('dual_sided.txt', 'w')
def say(*a): print(*a, flush=True); print(*a, file=out)
def core(v, w=2.5):
    v = v[np.isfinite(v)]; m = np.median(v); s = 1.4826 * np.median(np.abs(v - m))
    for _ in range(6): g = np.abs(v - m) < w * s; m, s = v[g].mean(), v[g].std()
    return m, s / 0.9
PX = {4: -80., 3: -30., 6: 30., 5: 80.}; PADS = (4, 3, 6, 5)
# 2017 note: positions are 10.1 mm per index step (run table); mean asymmetry read off their Fig. 7; mean time differences from Fig. 9
NOTE_X = np.arange(-8, 9, 2) * 10.1
NOTE_R = np.array([-0.0072, -0.0071, -0.0045, -0.0009, 0.0002, 0.0016, 0.0031, 0.0042, 0.0066])
NOTE_DT = {-80.8: -0.229 * 0.5, 80.8: 4.502 * 0.5}
fit = uproot.open('/sphenix/tg/tg01/jets/tmengel/eic/kcu0/fit_full/kcu0_sfil_fit_full.root')['fits'].arrays(['trig', 'pad', 'ok', 't0'], library='np')
trig = fit['trig']; pad = fit['pad']
f = uproot.open('/sphenix/user/tmengel/eic/bic_testbeam/data/kcu0_3062.root'); ts0 = f['header']['ts0'].array(library='np')[trig]
ev = uproot.open('/sphenix/user/tmengel/eic/bic_testbeam/data/events_3062.root')['events'].arrays(['trig_ts', 'astro_x', 'astro_purity'], library='np')
order = np.argsort(ev['trig_ts']); es = ev['trig_ts'][order]; pos = np.clip(np.searchsorted(es, ts0), 0, len(es) - 1); hit = es[pos] == ts0
ax = np.full(len(trig), np.nan); ax[hit] = ev['astro_x'][order[pos][hit]] - 88.75; pur = np.zeros(len(trig)); pur[hit] = ev['astro_purity'][order[pos][hit]]
astro = np.isfinite(ax) & (pur >= 0.3)
def end(tree, br, ifit):
    W = f[tree][br + '_adc'].array(library='np')[trig].astype(float); n = len(W); good = W[:, 0] >= 0
    raw_max = W.max(axis=1); y = W - W[:, 1:4].mean(axis=1)[:, None]
    ipk = y[:, 4:16].argmax(axis=1) + 4; A = y[np.arange(n), ipk]
    k = ipk.copy()                                            # first sample at or above half height
    for d in range(1, 6):
        kk = np.clip(ipk - d, 0, 19); k = np.where(y[np.arange(n), kk] >= 0.5 * A, kk, k)
    k = np.clip(k, 1, 18); y0 = y[np.arange(n), k - 1]; y1 = y[np.arange(n), k]
    t50 = 25. * (k - 1 + np.clip((0.5 * A - y0) / np.maximum(y1 - y0, 1e-3), 0, 1))
    c = y[:, :-1] - 0.8 * y[:, 1:]                            # CFD, d = 1 sample
    tcfd = np.full(n, np.nan)
    for j in range(2, 18):                                    # first negative -> positive crossing at or after the rise
        m = np.isnan(tcfd) & (j >= k - 1) & (c[:, j] < 0) & (c[:, j + 1] >= 0)
        tcfd[m] = 25. * (j + (-c[m, j]) / (c[m, j + 1] - c[m, j]))
    ttpl = np.where(fit['ok'][:, ifit], fit['t0'][:, ifit], np.nan)   # already in ns
    sel = good & (A > 100) & (ipk >= 5) & (ipk <= 14) & (raw_max < 1000)
    return dict(A=A, t={'CFD': tcfd, 'ADC 50 %': t50, 'template': ttpl}, sel=sel)
res = {}
for L, tree in ((2, 'sfil2'), (3, 'sfil3')):
    S = end(tree, 'Sc3', 8 * (L - 1) + 2); N = end(tree, 'Nc3', 8 * (L - 1) + 6)
    m = S['sel'] & N['sel'] & np.isin(pad, PADS)
    def asym(g): return (g * S['A'] - N['A']) / (g * S['A'] + N['A'])
    # gain: mean asymmetry 0 at x = 0 (P3, P6 average), solved by bisection on the South scale
    lo, hi = 0.2, 5.
    for _ in range(50):
        g = 0.5 * (lo + hi); r0 = 0.5 * (asym(g)[m & (pad == 3)].mean() + asym(g)[m & (pad == 6)].mean())
        lo, hi = (g, hi) if r0 < 0 else (lo, g)
    r = asym(g); say('SFIL%d cell 3: %d showers (one paddle, neither end clipped); South gain factor %.3f' % (L, m.sum(), g))
    rp = {p: (r[m & (pad == p)].mean(), r[m & (pad == p)].std() / np.sqrt((m & (pad == p)).sum()), r[m & (pad == p)].std(), (m & (pad == p)).sum()) for p in PADS}
    for p in PADS: say('   P%d x %+4.0f: n %7d  mean asymmetry %+.4f +- %.4f  rms %.3f' % (p, PX[p], rp[p][3], rp[p][0], rp[p][1], rp[p][2]))
    sl = np.polyfit([PX[p] for p in PADS], [rp[p][0] for p in PADS], 1)[0]
    say('   asymmetry slope %.2e /mm -> span over 160 mm %.3f (2017 note: %.3f); lambda_eff = 1/slope = %.0f mm' % (sl, 160 * sl, NOTE_R[-1] - NOTE_R[0], 1 / sl))
    ab = np.arange(-35, 96, 10.); ac = 0.5 * (ab[1:] + ab[:-1]); am = []
    for a0, a1 in zip(ab[:-1], ab[1:]):
        mm = m & astro & (ax >= a0) & (ax < a1); am.append((r[mm].mean(), r[mm].std() / np.sqrt(max(mm.sum(), 1)), mm.sum()) if mm.sum() > 100 else (np.nan, np.nan, mm.sum()))
    say('   AstroPix x bins: ' + ', '.join('%+.0f: %+.4f (%d)' % (x, a[0], a[2]) for x, a in zip(ac, am)))
    tp = {}
    for meth in ('CFD', 'ADC 50 %', 'template'):
        dt = N['t'][meth] - S['t'][meth]; mt = m & np.isfinite(dt) & (np.abs(dt) < 60)
        st = {p: (np.mean(dt[mt & (pad == p)]), core(dt[mt & (pad == p)])[1], np.std(dt[mt & (pad == p)]), (mt & (pad == p)).sum()) for p in PADS}
        x = np.array([PX[p] for p in PADS]); y_ = np.array([core(dt[mt & (pad == p)])[0] for p in PADS]); c1, c0 = np.polyfit(x, y_, 1)
        say('   %-9s n %7d  N - S per paddle (core mean / core sigma / rms, ns): %s;  slope %.4f ns/mm -> v_eff %.0f mm/ns; P5 - P4 %.2f ns (2017: %.2f ns over 162 mm)' % (meth, mt.sum(),
            '  '.join('P%d %+.2f/%.2f/%.2f' % (p, y_[i], st[p][1], st[p][2]) for i, p in enumerate(PADS)), c1, 2 / c1, y_[3] - y_[0], NOTE_DT[80.8] - NOTE_DT[-80.8]))
        tp[meth] = dict(dt=dt, m=mt, y=y_, c0=c0, c1=c1, sig=[st[p][1] for p in PADS])
    res[L] = dict(r=r, m=m, rp=rp, ac=ac, am=am, tp=tp, g=g)
COL = {4: ROOT.kRed + 1, 5: ROOT.kBlue + 1}; keep = []
def draw_hist(vals, lo, hi, nb, col, style):
    h = ROOT.TH1D('h%d' % len(keep), '', nb, lo, hi); h.FillN(len(vals), np.ascontiguousarray(vals, dtype=float), ROOT.nullptr); h.Scale(1. / max(h.Integral(), 1))
    h.SetLineColor(col); h.SetLineWidth(3); h.SetLineStyle(style); keep.append(h); return h
# Fig 1: asymmetry distributions at the two extreme positions (their Fig. 6)
c = ROOT.TCanvas('c1', '', 1600, 800); c.Divide(2, 1)
for i, L in enumerate((2, 3)):
    c.cd(i + 1); R = res[L]; hs = [draw_hist(R['r'][R['m'] & (pad == p)], -0.6, 0.6, 96, COL[p], 1) for p in (4, 5)]
    fr = ROOT.gPad.DrawFrame(-0.6, 0, 0.6, 1, 'SFIL%d cell 3;asymmetry (S_{S} - S_{N}) / (S_{S} + S_{N});fraction of showers' % L); headroom(fr, 0, max(h.GetMaximum() for h in hs), 0.35); keep.append(fr)
    lg = legend(0.20, 0.72, 0.60, 0.86)
    for h, p in zip(hs, (4, 5)): h.Draw('hist same'); lg.AddEntry(h, 'P%d (x = %+.0f mm), mean %+.3f' % (p, PX[p], R['rp'][p][0]), 'l')
    lg.Draw(); keep.append(lg)
finish(c); c.SaveAs('fig_asym_dist.png'); c.SaveAs('fig_asym_dist.pdf')
# Fig 2: mean asymmetry against position (their Fig. 7), paddles + AstroPix x
c = ROOT.TCanvas('c2', '', 1100, 850)
fr = ROOT.gPad.DrawFrame(-100, -0.12, 110, 0.12, 'mean light asymmetry against beam position;x along the fibres [mm];mean asymmetry'); headroom(fr, -0.12, 0.10, 0.30); keep.append(fr)
lg = legend(0.20, 0.72, 0.62, 0.86, size=0.036)
for L, mk, cl in ((2, 20, ROOT.kRed + 1), (3, 21, ROOT.kBlue + 1)):
    R = res[L]; x = np.array([PX[p] for p in PADS], float); y_ = np.array([R['rp'][p][0] for p in PADS]); e = np.array([R['rp'][p][1] for p in PADS])
    g = ROOT.TGraphErrors(4, x, y_, np.zeros(4), e); g.SetMarkerStyle(mk); g.SetMarkerSize(1.8); g.SetMarkerColor(cl); g.SetLineColor(cl); g.Draw('p same'); keep.append(g); lg.AddEntry(g, 'SFIL%d, paddles' % L, 'p')
    ok = np.array([np.isfinite(a[0]) for a in R['am']]); ga = ROOT.TGraphErrors(int(ok.sum()), R['ac'][ok].astype(float), np.array([a[0] for a in R['am']])[ok], np.zeros(int(ok.sum())), np.array([a[1] for a in R['am']])[ok])
    ga.SetMarkerStyle(mk + 4); ga.SetMarkerSize(1.4); ga.SetMarkerColor(cl); ga.SetLineColor(cl); ga.Draw('p same'); keep.append(ga); lg.AddEntry(ga, 'SFIL%d, AstroPix x' % L, 'p')
lg.Draw(); finish(c); c.SaveAs('fig_asym_pos.png'); c.SaveAs('fig_asym_pos.pdf')
# Fig 3: N - S time distributions at the extreme positions (their Fig. 9), one panel per method, SFIL2
c = ROOT.TCanvas('c3', '', 2400, 800); c.Divide(3, 1)
for i, meth in enumerate(('CFD', 'ADC 50 %', 'template')):
    c.cd(i + 1); T = res[2]['tp'][meth]; hs = [draw_hist(T['dt'][T['m'] & (pad == p)], -15, 15, 120, COL[p], 1) for p in (4, 5)]
    fr = ROOT.gPad.DrawFrame(-15, 0, 15, 1, 'SFIL2 cell 3, %s;t_{N} - t_{S} [ns];fraction of showers' % meth); headroom(fr, 0, max(h.GetMaximum() for h in hs), 0.35); keep.append(fr)
    lg = legend(0.20, 0.72, 0.62, 0.86, size=0.040)
    for h, p, k_ in zip(hs, (4, 5), (0, 3)): h.Draw('hist same'); lg.AddEntry(h, 'P%d: mean %+.2f, #sigma %.2f ns' % (p, T['y'][k_], T['sig'][k_]), 'l')
    lg.Draw(); keep.append(lg)
finish(c); c.SaveAs('fig_dt_dist.png'); c.SaveAs('fig_dt_dist.pdf')
# Fig 4: mean N - S time against position, relative to x = 0, all methods
c = ROOT.TCanvas('c4', '', 1100, 850)
fr = ROOT.gPad.DrawFrame(-100, -1.8, 110, 1.8, 'mean time difference against beam position;x along the fibres [mm];#LTt_{N} - t_{S}#GT - offset at x = 0 [ns]'); headroom(fr, -1.8, 1.4, 0.30); keep.append(fr)
lg = legend(0.20, 0.70, 0.75, 0.86, size=0.034)
x = np.array([PX[p] for p in PADS], float)
for meth, mk, cl in (('CFD', 20, ROOT.kRed + 1), ('ADC 50 %', 21, ROOT.kBlue + 1), ('template', 22, ROOT.kGreen + 2)):
    T = res[2]['tp'][meth]; g = ROOT.TGraph(4, x, (T['y'] - T['c0']).astype(float)); g.SetMarkerStyle(mk); g.SetMarkerSize(1.8); g.SetMarkerColor(cl); g.SetLineColor(cl); g.Draw('pl same'); keep.append(g)
    lg.AddEntry(g, 'SFIL2 %s: %.0f mm/ns' % (meth, 2 / T['c1']), 'pl')
lg.Draw(); finish(c); c.SaveAs('fig_dt_pos.png'); c.SaveAs('fig_dt_pos.pdf')
say('wrote fig_asym_dist, fig_asym_pos, fig_dt_dist, fig_dt_pos (.png/.pdf)')
