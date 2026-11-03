"""Walk-corrected North - South time from the template-fit t0 (kcu0_sfil_fit_full.root: per-channel 20 % crossing of the
aligned template, clipped samples masked), combined with the amplitude-asymmetry position.
Per layer: dt = t0_N - t0_S.  Joint linear regression on paddle-tagged showers
    dt = c0 + c1 x_paddle + a ln(A_N/A_S)
gives the fibre slope c1 (= 2 / v_eff) and the residual walk coefficient a; x_time = (dt - a ln(A_N/A_S) - c0) / c1.
x_amp = -(lambda_eff/2)(ln(A_N/A_S) - ln g_N/g_S) for SFIL2, SFIL3 (measured 1329 / 1646 mm, 0.851 / 1.032).
    python3 event_x_walk.py -> event_x_walk.txt/.png, event_x_walk.npz
"""
import numpy as np, uproot, sys, ROOT
sys.path.insert(0, '/sphenix/user/tmengel/eic/bic_testbeam/common'); from plotstyle import setup, legend, label
setup(); ROOT.gROOT.SetBatch(True); out = open('event_x_walk.txt', 'w')
def say(*a): print(*a, flush=True); print(*a, file=out)
def medres(arr, sel, pad, PX): return {p: (int(np.median(arr[sel & (pad == p)] - PX[p])) if (sel & (pad == p)).sum() > 20 else None) for p in (4, 3, 6, 5)}
t = uproot.open('/sphenix/tg/tg01/jets/tmengel/eic/kcu0/fit_full/kcu0_sfil_fit_full.root')['fits']
a = t.arrays(['trig', 'pad', 'ok', 'afit', 'amp', 't0'], library='np'); trig, pad, ok, afit, amp, t0 = a['trig'], a['pad'], a['ok'], a['afit'], a['amp'], a['t0']
D = np.load('/sphenix/user/tmengel/eic/bic_testbeam/sim/analysis/data_fit_layers.npz'); assert np.array_equal(D['trig'], trig); S = D['S'].astype(float); shower = S.sum(axis=1) > 600
PX = {4: -80., 3: -30., 6: 30., 5: 80.}; xpad = np.array([PX.get(int(p), np.nan) for p in pad]); haspad = np.isfinite(xpad)
f = uproot.open('/sphenix/user/tmengel/eic/bic_testbeam/data/kcu0_3062.root'); ts0 = f['header']['ts0'].array(library='np')
ev = uproot.open('/sphenix/user/tmengel/eic/bic_testbeam/data/events_3062.root')['events'].arrays(['trig_ts', 'astro_x', 'astro_purity'], library='np')
ets = ev['trig_ts']; order = np.argsort(ets); ets_s = ets[order]; tt = ts0[trig]; pos = np.clip(np.searchsorted(ets_s, tt), 0, len(ets_s) - 1); hit = ets_s[pos] == tt; row = np.where(hit, order[pos], -1)
ax = np.full(len(trig), np.nan); ax[hit] = ev['astro_x'][row[hit]] - 88.75; pur = np.zeros(len(trig)); pur[hit] = ev['astro_purity'][row[hit]]
LAM = {1: 1329., 2: 1646.}; GNS = {1: 0.851, 2: 1.032}
xt = np.full((len(trig), 4), np.nan); xa = np.full((len(trig), 4), np.nan); info = {}
say('showers %d (with paddle %d)' % (shower.sum(), (shower & haspad).sum()))
for L in range(4):
    iS, iN = 8 * L + 2, 8 * L + 6
    m = shower & haspad & ok[:, iS] & ok[:, iN] & (afit[:, iS] > 100) & (afit[:, iN] > 100)
    dt = t0[:, iN] - t0[:, iS]; la = np.log(np.maximum(afit[:, iN], 1) / np.maximum(afit[:, iS], 1))
    med = np.median(dt[m]); g = m & (np.abs(dt - med) < 15)
    X = np.stack([np.ones(g.sum()), xpad[g], la[g]], axis=1); coef, *_ = np.linalg.lstsq(X, dt[g], rcond=None); c0, c1, aw = coef
    X2 = np.stack([np.ones(g.sum()), xpad[g]], axis=1); coef2, *_ = np.linalg.lstsq(X2, dt[g], rcond=None)
    # the explicit walk term is degenerate with x (ln(A_N/A_S) is itself the attenuation position), so the template t0 is used
    # without it: slope and offset from the two-parameter fit; the three-parameter fit is reported for the record only
    c0, c1 = coef2; aw = coef[2]; xt_L = (dt - c0) / c1; xt[g, L] = xt_L[g]
    clipped = (amp[:, iS] > 900) | (amp[:, iN] > 900)
    r_all = np.std(np.clip(xt_L[g] - xpad[g], -400, 400)); r_unc = np.std(np.clip(xt_L[g & ~clipped] - xpad[g & ~clipped], -400, 400)); r_cl = np.std(np.clip(xt_L[g & clipped] - xpad[g & clipped], -400, 400)) if (g & clipped).sum() > 100 else np.nan
    info[L] = (c0, c1, aw, r_all)
    say('  SFIL%d: n %7d; without walk term slope %.4f ns/mm (v_eff %.0f); three-parameter fit slope %.4f ns/mm (v_eff %.0f) with a = %.2f ns per unit ln(A_N/A_S) (degenerate with x, not applied); offset %.2f ns; x_time rms vs paddle: all %.0f, unclipped %.0f, clipped %.0f mm; per-paddle median residual %s' % (L + 1, g.sum(), coef2[1], 2 / coef2[1], coef[1], 2 / coef[1], aw, c0, r_all, r_unc, r_cl, medres(xt_L, g, pad, PX)))
    if L in LAM:
        ma = shower & ok[:, iS] & ok[:, iN] & (afit[:, iS] > 50) & (afit[:, iN] > 50) & (amp[:, iS] < 900) & (amp[:, iN] < 900)
        xa[ma, L] = -(LAM[L] / 2) * (la[ma] - np.log(GNS[L]))
        mm = ma & haspad; say('  SFIL%d amplitude handle: n %d, rms vs paddle %.0f mm, medians %s' % (L + 1, mm.sum(), np.std(np.clip(xa[mm, L] - xpad[mm], -400, 400)), medres(xa[:, L], mm, pad, PX)))
# combination with inverse-variance weights from the per-paddle rms
w = {}
for L in range(4):
    for arr, nm in ((xt, 't'), (xa, 'a')):
        m = shower & haspad & np.isfinite(arr[:, L])
        if m.sum() > 1000: w[(nm, L)] = 1. / max(30., np.std(np.clip(arr[m, L] - xpad[m], -400, 400))) ** 2
num = np.zeros(len(trig)); den = np.zeros(len(trig))
for (nm, L), wt in w.items():
    arr = xt if nm == 't' else xa; m = np.isfinite(arr[:, L]); num[m] += wt * arr[m, L]; den[m] += wt
xc = np.where(den > 0, num / np.maximum(den, 1e-12), np.nan)
say('\nweights (1/rms^2, rms in mm): ' + ', '.join('%s SFIL%d %.0f' % ('time' if nm == 't' else 'amp', L + 1, 1 / np.sqrt(wt)) for (nm, L), wt in sorted(w.items())))
say('combined x: defined for %d of %d showers' % ((shower & np.isfinite(xc)).sum(), shower.sum()))
for p in (4, 3, 6, 5):
    m = shower & np.isfinite(xc) & (pad == p); v = xc[m] - PX[p]; say('  P%d (x %+.0f): n %7d  median residual %+.0f mm, rms %.0f mm, 16-84 %% %+.0f / %+.0f' % (p, PX[p], m.sum(), np.median(v), np.std(np.clip(v, -400, 400)), np.percentile(v, 16), np.percentile(v, 84)))
# per-paddle resolution with the paddle width removed: paddle 23.8 mm wide -> sigma 6.9 mm, negligible
m = shower & np.isfinite(xc) & np.isfinite(ax) & (pur >= 0.3); v = xc[m] - ax[m]
say('  against the AstroPix x (%d showers, purity >= 0.3): median %+.0f mm, rms %.0f mm, corr %.3f' % (m.sum(), np.median(v), np.std(np.clip(v, -400, 400)), np.corrcoef(xc[m], ax[m])[0, 1]))
# time-only and amplitude-only combos for the record
for nm in ('t', 'a'):
    num = np.zeros(len(trig)); den = np.zeros(len(trig))
    for (n_, L), wt in w.items():
        if n_ != nm: continue
        arr = xt if nm == 't' else xa; mm = np.isfinite(arr[:, L]); num[mm] += wt * arr[mm, L]; den[mm] += wt
    xo = np.where(den > 0, num / np.maximum(den, 1e-12), np.nan); mm = shower & haspad & np.isfinite(xo)
    say('  %s-only combination: n %d, rms vs paddle %.0f mm, medians %s' % ('time' if nm == 't' else 'amplitude', mm.sum(), np.std(np.clip(xo[mm] - xpad[mm], -400, 400)), medres(xo, mm, pad, PX)))
np.savez('event_x_walk.npz', trig=trig, x_time=xt, x_amp=xa, x=xc, xpad=xpad, astro_x=ax, coef={L: info[L] for L in info}.__repr__())
c = ROOT.TCanvas('c', '', 2000, 700); c.Divide(3, 1); keep = []
c.cd(1); fr = ROOT.gPad.DrawFrame(-300, 0, 300, 0.12, ';combined x per shower [mm];fraction'); keep.append(fr); lg = legend(0.6, 0.65, 0.88, 0.88)
for p, colr in zip((4, 3, 6, 5), (ROOT.kBlue + 1, ROOT.kGreen + 2, ROOT.kOrange + 7, ROOT.kRed + 1)):
    m = shower & np.isfinite(xc) & (pad == p); h = ROOT.TH1F('h%d' % p, '', 60, -300, 300); [h.Fill(v) for v in xc[m]]; h.Scale(1. / max(1, h.GetEntries())); h.SetLineColor(colr); h.SetLineWidth(2); h.Draw('hist same'); keep.append(h); lg.AddEntry(h, 'P%d (x %+.0f mm)' % (p, PX[p]), 'l')
fr.SetMaximum(max(k.GetMaximum() for k in keep[-4:]) * 1.3); lg.Draw(); keep.append(lg)
c.cd(2); m = shower & np.isfinite(xc) & np.isfinite(ax) & (pur >= 0.3); h2 = ROOT.TH2F('h2', ';AstroPix x [mm];combined x [mm]', 40, -60, 100, 60, -300, 300); [h2.Fill(a_, b_) for a_, b_ in zip(ax[m], xc[m])]; h2.Draw('col'); keep.append(h2)
c.cd(3); fr3 = ROOT.gPad.DrawFrame(-100, -100, 100, 100, ';paddle x [mm];median x per shower [mm]'); keep.append(fr3)
for arr, mk, colr, nm in ((xc, 20, ROOT.kBlack, 'combined'), (xt[:, 1], 24, ROOT.kBlue + 1, 'SFIL2 time'), (xa[:, 1], 25, ROOT.kRed + 1, 'SFIL2 amplitude')):
    g = ROOT.TGraph(4, np.array([-80., -30., 30., 80.]), np.array([np.median(arr[shower & np.isfinite(arr) & (pad == p)]) for p in (4, 3, 6, 5)])); g.SetMarkerStyle(mk); g.SetMarkerColor(colr); g.Draw('p same'); keep.append(g)
ln = ROOT.TLine(-100, -100, 100, 100); ln.SetLineStyle(2); ln.Draw(); keep.append(ln); lg3 = legend(0.15, 0.68, 0.55, 0.88); [lg3.AddEntry(keep[-4 + i], nm, 'p') for i, nm in enumerate(('combined', 'SFIL2 time', 'SFIL2 amplitude'))]; lg3.Draw(); keep.append(lg3)
c.SaveAs('event_x_walk.png'); c.SaveAs('event_x_walk.pdf'); say('wrote event_x_walk.png')
