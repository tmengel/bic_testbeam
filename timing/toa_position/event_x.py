"""Per-event position along the fibres for every shower, from two data-only handles on the beam cell of each layer:
  x_amp  = -(lambda_eff / 2) (ln(A_N/A_S) - ln(g_N/g_S))      lambda_eff and g_N/g_S from astro_x_attenuation (SFIL2, SFIL3)
  x_time = (v_eff / 2) (dt_NS - dt_0)                          dt_NS from the ADC 50 % crossings, v_eff = 158 mm/ns, dt_0 per layer
Calibrated and validated against the four paddle positions (-80 / -30 / +30 / +80 mm) and, where present, the AstroPix x.
Output: event_x.npz (per fit-cache shower: x_amp per layer, x_time per layer, combined x) and event_x.txt/.png.
    python3 event_x.py [nmax]
"""
import numpy as np, uproot, sys, ROOT
sys.path.insert(0, '/sphenix/user/tmengel/eic/bic_testbeam/common'); from plotstyle import setup, legend, label
setup(); ROOT.gROOT.SetBatch(True); out = open('event_x.txt', 'w')
def say(*a): print(*a, flush=True); print(*a, file=out)
NMAX = int(sys.argv[1]) if len(sys.argv) > 1 else None
D = np.load('/sphenix/user/tmengel/eic/bic_testbeam/sim/analysis/data_fit_layers.npz'); trig = D['trig']; S3, N3, S, pad = D['S3'].astype(float), D['N3'].astype(float), D['S'].astype(float), D['pad']
if NMAX: keep = trig < NMAX; trig, S3, N3, S, pad = trig[keep], S3[keep], N3[keep], S[keep], pad[keep]
sat, satN = D['satS'][:len(trig)] if NMAX else D['satS'], D['satN'][:len(trig)] if NMAX else D['satN']
f = uproot.open('/sphenix/user/tmengel/eic/bic_testbeam/data/kcu0_3062.root'); ts0 = f['header']['ts0'].array(library='np', entry_stop=NMAX)
ev = uproot.open('/sphenix/user/tmengel/eic/bic_testbeam/data/events_3062.root')['events'].arrays(['trig_ts', 'astro_x', 'astro_purity'], library='np')
ets = ev['trig_ts']; order = np.argsort(ets); ets_s = ets[order]; t_ = ts0[trig]; pos = np.clip(np.searchsorted(ets_s, t_), 0, len(ets_s) - 1); hit = ets_s[pos] == t_; row = np.where(hit, order[pos], -1)
ax = np.full(len(trig), np.nan); ax[hit] = ev['astro_x'][row[hit]] - 88.75; pur = np.zeros(len(trig)); pur[hit] = ev['astro_purity'][row[hit]]
PX = {4: -80., 3: -30., 6: 30., 5: 80.}; xpad = np.array([PX.get(int(p), np.nan) for p in pad])
shower = S.sum(axis=1) > 600
# ---- ADC crossing times per layer and end (cell 3), for the fit-cache triggers
def crossing(tree, br):
    a = f[tree].arrays([br + '_adc'], entry_stop=NMAX, library='np')[br + '_adc'].astype(float)
    good = a[:, 0] >= 0; ped = a[:, 1:4].mean(axis=1); y = a - ped[:, None]; ipk = y[:, 4:13].argmax(axis=1) + 4; A = y[np.arange(len(y)), ipk]; half = 0.5 * A; k = ipk.copy()
    for d in range(1, 5):
        kk = np.clip(ipk - d, 0, 19); k = np.where(y[np.arange(len(y)), kk] >= half, kk, k)
    k = np.clip(k, 1, 19); y0 = y[np.arange(len(y)), k - 1]; y1 = y[np.arange(len(y)), k]; frac = np.clip((half - y0) / np.maximum(y1 - y0, 1e-3), 0, 1)
    t = np.where(good & (A > 100) & (ipk >= 6) & (ipk <= 10), 25. * (k - 1 + frac), np.nan); return t[trig], A[trig]
LAM = {1: 1329., 2: 1646.}; GNS = {1: 0.851, 2: 1.032}; VEFF = 158.
xa = np.full((len(trig), 4), np.nan); xt = np.full((len(trig), 4), np.nan)
say('calibration of each handle on the paddle positions (showers, beam cell unsaturated at both ends where amplitude is used):')
for L in range(4):
    tS, AS = crossing('sfil%d' % (L + 1), 'Sc3'); tN, AN = crossing('sfil%d' % (L + 1), 'Nc3'); dt = tN - tS
    ok_t = shower & np.isfinite(dt) & np.isfinite(xpad)
    # per-layer offset and slope from the paddle medians, then x_time with the fibre speed fixed
    med = {p: np.median(dt[ok_t & (pad == p)]) for p in (4, 3, 6, 5) if (ok_t & (pad == p)).sum() > 500}
    if len(med) >= 3:
        xs = np.array([PX[p] for p in med]); ys = np.array([med[p] for p in med]); pf = np.polyfit(xs, ys, 1); dt0 = pf[1]
        xt[:, L] = (dt - dt0) * VEFF / 2
        res = {p: np.std(np.clip(xt[ok_t & (pad == p), L] - PX[p], -400, 400)) for p in med}
        say('  SFIL%d time: fitted slope %.4f ns/mm (v_eff %.0f), offset %.2f ns; x_time - x_paddle rms per paddle %s mm; usable %d' % (L + 1, pf[0], 2 / pf[0], dt0, {p: int(v) for p, v in res.items()}, ok_t.sum()))
    if L in LAM:
        ok_a = shower & (S3[:, L] > 50) & (N3[:, L] > 50) & (sat[:, L] == 0) & (satN[:, L] == 0)
        xa[ok_a, L] = -(LAM[L] / 2) * (np.log(N3[ok_a, L] / S3[ok_a, L]) - np.log(GNS[L]))
        res = {p: np.std(np.clip(xa[ok_a & (pad == p), L] - PX[p], -400, 400)) for p in (4, 3, 6, 5)}; bias = {p: np.median(xa[ok_a & (pad == p), L]) - PX[p] for p in (4, 3, 6, 5)}
        say('  SFIL%d amplitude: x_amp - x_paddle median per paddle %s mm, rms %s mm; usable %d' % (L + 1, {p: int(v) for p, v in bias.items()}, {p: int(v) for p, v in res.items()}, ok_a.sum()))
# ---- combination: inverse-variance weights from the per-paddle rms of each handle (global numbers)
w = {}
for L in range(4):
    for arr, nm in ((xt, 't'), (xa, 'a')):
        m = shower & np.isfinite(arr[:, L]) & np.isfinite(xpad)
        if m.sum() > 1000: w[(nm, L)] = 1. / max(20., np.std(np.clip(arr[m, L] - xpad[m], -400, 400))) ** 2
num = np.zeros(len(trig)); den = np.zeros(len(trig))
for (nm, L), wt in w.items():
    arr = xt if nm == 't' else xa; m = np.isfinite(arr[:, L]); num[m] += wt * arr[m, L]; den[m] += wt
xc = np.where(den > 0, num / np.maximum(den, 1e-12), np.nan)
say('\ncombined x (weights 1/rms^2 of each handle): defined for %d of %d showers' % ((shower & np.isfinite(xc)).sum(), shower.sum()))
for p in (4, 3, 6, 5):
    m = shower & np.isfinite(xc) & (pad == p); v = xc[m] - PX[p]; say('  P%d (x %+.0f): n %7d  median residual %+.0f mm, rms %.0f mm, 16-84 %% %+.0f / %+.0f' % (p, PX[p], m.sum(), np.median(v), np.std(np.clip(v, -400, 400)), np.percentile(v, 16), np.percentile(v, 84)))
m = shower & np.isfinite(xc) & np.isfinite(ax) & (pur >= 0.3); v = xc[m] - ax[m]
say('  against the AstroPix x (%d showers): median %+.0f mm, rms %.0f mm, corr %.3f' % (m.sum(), np.median(v), np.std(np.clip(v, -400, 400)), np.corrcoef(xc[m], ax[m])[0, 1]))
np.savez('event_x.npz', trig=trig, x_amp=xa, x_time=xt, x=xc, xpad=xpad, astro_x=ax)
c = ROOT.TCanvas('c', '', 2000, 700); c.Divide(3, 1); keep = []
c.cd(1); fr = ROOT.gPad.DrawFrame(-300, 0, 300, 0.12, ';combined x per shower [mm];fraction'); keep.append(fr); lg = legend(0.6, 0.65, 0.88, 0.88)
for p, colr in zip((4, 3, 6, 5), (ROOT.kBlue + 1, ROOT.kGreen + 2, ROOT.kOrange + 7, ROOT.kRed + 1)):
    m = shower & np.isfinite(xc) & (pad == p); h = ROOT.TH1F('h%d' % p, '', 60, -300, 300); [h.Fill(v) for v in xc[m]]; h.Scale(1. / max(1, h.GetEntries())); h.SetLineColor(colr); h.SetLineWidth(2); h.Draw('hist same'); keep.append(h); lg.AddEntry(h, 'P%d (x %+.0f mm)' % (p, PX[p]), 'l')
fr.SetMaximum(max(k.GetMaximum() for k in keep[-4:]) * 1.3); lg.Draw(); keep.append(lg)
c.cd(2); m = shower & np.isfinite(xc) & np.isfinite(ax) & (pur >= 0.3); h2 = ROOT.TH2F('h2', ';AstroPix x [mm];combined x [mm]', 40, -60, 100, 60, -300, 300); [h2.Fill(a, b) for a, b in zip(ax[m], xc[m])]; h2.Draw('col'); keep.append(h2)
c.cd(3); fr3 = ROOT.gPad.DrawFrame(-100, -100, 100, 100, ';paddle x [mm];median combined x [mm]'); keep.append(fr3)
g = ROOT.TGraph(4, np.array([-80., -30., 30., 80.]), np.array([np.median(xc[shower & np.isfinite(xc) & (pad == p)]) for p in (4, 3, 6, 5)])); g.SetMarkerStyle(20); g.Draw('p same'); keep.append(g); ln = ROOT.TLine(-100, -100, 100, 100); ln.SetLineStyle(2); ln.Draw(); keep.append(ln)
c.SaveAs('event_x.png'); c.SaveAs('event_x.pdf'); say('wrote event_x.png')
