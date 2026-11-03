"""SFIL2 cell 3 TOA without the ADC in the timing (replaces toa_sfil2/ for the TOA conclusions).
The ADC is used only to find the pulse (amplitude 150-950, peak in samples 5-14, no clipping) and the 25 ns window that holds
it (the 50 % crossing sample k). The TOA time of an end is
    t = 25 s + T(code)
with s the sample the word sits on (the word on k, else the one on k+1) and T the code-density calibration: the unwrapped
code distribution of pulse words, built on even triggers, mapped onto 0-25 ns assuming the pulse phase is uniform (the beam
is not synchronous with the 40 MHz clock). No walk, baseline or ADC-window correction. Every result is on odd triggers.
The ADC 50 % crossing is computed alongside as a separate, independent estimator for comparison only.
    python3 toa_clean.py -> toa_clean.txt, toa_clean.npz"""
import numpy as np, uproot
out = open('toa_clean.txt', 'w')
def say(*a): print(*a, flush=True); print(*a, file=out)
NMAX = 1200000
f = uproot.open('/sphenix/user/tmengel/eic/bic_testbeam/data/kcu0_3062.root')
ts0 = f['header']['ts0'].array(library='np', entry_stop=NMAX)
ev = uproot.open('/sphenix/user/tmengel/eic/bic_testbeam/data/events_3062.root')['events'].arrays(['trig_ts', 'pad', 'astro_x', 'astro_purity'], library='np')
order = np.argsort(ev['trig_ts']); es = ev['trig_ts'][order]; pos = np.clip(np.searchsorted(es, ts0), 0, len(es) - 1); hit = es[pos] == ts0
row = order[pos]; pad = np.where(hit, ev['pad'][row], 0); ax = np.where(hit, ev['astro_x'][row] - 88.75, np.nan); pur = np.where(hit, ev['astro_purity'][row], 0.)
PB = {4: 2, 3: 1, 6: 8, 5: 4}; PX = {4: -80., 3: -30., 6: 30., 5: 80.}
NT = len(ts0); even = (np.arange(NT) % 2) == 0
def end(br):
    a = f['sfil2'].arrays([br + '_adc', br + '_toa'], entry_stop=NMAX, library='np'); W = a[br + '_adc'].astype(float); T = a[br + '_toa'].astype(int)
    n = len(W); good = W[:, 0] >= 0; y = W - W[:, 1:4].mean(axis=1)[:, None]
    ipk = y[:, 4:16].argmax(axis=1) + 4; A = y[np.arange(n), ipk]; k = ipk.copy()
    for d in range(1, 6):
        kk = np.clip(ipk - d, 0, 19); k = np.where(y[np.arange(n), kk] >= 0.5 * A, kk, k)
    k = np.clip(k, 1, 18); y0 = y[np.arange(n), k - 1]; y1 = y[np.arange(n), k]
    tc = 25. * (k - 1 + np.clip((0.5 * A - y0) / np.maximum(y1 - y0, 1e-3), 0, 1))
    sel = good & (A > 150) & (A < 950) & (ipk >= 5) & (ipk <= 14) & (W.max(axis=1) < 1000)
    w0 = T[np.arange(n), k]; w1 = T[np.arange(n), k + 1]
    code = np.where(w0 > 0, w0, np.where(w1 > 0, w1, -1)); s = np.where(w0 > 0, k, k + 1); has = sel & (code > 0)
    nop = good & (A < 15); bg = T[nop][:, 1:19].ravel(); bg = bg[bg > 0]
    ped = W[:, 1:4].mean(axis=1); base = ped - np.median(ped[good])
    npre = ((T > 0) & (np.arange(20)[None, :] >= 1) & (np.arange(20)[None, :] < k[:, None])).sum(axis=1)   # words before the pulse
    return dict(npre=npre, base=base, A=A, tc=tc, sel=sel, has=has, code=code, s=s, k=k, both=sel & (w0 > 0) & (w1 > 0), bg=bg, bgrate=(T[nop][:, 1:19] > 0).mean(),
                pre=(T[np.arange(n), k - 1] > 0)[sel].mean())
E = {'S': end('Sc3'), 'N': end('Nc3')}
for e, P in E.items():
    h = np.bincount(P['code'][P['has'] & even], minlength=1024)[:1024].astype(float); w64 = np.add.reduceat(h, np.arange(0, 1024, 64))
    c0 = 0
    if w64.min() < 0.02 * w64.mean():                                 # an empty code band = where the counter wraps
        lo = np.argmin(w64) * 64; z = np.where(h[lo:lo + 64] == 0)[0]; c0 = lo + (z[0] if len(z) else 0)
    hu = np.roll(h, -c0); cdf = np.cumsum(hu) / hu.sum(); tab = 25. * (cdf - 0.5 * hu / hu.sum())
    P['c0'] = c0; P['tab'] = tab; P['hu'] = hu
    t = np.where(P['has'], 25. * P['s'] + tab[(P['code'] - c0) % 1024], np.nan); P['t'] = t
    say('%s end: %d clean pulses; TOA word on k or k+1 for %.2f (both %.3f); word on k-1 (before the pulse) %.2f; no-pulse word rate %.2f per sample; unwrap c0 %d' %
        (e, P['sel'].sum(), P['has'][P['sel']].mean(), P['both'][P['sel']].mean(), P['pre'], P['bgrate'], c0))
both = E['S']['has'] & E['N']['has'] & ~even & np.isin(pad, list(PB.values()))
both_adc = E['S']['sel'] & E['N']['sel'] & ~even & np.isin(pad, list(PB.values()))
say('odd-trigger showers with a clean pulse at both ends and one paddle: %d; with a TOA word at both ends: %d (%.2f)' % (both_adc.sum(), both.sum(), both.sum() / both_adc.sum()))
dT0 = E['N']['t'] - E['S']['t']; dA = E['N']['tc'] - E['S']['tc']
# amplitude and baseline correction of the TOA N - S, using ADC amplitudes and pedestals only (no ADC time): regression of the
# within-paddle residual on ln A_S, ln A_N, baseline_S, baseline_N, fitted on even triggers inside +-10 ns, applied to odd triggers
padok = np.isin(pad, list(PB.values())); bothall = E['S']['has'] & E['N']['has'] & padok
pm = {p: np.median(dT0[bothall & even & (pad == PB[p])]) for p in PX}
res0 = dT0 - np.select([pad == PB[p] for p in PX], [pm[p] for p in PX], np.nan)
X = np.stack([np.ones(NT), np.log(np.maximum(E['S']['A'], 1)), np.log(np.maximum(E['N']['A'], 1)), E['S']['base'], E['N']['base']], axis=1)
fm = bothall & even & (np.abs(res0) < 10); beta = np.linalg.lstsq(X[fm], res0[fm], rcond=None)[0]
dT = dT0 - X[:, 1:] @ beta[1:]
say('amplitude/baseline correction (even triggers): %.2f ns per ln A_S, %.2f per ln A_N, %.3f ns per ADC of baseline S, %.3f N' % tuple(beta[1:]))
say('raw TOA times: corr(t_N, t_S) %.3f over %d showers (spread of t_S %.1f ns rms, of t_N - t_S %.1f ns rms)' % (np.corrcoef(E['N']['t'][both], E['S']['t'][both])[0, 1], both.sum(), np.std(E['S']['t'][both]), np.std(dT[both])))
# N - S per paddle: Gaussian core + flat background, fitted by EM on +-25 ns
def gfit(v, lo=-25., hi=25.):
    v = v[(v > lo) & (v < hi)]; m, s, fsig = np.median(v), 2.5, 0.8
    for _ in range(200):
        g = fsig * np.exp(-0.5 * ((v - m) / s) ** 2) / (s * np.sqrt(2 * np.pi)); b = (1 - fsig) / (hi - lo); r = g / (g + b)
        fsig = r.mean(); m = (r * v).sum() / r.sum(); s = np.sqrt((r * (v - m) ** 2).sum() / r.sum())
    return m, s, fsig, len(v)
R = {}
for lab, d, m in (('TOA raw', dT0, both), ('TOA', dT, both), ('ADC', dA, both)):
    fits = {p: gfit(d[m & (pad == PB[p])]) for p in PX}
    x = np.array([PX[p] for p in (4, 3, 6, 5)]); mu = np.array([fits[p][0] for p in (4, 3, 6, 5)]); c1, c0 = np.polyfit(x, mu, 1)
    sig = np.array([fits[p][1] for p in (4, 3, 6, 5)]); sep = (mu[3] - mu[0]) / sig.mean()
    say('%s N - S, same %d showers: per paddle P4/P3/P6/P5 mean %s ns, core sigma %s ns, core fraction %s; slope %.4f ns/mm -> v_eff %.0f mm/ns; sigma_x %s mm; P5 vs P4 separated by %.2f sigma' %
        (lab, m.sum(), ' / '.join('%+.2f' % v for v in mu), ' / '.join('%.2f' % v for v in sig), ' / '.join('%.2f' % fits[p][2] for p in (4, 3, 6, 5)), c1, 2 / c1, ' / '.join('%.0f' % (v / c1) for v in sig), sep))
    R[lab] = dict(mu=mu, sig=sig, frac=np.array([fits[p][2] for p in (4, 3, 6, 5)]), c0=c0, c1=c1)
# within-paddle correlation of the two independent estimators (shared light timing)
rT = np.full(NT, np.nan); rA = np.full(NT, np.nan)
for i, p in enumerate((4, 3, 6, 5)):
    mm = both & (pad == PB[p]); rT[mm] = dT[mm] - R['TOA']['mu'][i]; rA[mm] = dA[mm] - R['ADC']['mu'][i]
core = both & (np.abs(rT) < 2 * R['TOA']['sig'].mean()) & (np.abs(rA) < 2 * R['ADC']['sig'].mean())
say('within a paddle, TOA and ADC N - S residuals (both inside 2 sigma, %d showers): correlation %.2f' % (core.sum(), np.corrcoef(rT[core], rA[core])[0, 1]))
# per-end QA against the ADC crossing (diagnostic only, never used in the calibration)
for e, P in E.items():
    r = (P['t'] - P['tc'])[P['has'] & ~even]; r = r - np.median(r)
    say('%s end, diagnostic only: TOA - ADC crossing: within 5 ns %.2f, beyond 12.5 ns %.2f; core sigma %.2f ns' % (e, (np.abs(r) < 5).mean(), (np.abs(r) > 12.5).mean(), gfit(r)[1]))
# AstroPix
ma = both & np.isfinite(ax) & (pur >= 0.3)
for lab, d in (('TOA', dT), ('ADC', dA)):
    xx = (d - R[lab]['c0']) / R[lab]['c1']
    say('%s vs AstroPix (%d showers): mean x per paddle P3/P6/P5 %s mm (AstroPix %s); per shower rms of x - x_AstroPix inside +-400 mm %.0f mm' % (lab, ma.sum(),
        ' / '.join('%+.0f' % np.median(xx[ma & (pad == PB[p])]) for p in (3, 6, 5)), ' / '.join('%+.0f' % np.median(ax[ma & (pad == PB[p])]) for p in (3, 6, 5)),
        np.std(np.clip((xx - ax)[ma], -400, 400))))
np.savez('toa_clean.npz', dT0=dT0, beta=beta, pad=pad, even=even, ax=ax, pur=pur, both=both, dT=dT, dA=dA, **{k + e: E[e][k] for e in E for k in ('t', 'tc', 'A', 'has', 'code', 's', 'k', 'tab', 'hu', 'bg', 'npre', 'both', 'base')},
         c0S=E['S']['c0'], c0N=E['N']['c0'], **{'%s_%s' % (lab.replace(' ', '_'), k): v for lab in R for k, v in R[lab].items()})
say('wrote toa_clean.npz')
