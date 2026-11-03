"""Second pass on the SFIL2 cell 3 KCU0 TOA: pedestal-word background subtracted from the code density, word chosen among
samples k, k+1 by consistency with the ADC crossing, walk in ln A, resolution per amplitude bin; TOA alone and TOA+ADC combined.
    python3 toa_cdc2.py [nmax] [cut_ns] -> toa_cdc2.txt/.png/.npz"""
import numpy as np, uproot, sys, ROOT
sys.path.insert(0, '/sphenix/user/tmengel/eic/bic_testbeam/common'); from plotstyle import setup, legend
setup(); ROOT.gROOT.SetBatch(True); out = open('toa_cdc2.txt', 'w')
def say(*a): print(*a, flush=True); print(*a, file=out)
NMAX = int(sys.argv[1]) if len(sys.argv) > 1 else 1200000
CUT = float(sys.argv[2]) if len(sys.argv) > 2 else 6.
f = uproot.open('/sphenix/user/tmengel/eic/bic_testbeam/data/kcu0_3062.root')
ts0 = f['header']['ts0'].array(library='np', entry_stop=NMAX)
ev = uproot.open('/sphenix/user/tmengel/eic/bic_testbeam/data/events_3062.root')['events'].arrays(['trig_ts', 'pad'], library='np')
ets = ev['trig_ts']; order = np.argsort(ets); ets_s = ets[order]; pos = np.clip(np.searchsorted(ets_s, ts0), 0, len(ets_s) - 1); hit = ets_s[pos] == ts0
pad = np.zeros(len(ts0), np.int16); pad[hit] = ev['pad'][order[pos][hit]]
PB = {4: 2, 3: 1, 6: 8, 5: 4}; PX = {4: -80., 3: -30., 6: 30., 5: 80.}
def msig(v): v = v[np.isfinite(v)]; return 1.4826 * np.median(np.abs(v - np.median(v)))
def core(v, w=2.5):
    v = v[np.isfinite(v)]
    if len(v) < 30: return np.nan, np.nan
    m, s = np.median(v), msig(v)
    for _ in range(6): g = np.abs(v - m) < w * s; m, s = v[g].mean(), v[g].std()
    return m, s / 0.9   # 2.5-sigma truncation of a Gaussian -> std x 0.94; keep conservative /0.9
def load(br):
    a = f['sfil2'].arrays([br + '_adc', br + '_toa'], entry_stop=NMAX, library='np'); W = a[br + '_adc'].astype(float); T = a[br + '_toa'].astype(int)
    n = len(W); good = W[:, 0] >= 0; y = W - W[:, 1:4].mean(axis=1)[:, None]
    ipk = y[:, 4:16].argmax(axis=1) + 4; A = y[np.arange(n), ipk]; half = 0.5 * A; k = ipk.copy()
    for d in range(1, 6):
        kk = np.clip(ipk - d, 0, 19); k = np.where(y[np.arange(n), kk] >= half, kk, k)
    k = np.clip(k, 1, 18); y0 = y[np.arange(n), k - 1]; y1 = y[np.arange(n), k]; tc = 25. * (k - 1 + np.clip((half - y0) / np.maximum(y1 - y0, 1e-3), 0, 1))
    sel = good & (A > 150) & (A < 950) & (ipk >= 5) & (ipk <= 14)
    w0 = T[np.arange(n), k]; w1 = T[np.arange(n), k + 1]
    noise = T[good][:, 1:4].ravel(); noise = noise[noise > 0]; rate = (T[good][:, 1:4] > 0).mean()
    return dict(A=A, tc=tc, k=k, w0=w0, w1=w1, sel=sel, noise=noise, rate=rate)
def cal_table(codes, noise, nscale):
    occ = np.bincount(codes, minlength=1024).astype(float); bg = np.bincount(noise, minlength=1024).astype(float); bg *= nscale / max(bg.sum(), 1)
    sig = np.clip(occ - bg, 0, None); cdf = np.cumsum(sig) / sig.sum(); return 25. * (cdf - 0.5 * sig / sig.sum())
def mode_table(codes, phase, c0):
    """ADC-referenced: per 16-code bin the mode of the crossing phase (robust to flat noise words), interpolated."""
    cb = ((codes - c0) % 1024) // 16; xs, ys = [], []
    for b in range(64):
        m = cb == b
        if m.sum() < 200: continue
        h, e = np.histogram(phase[m], bins=np.arange(-40, 40.01, 0.5)); h = np.convolve(h, np.ones(5) / 5, 'same'); xs.append(b * 16 + 8); ys.append(e[h.argmax()] + 0.25)
    return np.interp(np.arange(1024), xs, ys) if len(xs) > 3 else None
def calibrate(P, name):
    sel = P['sel']; n = len(P['A'])
    cand = [(P['w0'], P['k']), (P['w1'], P['k'] + 1)]
    codes_raw = np.concatenate([w[sel & (w > 0)] for w, s in cand]); sams = np.concatenate([s[sel & (w > 0)] for w, s in cand]); tcs = np.concatenate([P['tc'][sel & (w > 0)] for w, s in cand])
    nb0 = P['rate'] * 2 * sel.sum(); best = None
    for c0 in range(0, 1024, 16):
        for fs in (0., 0.25, 0.5, 0.75, 1.0):
            tab = cal_table((codes_raw - c0) % 1024, (P['noise'] - c0) % 1024, fs * nb0)
            r = 25. * sams + tab[(codes_raw - c0) % 1024] - tcs; m0 = np.median(r); h, e = np.histogram(r, bins=np.arange(-40, 40.01, 0.5)); m0 = e[np.convolve(h, np.ones(5), 'same').argmax()] + 0.25
            sc = np.mean(np.abs(r - m0) < 3.)
            if best is None or sc > best[0]: best = ('cdf', c0, tab, m0, fs) if False else (sc, c0, tab, m0, 'code density, noise scale %.2f' % fs)
    # ADC-referenced mode table at the best unwrap
    c0 = best[1]; ph = tcs - 25. * sams; mt = mode_table(codes_raw, ph, c0)
    if mt is not None:
        r = 25. * sams + mt[(codes_raw - c0) % 1024] - tcs; sc2 = np.mean(np.abs(r) < 3.)
        if sc2 > best[0]: best = (sc2, c0, mt, 0., 'ADC-referenced mode per 16 codes')
    sc, c0, tab, m0, how = best
    t = np.full(n, np.nan); r_best = np.full(n, np.inf)
    for w, s in cand:
        tt = np.where(sel & (w > 0), 25. * s + tab[(w - c0) % 1024] - m0, np.nan); rr = np.abs(tt - P['tc'])
        better = np.isfinite(tt) & (rr < r_best); t[better] = tt[better]; r_best[better] = rr[better]
    ok = sel & (r_best < CUT)
    la = np.log(np.maximum(P['A'], 1)); bins = np.linspace(np.log(150), np.log(950), 11); xb, yb = [], []
    for lo_, hi_ in zip(bins[:-1], bins[1:]):
        m = ok & (la >= lo_) & (la < hi_)
        if m.sum() > 300: xb.append(np.median(la[m])); yb.append(core((t - P['tc'])[m])[0])
    tw = t - np.interp(la, xb, yb)
    m_, s_ = core((tw - P['tc'])[ok])
    say('%s: pulses %d; pedestal word rate %.3f/sample; unwrap c0 %d; table: %s; frac of words within 3 ns %.2f; kept (|r| < %.0f ns) %d = %.2f of pulses; TOA - crossing core sigma %.2f ns; walk span %.1f ns' %
        (name, sel.sum(), P['rate'], c0, how, sc, CUT, ok.sum(), ok.sum() / sel.sum(), s_, max(yb) - min(yb)))
    # diagnostics: TOA - crossing per code range and per amplitude (is the error local in code, or everywhere?)
    cw = np.full(n, -1)
    for w, s_ in cand:
        tt = np.where(sel & (w > 0), 25. * s_ + tab[(w - c0) % 1024] - m0, np.nan); hitm = np.isfinite(tt) & (np.abs(tt - t) < 1e-9); cw[hitm] = ((w - c0) % 1024)[hitm]
    rr = tw - P['tc']; parts = []
    for b in range(0, 1024, 128):
        m = ok & (cw >= b) & (cw < b + 128); parts.append('%d-%d: %.2f (%d)' % (b, b + 127, core(rr[m])[1], m.sum()))
    say('   per unwrapped-code range, core sigma TOA - crossing [ns]: ' + ';  '.join(parts))
    parts = []
    for lo_, hi_ in ((150, 250), (250, 400), (400, 600), (600, 950)):
        m = ok & (P['A'] >= lo_) & (P['A'] < hi_); parts.append('A %d-%d: %.2f' % (lo_, hi_, core(rr[m])[1]))
    say('   per amplitude: ' + ';  '.join(parts))
    P.update(t=tw, ok=ok, c0=c0, tab=tab - (tab[0] if how.startswith('ADC') else 0), xb=xb, yb=yb, how=how)
    return P
S = calibrate(load('Sc3'), 'SFIL2 S c3'); N = calibrate(load('Nc3'), 'SFIL2 N c3')
both = S['ok'] & N['ok']; dT = N['t'] - S['t']; dA = N['tc'] - S['tc']
say('\nboth ends with a good TOA: %d of %d showers with both pulses (%.2f)' % (both.sum(), (S['sel'] & N['sel']).sum(), both.sum() / (S['sel'] & N['sel']).sum()))
# remove the paddle means (position) with each estimator's own means
def resid(d, m):
    r = np.full(len(d), np.nan); means = []
    for p in (4, 3, 6, 5):
        mm = m & (pad == PB[p]); mu = core(d[mm])[0]; r[mm] = d[mm] - mu; means.append(mu)
    return r, means
rT, mT = resid(dT, both); rA, mA = resid(dA, both)
for lab, mm in (('TOA', mT), ('ADC', mA)):
    pf = np.polyfit([PX[p] for p in (4, 3, 6, 5)], mm, 1); say('%s: paddle means %s ns, slope %.4f ns/mm, v_eff %.0f mm/ns' % (lab, ' / '.join('%.2f' % v for v in mm), pf[0], 2 / pf[0]))
# per amplitude bin of the smaller of the two pulses: sigma of TOA N-S, ADC N-S, their combination, and the per-end decomposition
Amin = np.minimum(S['A'], N['A']); edges = [150, 250, 350, 450, 600, 950]
say('\nper bin of min(A_S, A_N): core sigma of N - S [ns] (TOA / ADC / mean of the two / weighted), and the three-cornered split')
rows = []
for lo, hi in zip(edges[:-1], edges[1:]):
    m = both & (Amin >= lo) & (Amin < hi) & np.isfinite(rT) & np.isfinite(rA)
    sT = core(rT[m])[1]; sA = core(rA[m])[1]; sD = core((rT - rA)[m])[1]
    # TOA and ADC errors independent: var(T-A) = vT + vA  ->  vT = (sT^2 + sD^2 - sA^2)/2
    vT = (sT ** 2 + sD ** 2 - sA ** 2) / 2; vA = (sA ** 2 + sD ** 2 - sT ** 2) / 2
    wT = 1 / max(vT, 1e-3); wA = 1 / max(vA, 1e-3); comb = (wT * rT + wA * rA) / (wT + wA); sC = core(comb[m])[1]
    rows.append((lo, hi, m.sum(), sT, sA, sD, np.sqrt(max(vT, 0)), np.sqrt(max(vA, 0)), sC))
    say('  A %3d-%3d n %6d: TOA %.2f  ADC %.2f  TOA-ADC %.2f  -> TOA-only error %.2f, ADC error %.2f (N-S); combined %.2f ns = %.0f mm' % (lo, hi, m.sum(), sT, sA, sD, np.sqrt(max(vT, 0)), np.sqrt(max(vA, 0)), sC, sC / 0.0128))
np.savez('toa_cdc2.npz', ts0=ts0, pad=pad, tS=S['t'], tN=N['t'], okS=S['ok'], okN=N['ok'], tcS=S['tc'], tcN=N['tc'], AS=S['A'], AN=N['A'], rows=np.array(rows))
c = ROOT.TCanvas('c', '', 2400, 800); c.Divide(3, 1); keep = []
for i, (P, nm) in enumerate(((S, 'South'), (N, 'North'))):
    c.cd(i + 1); fr = ROOT.gPad.DrawFrame(0, -3, 1024, 28, 'SFIL2 c3 %s;unwrapped TOA code;phase in the 25 ns frame [ns]' % nm); keep.append(fr)
    g = ROOT.TGraph(1024, np.arange(1024, dtype=float), (P['tab'] - P['tab'].min()).astype(float)); g.SetLineWidth(3); g.SetLineColor(ROOT.kRed + 1); g.Draw('l same'); keep.append(g)
    g2 = ROOT.TGraph(1024, np.arange(1024, dtype=float), np.arange(1024) * 25. / 1024); g2.SetLineStyle(2); g2.Draw('l same'); keep.append(g2)
    lg = legend(0.18, 0.76, 0.7, 0.88); lg.AddEntry(g, P['how'], 'l'); lg.AddEntry(g2, 'nominal 24.4 ps per code', 'l'); lg.Draw(); keep.append(lg)
c.cd(3); rr = np.array(rows); xm = 0.5 * (rr[:, 0] + rr[:, 1]); fr = ROOT.gPad.DrawFrame(100, 0, 1000, max(rr[:, 3].max(), rr[:, 4].max()) * 1.4, 'N - S, beam cell of SFIL2;min(A_{S}, A_{N}) [ADC];core #sigma(t_{N} - t_{S}) [ns]'); keep.append(fr)
lg = legend(0.45, 0.72, 0.9, 0.88)
for col_, mk, cl, lab in ((3, 20, ROOT.kRed + 1, 'TOA'), (4, 21, ROOT.kBlue + 1, 'ADC 50 % crossing'), (8, 22, ROOT.kBlack, 'combined')):
    g = ROOT.TGraph(len(xm), xm.astype(float), rr[:, col_].astype(float)); g.SetMarkerStyle(mk); g.SetMarkerColor(cl); g.SetLineColor(cl); g.Draw('pl same'); keep.append(g); lg.AddEntry(g, lab, 'p')
lg.Draw(); keep.append(lg)
c.SaveAs('toa_cdc2.png'); c.SaveAs('toa_cdc2.pdf'); say('wrote toa_cdc2.png')
