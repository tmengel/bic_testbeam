"""SFIL2 cell 3 KCU0 TOA as a sub-sample TDC.
Per end: pulse word = TOA word on the ADC-crossing sample k or k+1 (earliest); code unwrapped circularly at the offset c0 that
minimises the residual against the ADC crossing; code -> time by code density (uniform phase over one 25 ns frame),
t_TOA = 25 s_word + 25 CDF(code'); outliers (noise words) removed by |t_TOA - t_cross - med| < CUT; amplitude walk removed
against the ADC 50 % crossing in ln A bins.  Then N - S per paddle (x = -80/-30/+30/+80 mm) for TOA and ADC on the same showers.
    python3 toa_cdc.py [nmax] [cut_ns] -> toa_cdc.txt/.png/.npz"""
import numpy as np, uproot, sys, ROOT
sys.path.insert(0, '/sphenix/user/tmengel/eic/bic_testbeam/common'); from plotstyle import setup, legend
setup(); ROOT.gROOT.SetBatch(True); out = open('toa_cdc.txt', 'w')
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
def core(v, w=3.):
    v = v[np.isfinite(v)]; m, s = np.median(v), msig(v)
    for _ in range(5): g = np.abs(v - m) < w * s; m, s = v[g].mean(), v[g].std()
    return m, s
def pulses(tree, br, neg=False, lo=200, hi=900):
    a = f[tree].arrays([br + '_adc', br + '_toa'], entry_stop=NMAX, library='np'); W = a[br + '_adc'].astype(float); T = a[br + '_toa'].astype(int)
    n = len(W); good = W[:, 0] >= 0; y = W - W[:, 1:4].mean(axis=1)[:, None]
    if neg: y = -y
    ipk = y[:, 4:16].argmax(axis=1) + 4; A = y[np.arange(n), ipk]; half = 0.5 * A; k = ipk.copy()
    for d in range(1, 6):
        kk = np.clip(ipk - d, 0, 19); k = np.where(y[np.arange(n), kk] >= half, kk, k)
    k = np.clip(k, 1, 18); y0 = y[np.arange(n), k - 1]; y1 = y[np.arange(n), k]; tc = 25. * (k - 1 + np.clip((half - y0) / np.maximum(y1 - y0, 1e-3), 0, 1))
    s = np.full(n, -1); code = np.full(n, -1)
    for d in (0, 1):
        j = k + d; w = T[np.arange(n), j]; m = (w > 0) & (s < 0); s[m] = j[m]; code[m] = w[m]
    sel = good & (A > lo) & (A < hi) & (ipk >= 5) & (ipk <= 14)
    return dict(A=A, tc=tc, s=s, code=code, sel=sel, has=sel & (s >= 0), noise=(T[good][:, 1:4] > 0).mean())
def calibrate(P, name):
    h = P['has']; best = None
    for c0 in range(0, 1024, 16):
        cp = (P['code'][h] - c0) % 1024; occ = np.bincount(cp, minlength=1024).astype(float); cdf = np.cumsum(occ) / occ.sum(); tab = 25. * (cdf - 0.5 * occ / occ.sum())
        r = 25. * P['s'][h] + tab[cp] - P['tc'][h]; r = r - 25. * np.round((r - np.median(r)) / 25.) * 0   # no wrapping of frames
        sc = msig(r)
        if best is None or sc < best[0]: best = (sc, c0, tab)
    sc, c0, tab = best; cp = (P['code'] - c0) % 1024
    t = np.where(P['has'], 25. * P['s'] + tab[np.clip(cp, 0, 1023)], np.nan); r = t - P['tc']; med = np.nanmedian(r[P['has']])
    ok = P['has'] & (np.abs(r - med) < CUT)
    # walk vs ADC crossing in ln A bins (TOA leading edge vs 50 % crossing)
    la = np.log(np.maximum(P['A'], 1)); bins = np.linspace(np.log(200), np.log(900), 9); xb, yb = [], []
    for lo_, hi_ in zip(bins[:-1], bins[1:]):
        m = ok & (la >= lo_) & (la < hi_)
        if m.sum() > 200: xb.append(np.median(la[m])); yb.append(np.median(r[m]))
    walk = np.interp(la, xb, yb); tw = t - walk
    m_, s_ = core((tw - P['tc'])[ok])
    say('%s: pulses %d, with a word on k/k+1 %d (%.2f), pedestal word rate %.3f/sample; unwrap offset c0 = %d; TOA - crossing MAD-sigma %.2f ns (before cut), kept |r| < %.0f ns: %d (%.2f of words); after walk core sigma %.2f ns; walk span %.2f ns over A 200-900' %
        (name, P['sel'].sum(), P['has'].sum(), P['has'].sum() / P['sel'].sum(), P['noise'], c0, sc, CUT, ok.sum(), ok.sum() / P['has'].sum(), s_, max(yb) - min(yb)))
    P.update(t=tw, ok=ok, c0=c0, tab=tab, walk=(xb, yb))
    return P
S = calibrate(pulses('sfil2', 'Sc3'), 'SFIL2 S c3'); N = calibrate(pulses('sfil2', 'Nc3'), 'SFIL2 N c3')
both = S['ok'] & N['ok']; dT = N['t'] - S['t']; dA = N['tc'] - S['tc']
say('\nboth ends with a good TOA: %d of %d showers with both pulses (%.2f)' % (both.sum(), (S['sel'] & N['sel']).sum(), both.sum() / (S['sel'] & N['sel']).sum()))
res = {}
for lab, d in (('TOA', dT), ('ADC crossing', dA)):
    xs, ms, ss = [], [], []
    for p in (4, 3, 6, 5):
        m = both & (pad == PB[p]); mu, sg = core(d[m]); xs.append(PX[p]); ms.append(mu); ss.append(sg)
    pf = np.polyfit(xs, ms, 1); res[lab] = (xs, ms, ss, pf)
    say('%-12s N - S per paddle P4/P3/P6/P5: mean %s ns, core sigma %s ns; slope %.4f ns/mm -> v_eff %.0f mm/ns; position sigma per shower %s mm' %
        (lab, ' / '.join('%.2f' % v for v in ms), ' / '.join('%.2f' % v for v in ss), pf[0], 2 / pf[0], ' / '.join('%.0f' % (v / pf[0]) for v in ss)))
# correlation of the two N - S estimates after removing the paddle means
rT = np.full(len(dT), np.nan); rA = rT.copy()
for i_, p in enumerate((4, 3, 6, 5)):
    m = both & (pad == PB[p]); rT[m] = dT[m] - res['TOA'][1][i_]; rA[m] = dA[m] - res['ADC crossing'][1][i_]
g = both & np.isfinite(rT) & (np.abs(rT) < 10) & (np.abs(rA) < 10); say('residual correlation TOA vs ADC N - S: %.3f; rms of their difference %.2f ns' % (np.corrcoef(rT[g], rA[g])[0, 1], np.std(rT[g] - rA[g])))
np.savez('toa_cdc.npz', ts0=ts0, pad=pad, tS=S['t'], tN=N['t'], okS=S['ok'], okN=N['ok'], tcS=S['tc'], tcN=N['tc'], AS=S['A'], AN=N['A'], c0S=S['c0'], c0N=N['c0'], tabS=S['tab'], tabN=N['tab'])
# figure
c = ROOT.TCanvas('c', '', 2400, 1600); c.Divide(3, 2); keep = []
for i, (P, nm) in enumerate(((S, 'SFIL2 South c3'), (N, 'SFIL2 North c3'))):
    c.cd(1 + 3 * i); h = ROOT.TH2F('hc%d' % i, '%s;unwrapped TOA code;t_{cross} - 25 s_{word} [ns]' % nm, 128, 0, 1024, 100, -30, 40)
    cp = (P['code'] - P['c0']) % 1024; m = P['has']
    for a_, b_ in zip(cp[m][:300000], (P['tc'] - 25. * P['s'])[m][:300000]): h.Fill(a_, b_)
    ROOT.gPad.SetLogz(); ROOT.gPad.SetRightMargin(0.14); h.Draw('colz'); keep.append(h)
    g_ = ROOT.TGraph(1024, np.arange(1024, dtype=float), (P['tab'] + np.nanmedian((P['tc'] - P['t'])[P['ok']])).astype(float)); g_.SetLineColor(ROOT.kRed + 1); g_.SetLineWidth(3); g_.Draw('l same'); keep.append(g_)
    c.cd(2 + 3 * i); fr = ROOT.gPad.DrawFrame(-10, 0, 10, 1, '%s;t_{TOA} - t_{cross} [ns];fraction' % nm); keep.append(fr)
    hh = ROOT.TH1F('hr%d' % i, '', 160, -10, 10); [hh.Fill(v) for v in (P['t'] - P['tc'])[P['ok']] - np.median((P['t'] - P['tc'])[P['ok']])]; hh.Scale(1. / hh.GetEntries()); hh.SetLineWidth(2); hh.Draw('hist same'); fr.SetMaximum(hh.GetMaximum() * 1.3); keep.append(hh)
c.cd(3); fr = ROOT.gPad.DrawFrame(-8, 0, 8, 1, ';(t_{N} - t_{S}) - paddle mean [ns];fraction'); keep.append(fr); lg = legend(0.6, 0.72, 0.9, 0.88)
for lab, d, col in (('TOA', rT, ROOT.kRed + 1), ('ADC crossing', rA, ROOT.kBlue + 1)):
    hh = ROOT.TH1F('hd' + lab[:3], '', 160, -8, 8); [hh.Fill(v) for v in d[both & np.isfinite(d)]]; hh.Scale(1. / hh.GetEntries()); hh.SetLineColor(col); hh.SetLineWidth(2); hh.Draw('hist same'); keep.append(hh); lg.AddEntry(hh, lab, 'l')
fr.SetMaximum(max(k.GetMaximum() for k in keep[-2:]) * 1.3); lg.Draw(); keep.append(lg)
c.cd(6); xs, ms, ss, pf = res['TOA']; xa, ma, sa, pa = res['ADC crossing']; fr = ROOT.gPad.DrawFrame(-100, min(ms + ma) - 1, 100, max(ms + ma) + 1.5, ';paddle position x [mm];mean t_{N} - t_{S} [ns]'); keep.append(fr)
g1 = ROOT.TGraph(4, np.array(xs), np.array(ms)); g1.SetMarkerStyle(20); g1.SetMarkerColor(ROOT.kRed + 1); g1.Draw('p same'); keep.append(g1)
g2 = ROOT.TGraph(4, np.array(xa), np.array(ma)); g2.SetMarkerStyle(21); g2.SetMarkerColor(ROOT.kBlue + 1); g2.Draw('p same'); keep.append(g2)
lg = legend(0.18, 0.74, 0.7, 0.88); lg.AddEntry(g1, 'TOA: v_{eff} %.0f mm/ns' % (2 / pf[0]), 'p'); lg.AddEntry(g2, 'ADC crossing: v_{eff} %.0f mm/ns' % (2 / pa[0]), 'p'); lg.Draw(); keep.append(lg)
c.SaveAs('toa_cdc.png'); c.SaveAs('toa_cdc.pdf'); say('wrote toa_cdc.png')
