"""Position along the fibres from the North - South TOA difference on the channels whose TOA words sit on the pulse
(SFIL2 cell 3: ch 29 S / 173 N; SFIL4 cell 3: ch 0 S / 144 N), calibrated against the AstroPix x of the matched showers.
Valid TOA = a TOA word on samples ipk-2..ipk of a pulse > 200 ADC (earliest such word).  Each channel's TOA code is
self-calibrated to ns from its own code distribution (uniform occupancy over the 25 ns window), t = 25 s + TOA_ns.
    python3 toa_position.py [nmax]  -> toa_position.txt/.png
"""
import numpy as np, uproot, sys, ROOT
sys.path.insert(0, '/sphenix/user/tmengel/eic/bic_testbeam/common'); from plotstyle import setup, legend, label
setup(); ROOT.gROOT.SetBatch(True); out = open('toa_position.txt', 'w')
def say(*a): print(*a, flush=True); print(*a, file=out)
NMAX = int(sys.argv[1]) if len(sys.argv) > 1 else None
f = uproot.open('/sphenix/user/tmengel/eic/bic_testbeam/data/kcu0_3062.root')
ts0 = f['header']['ts0'].array(library='np', entry_stop=NMAX)
ev = uproot.open('/sphenix/user/tmengel/eic/bic_testbeam/data/events_3062.root')['events'].arrays(['trig_ts', 'pad', 'astro_x', 'astro_purity', 'n_astro'], library='np')
ets = ev['trig_ts']; order = np.argsort(ets); ets_s = ets[order]; pos = np.clip(np.searchsorted(ets_s, ts0), 0, len(ets_s) - 1); hit = ets_s[pos] == ts0; row = np.where(hit, order[pos], -1)
ax = np.full(len(ts0), np.nan); ax[hit] = ev['astro_x'][row[hit]]; pad = np.zeros(len(ts0), np.int16); pad[hit] = ev['pad'][row[hit]]; pur = np.zeros(len(ts0)); pur[hit] = ev['astro_purity'][row[hit]]
say('triggers %d, joined to the event table %d, with an AstroPix x %d (purity > 0.7: %d)' % (len(ts0), hit.sum(), np.isfinite(ax).sum(), (np.isfinite(ax) & (pur > 0.7)).sum()))
def valid_toa(tree, br):
    a = f[tree].arrays([br + '_adc', br + '_toa'], entry_stop=NMAX, library='np'); W = a[br + '_adc'].astype(float); T = a[br + '_toa']
    good = W[:, 0] >= 0; ped = W[:, 1:4].mean(axis=1); y = W - ped[:, None]; ipk = y[:, 4:13].argmax(axis=1) + 4; A = y[np.arange(len(y)), ipk]
    s = np.full(len(W), -1); code = np.full(len(W), -1)
    for d in (2, 1, 0):      # earliest first: ipk-2, then ipk-1, then ipk
        k = np.clip(ipk - d, 0, 19); w = T[np.arange(len(W)), k]; m = (w > 0) & (s < 0); s[m] = k[m]; code[m] = w[m]
    ok = good & (A > 200) & (ipk >= 6) & (ipk <= 10) & (s >= 0)
    # self-calibration: TOA_ns = 25 * CDF(code) over the valid words
    c = code[ok]; h = np.bincount(c, minlength=1024).astype(float); cdf = np.cumsum(h) / h.sum(); tns = 25. * (cdf - 0.5 * h / h.sum())   # bin-centre CDF
    t = np.where(ok, 25. * s + tns[np.clip(code, 0, 1023)], np.nan)
    return t, A, ok, h
pairs = {'SFIL2 c3': ('sfil2', 'Sc3', 'Nc3'), 'SFIL4 c3': ('sfil4', 'Sc3', 'Nc3')}
c = ROOT.TCanvas('c', '', 2400, 1600); c.Divide(3, 2); keep = []; ip = 0
for name, (tree, bS, bN) in pairs.items():
    tS, AS, okS, hS = valid_toa(tree, bS); tN, AN, okN, hN = valid_toa(tree, bN)
    both = okS & okN; dt = tN - tS
    say('\n%s: valid TOA South %d, North %d, both %d (%.2f of pulses with both ends > 200 ADC)' % (name, okS.sum(), okN.sum(), both.sum(), both.sum() / max(1, (okS | okN).sum())))
    m = both; d = dt[m]; say('  dt = tN - tS: median %.2f ns, 16-84 %% %.2f-%.2f, rms (clipped +-15) %.2f' % (np.median(d), np.percentile(d, 16), np.percentile(d, 84), d[np.abs(d - np.median(d)) < 15].std()))
    # time walk: dt against ln(AN/AS)
    la = np.log(AN[m] / AS[m]); g = np.abs(d - np.median(d)) < 15; pw = np.polyfit(la[g], d[g], 1); say('  walk: dt vs ln(AN/AS) slope %.2f ns per unit, intercept %.2f' % (pw[0], pw[1]))
    dtc = dt - pw[0] * np.log(np.maximum(AN, 1) / np.maximum(AS, 1))
    # against AstroPix x
    mx = m & np.isfinite(ax) & (pur > 0.7); x = ax[mx]
    for lab, dd in (('raw', dt[mx]), ('walk-corrected', dtc[mx])):
        g = np.abs(dd - np.median(dd)) < 15; p = np.polyfit(x[g], dd[g], 1); res = dd[g] - np.polyval(p, x[g]); r = np.corrcoef(x[g], dd[g])[0, 1]
        veff = 2. / p[0] if p[0] != 0 else np.nan
        say('  vs AstroPix x (%d showers, purity > 0.7), %s: slope %.4f ns/mm -> v_eff = 2/slope = %.0f mm/ns; intercept %.2f ns at x = 0; corr %.3f; residual rms %.2f ns -> position rms %.0f mm' % (g.sum(), lab, p[0], veff, p[1], r, res.std(), res.std() / abs(p[0]) / 2 if p[0] else np.nan))
    # per paddle medians (paddle x: P4 -80, P3 -30, P6 +30, P5 +80)
    PB = {4: 2, 3: 1, 6: 8, 5: 4}; PX = {4: -80, 3: -30, 6: 30, 5: 80}
    say('  per paddle median walk-corrected dt: ' + '  '.join('P%d (x %+d) %.2f ns (n %d)' % (p_, PX[p_], np.median(dtc[m & (pad == PB[p_])]), (m & (pad == PB[p_])).sum()) for p_ in (4, 3, 6, 5)))
    # profile in x bins
    edges = np.arange(-100, 101, 10); xc = 0.5 * (edges[1:] + edges[:-1]); prof = []; perr = []
    for lo, hi in zip(edges[:-1], edges[1:]):
        s_ = mx & (ax >= lo) & (ax < hi); v = dtc[s_]; v = v[np.abs(v - np.median(dtc[mx])) < 15] if s_.sum() > 20 else v
        prof.append(np.median(v) if len(v) > 20 else np.nan); perr.append(1.2533 * v.std() / np.sqrt(len(v)) if len(v) > 20 else np.nan)
    prof, perr = np.array(prof), np.array(perr)
    # figure: dt distribution; dt vs x scatter + profile; TOA code histograms
    c.cd(1 + 3 * ip); fr = ROOT.gPad.DrawFrame(-15, 0, 15, 1.0, '%s;t_{N} - t_{S} [ns];fraction' % name); keep.append(fr)
    h1 = ROOT.TH1F('hd%d' % ip, '', 120, -15, 15); [h1.Fill(v) for v in dtc[m]]; h1.Scale(1. / max(1, h1.GetEntries())); h1.SetLineWidth(2); h1.Draw('hist same'); fr.SetMaximum(h1.GetMaximum() * 1.3); keep.append(h1)
    c.cd(2 + 3 * ip); fr2 = ROOT.gPad.DrawFrame(-100, -12, 100, 12, '%s;AstroPix x [mm];t_{N} - t_{S} (walk-corrected) [ns]' % name); keep.append(fr2)
    h2 = ROOT.TH2F('hx%d' % ip, '', 40, -100, 100, 96, -12, 12); [h2.Fill(a_, b_) for a_, b_ in zip(x, dtc[mx])]; h2.Draw('col same'); keep.append(h2)
    fin = np.isfinite(prof); gp = ROOT.TGraphErrors(int(fin.sum()), xc[fin].astype(float), prof[fin], np.zeros(int(fin.sum())), perr[fin]); gp.SetMarkerStyle(20); gp.SetMarkerColor(ROOT.kRed + 1); gp.SetLineColor(ROOT.kRed + 1); gp.Draw('p same'); keep.append(gp)
    c.cd(3 + 3 * ip); fr3 = ROOT.gPad.DrawFrame(0, 0, 1024, max(hS.max(), hN.max()) * 1.2, '%s;TOA code;valid words' % name); keep.append(fr3)
    for hh, colr, nm in ((hS, ROOT.kBlue + 1, 'South ' + bS), (hN, ROOT.kRed + 1, 'North ' + bN)):
        g_ = ROOT.TGraph(1024, np.arange(1024, dtype=float), hh.astype(float)); g_.SetLineColor(colr); g_.Draw('l same'); keep.append(g_)
    lg = legend(0.5, 0.7, 0.88, 0.88); lg.AddEntry(keep[-2], 'South', 'l'); lg.AddEntry(keep[-1], 'North', 'l'); lg.Draw(); keep.append(lg)
    ip += 1
c.SaveAs('toa_position.png'); c.SaveAs('toa_position.pdf'); say('wrote toa_position.png')
