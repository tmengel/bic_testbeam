"""KCU0 TOA code calibration against the ADC rising edge, then the North - South time difference per paddle.
The TOA code occupancy on KCU0 is not uniform over the window (the codes bunch), so the self-calibration by occupancy
used on KCU1 is invalid here.  Instead, per channel, f(code) = median of (t_cross - 25 s_TOA) in code bins, where t_cross is
the 50 % crossing time from a linear interpolation of the two ADC samples around the crossing; t_TOA = 25 s_TOA + f(code).
The N - S difference per paddle (x = -80 / -30 / +30 / +80 mm) then gives 2/v_eff and the fixed N - S delay.
    python3 toa_adc_calib.py [nmax] -> toa_adc_calib.txt/.png
"""
import numpy as np, uproot, sys, ROOT
sys.path.insert(0, '/sphenix/user/tmengel/eic/bic_testbeam/common'); from plotstyle import setup, legend, label
setup(); ROOT.gROOT.SetBatch(True); out = open('toa_adc_calib.txt', 'w')
def say(*a): print(*a, flush=True); print(*a, file=out)
NMAX = int(sys.argv[1]) if len(sys.argv) > 1 else None
f = uproot.open('/sphenix/user/tmengel/eic/bic_testbeam/data/kcu0_3062.root')
ts0 = f['header']['ts0'].array(library='np', entry_stop=NMAX)
ev = uproot.open('/sphenix/user/tmengel/eic/bic_testbeam/data/events_3062.root')['events'].arrays(['trig_ts', 'pad', 'astro_x', 'astro_purity'], library='np')
ets = ev['trig_ts']; order = np.argsort(ets); ets_s = ets[order]; pos = np.clip(np.searchsorted(ets_s, ts0), 0, len(ets_s) - 1); hit = ets_s[pos] == ts0; row = np.where(hit, order[pos], -1)
pad = np.zeros(len(ts0), np.int16); pad[hit] = ev['pad'][row[hit]]; ax = np.full(len(ts0), np.nan); ax[hit] = ev['astro_x'][row[hit]]; pur = np.zeros(len(ts0)); pur[hit] = ev['astro_purity'][row[hit]]
PB = {4: 2, 3: 1, 6: 8, 5: 4}; PX = {4: -80., 3: -30., 6: 30., 5: 80.}
def channel(tree, br):
    a = f[tree].arrays([br + '_adc', br + '_toa'], entry_stop=NMAX, library='np'); W = a[br + '_adc'].astype(float); T = a[br + '_toa']
    good = W[:, 0] >= 0; ped = W[:, 1:4].mean(axis=1); y = W - ped[:, None]; ipk = y[:, 4:13].argmax(axis=1) + 4; A = y[np.arange(len(y)), ipk]
    # 50 % crossing by linear interpolation: first sample k <= ipk with y[k] >= A/2, between k-1 and k
    half = 0.5 * A; k = ipk.copy()
    for d in range(1, 5):
        kk = np.clip(ipk - d, 0, 19); below = y[np.arange(len(y)), kk] >= half; k = np.where(below, kk, k)
    k = np.clip(k, 1, 19); y0 = y[np.arange(len(y)), k - 1]; y1 = y[np.arange(len(y)), k]; frac = np.clip((half - y0) / np.maximum(y1 - y0, 1e-3), 0, 1); tcross = 25. * (k - 1 + frac)
    s = np.full(len(W), -1); code = np.full(len(W), -1)
    for d in (2, 1, 0):
        kk = np.clip(ipk - d, 0, 19); w = T[np.arange(len(W)), kk]; m = (w > 0) & (s < 0); s[m] = kk[m]; code[m] = w[m]
    ok = good & (A > 200) & (ipk >= 6) & (ipk <= 10) & (s >= 0) & (A < 900)
    # calibration table: median (tcross - 25 s) per code bin of 8
    nb = 128; cb = np.minimum(code // 8, nb - 1); tab = np.full(nb, np.nan); cnt = np.zeros(nb, int)
    for b in range(nb):
        m = ok & (cb == b); cnt[b] = m.sum()
        if cnt[b] >= 50: tab[b] = np.median(tcross[m] - 25. * s[m])
    fin = np.isfinite(tab); tabi = np.interp(np.arange(nb), np.where(fin)[0], tab[fin]) if fin.sum() > 2 else tab
    ttoa = np.where(ok, 25. * s + tabi[np.clip(cb, 0, nb - 1)], np.nan)
    resid = (ttoa - tcross)[ok]
    return ttoa, tcross, A, ok, tabi, cnt, resid
c = ROOT.TCanvas('c', '', 2400, 1600); c.Divide(3, 2); keep = []; ip = 0
for name, (tree, bS, bN) in {'SFIL2 c3': ('sfil2', 'Sc3', 'Nc3'), 'SFIL4 c3': ('sfil4', 'Sc3', 'Nc3')}.items():
    tS, cS, AS, okS, tabS, cntS, rS = channel(tree, bS); tN, cN, AN, okN, tabN, cntN, rN = channel(tree, bN)
    say('\n%s: calibrated against the ADC 50 %% crossing; TOA - crossing residual rms South %.2f ns, North %.2f ns (over %d / %d pulses); table span South %.1f-%.1f ns, North %.1f-%.1f ns' % (name, rS[np.abs(rS) < 15].std(), rN[np.abs(rN) < 15].std(), okS.sum(), okN.sum(), np.nanmin(tabS), np.nanmax(tabS), np.nanmin(tabN), np.nanmax(tabN)))
    both = okS & okN; dt = tN - tS; dc = cN - cS
    say('  both valid %d; N - S from TOA: median %.2f, clipped rms %.2f ns; from the ADC crossings alone: median %.2f, rms %.2f ns' % (both.sum(), np.median(dt[both]), dt[both][np.abs(dt[both] - np.median(dt[both])) < 15].std(), np.median(dc[both]), dc[both][np.abs(dc[both] - np.median(dc[both])) < 15].std()))
    xs, ms_toa, ms_adc, es = [], [], [], []
    for p_ in (4, 3, 6, 5):
        m = both & (pad == PB[p_]); v = dt[m]; v = v[np.abs(v - np.median(dt[both])) < 15]; w = dc[m]; w = w[np.abs(w - np.median(dc[both])) < 15]
        xs.append(PX[p_]); ms_toa.append(np.median(v)); ms_adc.append(np.median(w)); es.append(1.25 * v.std() / np.sqrt(len(v)))
        say('  P%d (x %+.0f): n %6d  N - S TOA median %.3f +- %.3f ns (rms %.2f)   ADC-crossing median %.3f ns (rms %.2f)' % (p_, PX[p_], m.sum(), np.median(v), 1.25 * v.std() / np.sqrt(len(v)), v.std(), np.median(w), w.std()))
    for lab, ms in (('TOA', ms_toa), ('ADC crossing', ms_adc)):
        p = np.polyfit(xs, ms, 1); say('  %s: slope %.4f ns/mm -> v_eff = 2/slope = %.0f mm/ns (expected 160-200 for scintillating fibres); offset at x = 0: %.2f ns' % (lab, p[0], 2 / p[0] if p[0] else np.nan, p[1]))
    # AstroPix x (whatever its coverage) for the continuous check
    mx = both & np.isfinite(ax) & (pur > 0.3); v = dt[mx]; g = np.abs(v - np.median(dt[both])) < 15
    if g.sum() > 200: p = np.polyfit(ax[mx][g], v[g], 1); say('  vs AstroPix x (%d showers, purity > 0.3, x %.0f-%.0f mm): slope %.4f ns/mm, corr %.3f' % (g.sum(), np.percentile(ax[mx], 5), np.percentile(ax[mx], 95), p[0], np.corrcoef(ax[mx][g], v[g])[0, 1]))
    c.cd(1 + 3 * ip); fr = ROOT.gPad.DrawFrame(0, -30, 1024, 30, '%s;TOA code;t_{cross} - 25 s_{TOA} [ns]' % name); keep.append(fr)
    for tab, colr, nm in ((tabS, ROOT.kBlue + 1, 'South'), (tabN, ROOT.kRed + 1, 'North')):
        g_ = ROOT.TGraph(128, (np.arange(128) * 8 + 4).astype(float), tab.astype(float)); g_.SetLineColor(colr); g_.SetLineWidth(2); g_.Draw('l same'); keep.append(g_)
    lg = legend(0.15, 0.72, 0.5, 0.88); lg.AddEntry(keep[-2], 'South', 'l'); lg.AddEntry(keep[-1], 'North', 'l'); lg.Draw(); keep.append(lg)
    c.cd(2 + 3 * ip); fr2 = ROOT.gPad.DrawFrame(-15, 0, 15, 0.1, '%s;t_{N} - t_{S} [ns];fraction' % name); keep.append(fr2); lg2 = legend(0.6, 0.65, 0.88, 0.88)
    for p_, colr in zip((4, 3, 6, 5), (ROOT.kBlue + 1, ROOT.kGreen + 2, ROOT.kOrange + 7, ROOT.kRed + 1)):
        m = both & (pad == PB[p_]); h = ROOT.TH1F('h%d%d' % (ip, p_), '', 120, -15, 15); [h.Fill(v) for v in dt[m] - np.median(dt[both])]; h.Scale(1. / max(1, h.GetEntries())); h.SetLineColor(colr); h.SetLineWidth(2); h.Draw('hist same'); keep.append(h); lg2.AddEntry(h, 'P%d (x %+.0f mm)' % (p_, PX[p_]), 'l')
    fr2.SetMaximum(max(k.GetMaximum() for k in keep[-4:]) * 1.3); lg2.Draw(); keep.append(lg2)
    c.cd(3 + 3 * ip); fr3 = ROOT.gPad.DrawFrame(-100, min(ms_toa + ms_adc) - 1, 100, max(ms_toa + ms_adc) + 1, '%s;paddle position x [mm];median t_{N} - t_{S} [ns]' % name); keep.append(fr3)
    g1 = ROOT.TGraphErrors(4, np.array(xs), np.array(ms_toa), np.zeros(4), np.array(es)); g1.SetMarkerStyle(20); g1.SetMarkerColor(ROOT.kRed + 1); g1.SetLineColor(ROOT.kRed + 1); g1.Draw('p same'); keep.append(g1)
    g2 = ROOT.TGraph(4, np.array(xs), np.array(ms_adc)); g2.SetMarkerStyle(24); g2.SetMarkerColor(ROOT.kBlue + 1); g2.Draw('p same'); keep.append(g2)
    lg3 = legend(0.15, 0.72, 0.6, 0.88); lg3.AddEntry(g1, 'TOA (ADC-calibrated codes)', 'p'); lg3.AddEntry(g2, 'ADC crossing only', 'p'); lg3.Draw(); keep.append(lg3)
    ip += 1
c.SaveAs('toa_adc_calib.png'); c.SaveAs('toa_adc_calib.pdf'); say('wrote toa_adc_calib.png')
