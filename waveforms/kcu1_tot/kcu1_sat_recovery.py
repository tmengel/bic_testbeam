"""SFIL5 South (KCU1 arrays B, C, D) amplitude recovery for saturated single-SiPM pulses.
Two handles, both data-only:
  * neighbour predictor P_k: the 16 SiPMs of an array view the same cell, so in unsaturated events the pair ratios
    r_jk = median(A_j / A_k) are stable; for a saturated channel k, P_k = median_j (A_j / r_jk) over the unsaturated
    channels j of the same array.  Validated on unsaturated channels (P_k / A_k) and used where there is no TOT.
  * TOT: for saturated pulses with a TOT word, the array-pooled profile F(TOT) = median P over TOT bins, scaled per
    channel (s_k), gives A(TOT).  The nsat = 0 pulses that carry a TOT word (A_fit 650-943) check the onset of F.
The template extrapolation (afit for nsat > 0) is NOT used: it returns less than the clipped peak.
    python3 kcu1_sat_recovery.py  -> kcu1_sat_recovery.txt, sfil5_recovered.npz (per events_3062 row), figures
"""
import numpy as np, uproot, sys, ROOT
sys.path.insert(0, '/sphenix/user/tmengel/eic/bic_testbeam/common'); from plotstyle import setup, legend, label
setup()
out = open('kcu1_sat_recovery.txt', 'w')
def say(*a): print(*a, flush=True); print(*a, file=out)
D = np.load('/sphenix/tg/tg01/jets/tmengel/eic/kcu1/fit_full/kcu1_fit_arrays.npz')
afit, amp, nsat, tot, ok, ts0 = D['afit'].astype(np.float32), D['amp'].astype(np.float32), D['nsat'], D['tot'], D['ok'], D['ts0']
n = len(ts0); CEIL = 943.
live = ok.sum(axis=0) > 20000
say('triggers %d; channels with > 20k good fits: %d (%s)' % (n, live.sum(), ','.join(str(k) for k in np.where(live)[0])))
ARR = {'B': range(16, 32), 'C': range(32, 48), 'D': range(48, 64)}
unsat = ok & (nsat == 0)
P = np.full((n, 64), np.nan, np.float32); NV = np.zeros((n, 64), np.int8)
R = np.full((64, 64), np.nan)      # R[k, j] = E[A_k | A_j] / A_j, slope through the origin, from the A_j range where k saturates in < 2 % of triggers
satk = nsat > 0; present = amp > -1
for name, ks in ARR.items():
    ks = [k for k in ks if live[k]]
    for k in ks:
        for j in ks:
            if j == k: continue
            base_j = unsat[:, j] & (afit[:, j] > 60) & present[:, k]
            # upper edge of the A_j range: highest 40-ADC bin where fewer than 2 % of the triggers have k saturated
            hi = 60
            for e in range(100, 960, 40):
                m = base_j & (afit[:, j] >= e - 40) & (afit[:, j] < e)
                if m.sum() < 300 or satk[m, k].mean() > 0.02: break
                hi = e
            m = base_j & (afit[:, j] < hi) & (afit[:, j] >= 60) & (unsat[:, k] | satk[:, k])
            mm = m & unsat[:, k]            # k measured; the saturated ones (< 2 %) are dropped, a negligible truncation
            if mm.sum() > 500: R[k, j] = afit[mm, k].mean() / afit[mm, j].mean()
    M = np.where(unsat[:, ks] & (afit[:, ks] > 30), afit[:, ks], np.nan).astype(np.float32)
    for k in ks:
        r = np.array([R[k, j] for j in ks], np.float32); Q = M * r[None, :]; Q[:, ks.index(k)] = np.nan; Q[:, np.isnan(r)] = np.nan
        NV[:, k] = np.sum(~np.isnan(Q), axis=1); P[:, k] = np.nanmedian(Q, axis=1)
    say('array %s: live channels %d; regression slopes defined for %d pairs; slope range %.2f-%.2f' % (name, len(ks), np.isfinite(R[np.ix_(ks, ks)]).sum(), np.nanmin(R[np.ix_(ks, ks)]), np.nanmax(R[np.ix_(ks, ks)])))
valid = NV >= 6
# ---- validation on unsaturated channels
say('\nPredictor check: median and 16-84 %% spread of A_fit,k / P_k in bins of P_k (unsaturated k; above P ~ 700 the saturated k are missing, so the ratio must drop)')
edges = [100, 200, 300, 400, 500, 600, 700, 800, 900, 1000, 1200, 1500, 2000]
for lo, hi in zip(edges[:-1], edges[1:]):
    m = unsat & valid & (P >= lo) & (P < hi) & live[None, :]; ms = satk & valid & (P >= lo) & (P < hi) & live[None, :]
    q = afit[m] / P[m]; say('  P %4d-%4d: n %8d (plus %7d saturated, %.1f %%)  median %.3f  16%% %.3f  84%% %.3f' % (lo, hi, m.sum(), ms.sum(), 100. * ms.sum() / max(1, m.sum() + ms.sum()), np.median(q), np.percentile(q, 16), np.percentile(q, 84)))
# per-channel bias of the predictor in 400-900 (should be ~1 by construction of the medians)
say('per-channel median A/P (P 300-600): ' + ' '.join('%d:%.2f' % (k, np.median(afit[unsat[:, k] & valid[:, k] & (P[:, k] > 300) & (P[:, k] < 600), k] / P[unsat[:, k] & valid[:, k] & (P[:, k] > 300) & (P[:, k] < 600), k])) for k in range(64) if live[k]))
# ---- TOT calibration: array-pooled F(TOT) from saturated pulses with a predictor, scaled per channel
tb = np.arange(0, 4200, 100); tc = 0.5 * (tb[1:] + tb[:-1])
sat = (nsat > 0); hastot = tot > 0
F = {}; S = np.ones(64); Fcnt = {}
say('\nTOT calibration: array-pooled median predictor per TOT bin (saturated pulses with a TOT word and >= 6 unsaturated neighbours)')
for name, ks in ARR.items():
    ks = [k for k in ks if live[k]]
    m = sat[:, ks] & hastot[:, ks] & valid[:, ks] & (tot[:, ks] < 4092); T = tot[:, ks][m].astype(float); A = P[:, ks][m]
    prof = np.full(len(tc), np.nan); cnt = np.zeros(len(tc), int)
    for b in range(len(tc)):
        s = (T >= tb[b]) & (T < tb[b + 1]); cnt[b] = s.sum()
        if cnt[b] >= 30: prof[b] = np.median(A[s])
    F[name] = prof; Fcnt[name] = cnt
    say('  array %s: %d saturated pulses with TOT and predictor; TOT bins populated: %s' % (name, m.sum(), ', '.join('%d:%.0f(%d)' % (tc[b], prof[b], cnt[b]) for b in range(len(tc)) if cnt[b] >= 30)))
    # per-channel scale: median of P / F(TOT) over that channel's saturated pulses
    for k in ks:
        mk = sat[:, k] & hastot[:, k] & valid[:, k]; T = tot[mk, k].astype(float); A = P[mk, k]
        fb = np.interp(T, tc[np.isfinite(prof)], prof[np.isfinite(prof)]); good = np.isfinite(fb) & (T >= tb[np.where(cnt >= 30)[0][0]]) & (T < tb[np.where(cnt >= 30)[0][-1] + 1])
        S[k] = np.median(A[good] / fb[good]) if good.sum() > 50 else 1.
    say('    per-channel scale s_k: ' + ' '.join('%d:%.2f' % (k, S[k]) for k in ks))
    # onset check: nsat = 0 pulses with a TOT word: median A_fit per TOT bin vs F
    m0 = unsat[:, ks] & hastot[:, ks]; T0 = tot[:, ks][m0].astype(float); A0 = afit[:, ks][m0]
    say('    onset check (nsat = 0 pulses with TOT, A_fit vs F): ' + ', '.join('TOT %d: A_fit %.0f (n %d) vs F %.0f' % (tc[b], np.median(A0[(T0 >= tb[b]) & (T0 < tb[b + 1])]), ((T0 >= tb[b]) & (T0 < tb[b + 1])).sum(), prof[b]) for b in range(len(tc)) if ((T0 >= tb[b]) & (T0 < tb[b + 1])).sum() >= 30 and np.isfinite(prof[b])))
    say('    TOT at the decoder maximum (>= 4092): %d of %d saturated-with-TOT pulses' % ((tot[:, ks][m] >= 4092).sum(), m.sum()))
# ---- rebuild
A = np.zeros((n, 64), np.float32); src = np.zeros((n, 64), np.int8)   # 0 none, 1 fit, 2 raw peak, 3 TOT, 4 predictor, 5 floor
for name, ks in ARR.items():
    ks = [k for k in ks if live[k]]; prof = F[name]; fin = np.isfinite(prof); lo_b, hi_b = tb[np.where(fin)[0][0]], tb[np.where(fin)[0][-1] + 1]
    for k in ks:
        f = unsat[:, k]; A[f, k] = afit[f, k]; src[f, k] = 1
        r = ~ok[:, k] & (nsat[:, k] == 0) & (amp[:, k] > 60); A[r, k] = amp[r, k]; src[r, k] = 2
        s = sat[:, k]
        t = s & hastot[:, k] & (tot[:, k] >= lo_b) & (tot[:, k] < hi_b); A[t, k] = S[k] * np.interp(tot[t, k].astype(float), tc[fin], prof[fin]); src[t, k] = 3
        p = s & ~t & valid[:, k]; A[p, k] = np.maximum(P[p, k], CEIL); src[p, k] = 4
        fl = s & ~t & ~valid[:, k]; A[fl, k] = np.maximum(amp[fl, k], CEIL); src[fl, k] = 5
ks_all = [k for ks in ARR.values() for k in ks if live[k]]
A5 = A[:, ks_all].sum(axis=1); A5old = np.where(unsat[:, ks_all], afit[:, ks_all], np.where(amp[:, ks_all] > 60, amp[:, ks_all], 0)).sum(axis=1)
nsat5 = sat[:, ks_all].sum(axis=1); ntot = (src[:, ks_all] == 3).sum(axis=1); npred = (src[:, ks_all] == 4).sum(axis=1); nfloor = (src[:, ks_all] == 5).sum(axis=1)
say('\nRebuild: channel-pulses by source over live B+C+D channels: fit %d, raw peak %d, TOT %d, predictor %d, floor %d' % tuple((src[:, ks_all] == v).sum() for v in (1, 2, 3, 4, 5)))
say('per-trigger: mean saturated channels %.2f; mean A5 old %.0f -> recovered %.0f (ratio %.3f); for triggers with >= 4 saturated: %.0f -> %.0f (%.3f)' % (nsat5.mean(), A5old.mean(), A5.mean(), A5.mean() / A5old.mean(), A5old[nsat5 >= 4].mean(), A5[nsat5 >= 4].mean(), A5[nsat5 >= 4].mean() / A5old[nsat5 >= 4].mean()))
# ---- join to the event table: KCU1 ts0 -> KCU0 clock (+1143665946, +0/+4 ticks) -> events trig_ts
ev = uproot.open('/sphenix/user/tmengel/eic/bic_testbeam/data/events_3062.root')['events'].arrays(['trig_ts', 'has_kcu1', 'sfil5_s'], library='np')
ets = ev['trig_ts']; order = np.argsort(ets); ets_s = ets[order]
row = np.full(n, -1, np.int64)
for off in (0, 4):
    t = ts0 + 1143665946 + off; pos = np.clip(np.searchsorted(ets_s, t), 0, len(ets_s) - 1); hit = (ets_s[pos] == t) & (row < 0); row[hit] = order[pos[hit]]
    say('join offset +%d: matched %d' % (off, hit.sum()))
hit = row >= 0
s5 = np.clip(ev['sfil5_s'][row[hit]].astype(float), 0, None); A5tab = s5[:, 18:35].sum(axis=1) + s5[:, 36:53].sum(axis=1) + s5[:, 54:71].sum(axis=1)
say('matched %d of %d fit triggers; event-table raw A5 vs old fit-sum A5 mean ratio %.3f' % (hit.sum(), n, A5tab.mean() / A5old[hit].mean()))
np.savez('/sphenix/user/tmengel/eic/bic_testbeam/sim/analysis/sfil5_recovered.npz', row=row[hit], A5=A5[hit], A5old=A5old[hit], A5tab=A5tab, nsat5=nsat5[hit], ntot=ntot[hit], npred=npred[hit], nfloor=nfloor[hit], ts0=ts0[hit])
np.savez('kcu1_tot_calibration.npz', tc=tc, F_B=F['B'], F_C=F['C'], F_D=F['D'], S=S, R=R)
say('wrote sfil5_recovered.npz (sim/analysis) and kcu1_tot_calibration.npz')
# ---- figures
c = ROOT.TCanvas('c', '', 2400, 1600); c.Divide(3, 2); keep = []
c.cd(1); fr = ROOT.gPad.DrawFrame(0, 0, 1500, 2.5, ';P_{k} (neighbour prediction) [ADC];A_{fit} / P_{k}'); keep.append(fr)
h2 = ROOT.TH2F('h2', '', 75, 0, 1500, 50, 0, 2.5); m = unsat & valid & live[None, :] & (afit > 60); [h2.Fill(x, y) for x, y in zip(P[m][::50], (afit[m] / P[m])[::50])]; h2.Draw('col same'); keep.append(h2); keep.append(label('fitted amplitude over neighbour prediction, unsaturated pulses', x=0.15, y=0.93, size=0.035))
for p, (name, col) in enumerate(zip('BCD', (ROOT.kRed + 1, ROOT.kBlue + 1, ROOT.kGreen + 2))):
    c.cd(2 + p); fr = ROOT.gPad.DrawFrame(0, 0, 4200, 5000, 'array %s;decoded TOT;amplitude [ADC]' % name); keep.append(fr)
    ks = [k for k in ARR[name] if live[k]]; m = sat[:, ks] & hastot[:, ks] & valid[:, ks]; T = tot[:, ks][m].astype(float); Am = P[:, ks][m]
    hh = ROOT.TH2F('h%s' % name, '', 84, 0, 4200, 100, 0, 5000); [hh.Fill(x, y) for x, y in zip(T[::20], Am[::20])]; hh.Draw('col same'); keep.append(hh)
    fin = np.isfinite(F[name]); g = ROOT.TGraph(int(fin.sum()), tc[fin].astype(float), F[name][fin].astype(float)); g.SetLineColor(ROOT.kBlack); g.SetLineWidth(3); g.Draw('l same'); keep.append(g)
    m0 = unsat[:, ks] & hastot[:, ks]; T0 = tot[:, ks][m0].astype(float); A0 = afit[:, ks][m0]
    g0 = ROOT.TGraph(len(T0[::5]), T0[::5].astype(float), A0[::5].astype(float)); g0.SetMarkerColor(ROOT.kRed); g0.SetMarkerStyle(20); g0.SetMarkerSize(0.3); g0.Draw('p same'); keep.append(g0)
    lg = legend(0.5, 0.15, 0.88, 0.32); lg.AddEntry(hh, 'saturated: predictor P_{k}', 'f'); lg.AddEntry(g, 'F(TOT) = median P per bin', 'l'); lg.AddEntry(g0, 'unsaturated with TOT: A_{fit}', 'p'); lg.Draw(); keep.append(lg)
c.cd(5); fr = ROOT.gPad.DrawFrame(0, 0, 25000, 40000, ';A5 old (fit + raw) [ADC];A5 recovered [ADC]'); keep.append(fr)
h5 = ROOT.TH2F('h5', '', 100, 0, 25000, 100, 0, 40000); [h5.Fill(x, y) for x, y in zip(A5old[::20], A5[::20])]; h5.Draw('col same'); keep.append(h5); ln = ROOT.TLine(0, 0, 25000, 25000); ln.SetLineStyle(2); ln.Draw(); keep.append(ln)
c.cd(6); fr = ROOT.gPad.DrawFrame(0, 0.5, 20, 3, ';saturated channels in the trigger;<A5 recovered / A5 old>'); keep.append(fr)
xs = np.arange(0, 20); ys = np.array([A5[nsat5 == x].mean() / A5old[nsat5 == x].mean() if (nsat5 == x).sum() > 20 else np.nan for x in xs]); f = np.isfinite(ys)
g = ROOT.TGraph(int(f.sum()), xs[f].astype(float), ys[f]); g.SetMarkerStyle(20); g.Draw('p same'); keep.append(g)
c.SaveAs('kcu1_sat_recovery.png'); c.SaveAs('kcu1_sat_recovery.pdf'); say('wrote kcu1_sat_recovery.png')
