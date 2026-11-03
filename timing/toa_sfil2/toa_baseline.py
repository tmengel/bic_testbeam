"""SFIL2 cell 3 TOA after the toa_cdc2 calibration: correct the per-trigger baseline offset (the TOA threshold acts on the
unsubtracted input) with a table in (baseline, amplitude) built against the ADC 50 % crossing on even triggers, applied to odd
triggers; then N - S per paddle for TOA, ADC and their inverse-variance combination, per amplitude.
    python3 toa_baseline.py -> toa_baseline.txt/.png"""
import numpy as np, uproot, sys, ROOT
sys.path.insert(0, '/sphenix/user/tmengel/eic/bic_testbeam/common'); from plotstyle import setup, legend
setup(); ROOT.gROOT.SetBatch(True); out = open('toa_baseline.txt', 'w')
def say(*a): print(*a, flush=True); print(*a, file=out)
z = np.load('toa_cdc2.npz'); NT = len(z['ts0']); pad = z['pad']
f = uproot.open('/sphenix/user/tmengel/eic/bic_testbeam/data/kcu0_3062.root')
PB = {4: 2, 3: 1, 6: 8, 5: 4}; PX = {4: -80., 3: -30., 6: 30., 5: 80.}
def core(v, w=2.5):
    v = v[np.isfinite(v)]
    if len(v) < 30: return np.nan, np.nan
    m = np.median(v); s = 1.4826 * np.median(np.abs(v - m))
    for _ in range(6): g = np.abs(v - m) < w * s; m, s = v[g].mean(), v[g].std()
    return m, s / 0.9
even = (np.arange(NT) % 2) == 0
BE = np.array([-200, -30, -20, -12, -6, 0, 6, 12, 20, 30, 200]); AE = np.array([150, 250, 350, 450, 600, 950])
T = {}
for end, br in (('S', 'Sc3'), ('N', 'Nc3')):
    W = f['sfil2'].arrays([br + '_adc'], entry_stop=NT, library='np')[br + '_adc'].astype(float)
    ped = W[:, 1:4].mean(axis=1); base = ped - np.median(ped[W[:, 0] >= 0]); A = z['A' + end]; ok = z['ok' + end]
    r = z['t' + end] - z['tc' + end]; r0 = core(r[ok])[0]; r = r - r0
    ib = np.clip(np.digitize(base, BE) - 1, 0, len(BE) - 2); ia = np.clip(np.digitize(A, AE) - 1, 0, len(AE) - 2)
    tab = np.zeros((len(BE) - 1, len(AE) - 1))
    for i in range(len(BE) - 1):
        for j in range(len(AE) - 1):
            m = ok & even & (ib == i) & (ia == j); tab[i, j] = core(r[m])[0] if m.sum() > 100 else np.nan
    # fill empty cells from the row/column mean
    for i in range(tab.shape[0]):
        if np.all(np.isnan(tab[i])): tab[i] = np.nanmean(tab, axis=0)
        tab[i] = np.where(np.isnan(tab[i]), np.nanmean(tab[i]), tab[i])
    corr = tab[ib, ia]; t = z['t' + end] - r0 - corr
    m_o = ok & ~even
    say('%s: baseline rms %.1f ADC; odd triggers: TOA - crossing core sigma %.2f ns before, %.2f ns after the (baseline, A) table' % (end, base[ok].std(), core(r[m_o])[1], core((t - z['tc' + end])[m_o])[1]))
    say('   table (rows baseline %s, columns A %s) [ns]:' % (list(BE), list(AE)))
    for i in range(tab.shape[0]): say('     ' + ' '.join('%6.2f' % v for v in tab[i]))
    T[end] = dict(t=t, ok=ok, tc=z['tc' + end], A=A, base=base)
both = T['S']['ok'] & T['N']['ok'] & ~even
dT = T['N']['t'] - T['S']['t']; dA = T['N']['tc'] - T['S']['tc']
def resid(d, m):
    r = np.full(len(d), np.nan); mus = []
    for p in (4, 3, 6, 5):
        mm = m & (pad == PB[p]); mu = core(d[mm])[0]; r[mm] = d[mm] - mu; mus.append(mu)
    return r, mus
rT, mT = resid(dT, both); rA, mA = resid(dA, both)
for lab, mm in (('TOA', mT), ('ADC', mA)):
    pf = np.polyfit([PX[p] for p in (4, 3, 6, 5)], mm, 1); say('%s N - S paddle means %s ns; slope %.4f ns/mm -> v_eff %.0f mm/ns' % (lab, ' / '.join('%.2f' % v for v in mm), pf[0], 2 / pf[0]))
pA = np.polyfit([PX[p] for p in (4, 3, 6, 5)], mA, 1)[0]
Amin = np.minimum(T['S']['A'], T['N']['A']); rows = []
say('\nodd triggers, both ends good; per bin of min(A_S, A_N): core sigma of N - S [ns]; per-estimator errors from the three-cornered split; inverse-variance combination')
for lo, hi in zip(AE[:-1], AE[1:]):
    m = both & (Amin >= lo) & (Amin < hi) & np.isfinite(rT) & np.isfinite(rA)
    sT = core(rT[m])[1]; sA = core(rA[m])[1]; sD = core((rT - rA)[m])[1]
    vT = max((sT ** 2 + sD ** 2 - sA ** 2) / 2, 1e-3); vA = max((sA ** 2 + sD ** 2 - sT ** 2) / 2, 1e-3); vC = max(sA ** 2 - vA, 0)
    comb = (rT / vT + rA / vA) / (1 / vT + 1 / vA); sC = core(comb[m])[1]
    rows.append((lo, hi, sT, sA, sC, np.sqrt(vT), np.sqrt(vA), np.sqrt(vC)))
    say('  A %3d-%3d n %6d: TOA %.2f  ADC %.2f  combined %.2f ns (%.0f mm)   [TOA-only %.2f, ADC-only %.2f, common %.2f]' % (lo, hi, m.sum(), sT, sA, sC, sC / pA, np.sqrt(vT), np.sqrt(vA), np.sqrt(vC)))
rr = np.array(rows); xm = 0.5 * (rr[:, 0] + rr[:, 1])
c = ROOT.TCanvas('c', '', 1600, 800); c.Divide(2, 1); keep = []
c.cd(1); fr = ROOT.gPad.DrawFrame(100, 0, 1000, rr[:, 2:5].max() * 1.6, 'SFIL2 cell 3, North - South;min(A_{S}, A_{N}) [ADC];core #sigma(t_{N} - t_{S}) [ns]'); fr.GetXaxis().SetNdivisions(505); keep.append(fr); lg = legend(0.36, 0.72, 0.86, 0.88)
for col_, mk, cl, lab in ((2, 20, ROOT.kRed + 1, 'TOA (baseline-corrected)'), (3, 21, ROOT.kBlue + 1, 'ADC 50 % crossing'), (4, 22, ROOT.kBlack, 'combined')):
    g = ROOT.TGraph(len(xm), xm.astype(float), rr[:, col_].astype(float)); g.SetMarkerStyle(mk); g.SetMarkerSize(1.6); g.SetMarkerColor(cl); g.SetLineColor(cl); g.Draw('pl same'); keep.append(g); lg.AddEntry(g, lab, 'pl')
lg.Draw(); keep.append(lg)
c.cd(2); fr = ROOT.gPad.DrawFrame(100, 0, 1000, rr[:, 5:8].max() * 1.6, 'three-cornered split;min(A_{S}, A_{N}) [ADC];contribution to #sigma(t_{N} - t_{S}) [ns]'); fr.GetXaxis().SetNdivisions(505); keep.append(fr); lg = legend(0.36, 0.72, 0.86, 0.88)
for col_, mk, cl, lab in ((5, 20, ROOT.kRed + 1, 'TOA only'), (6, 21, ROOT.kBlue + 1, 'ADC only'), (7, 22, ROOT.kBlack, 'common (light timing)')):
    g = ROOT.TGraph(len(xm), xm.astype(float), rr[:, col_].astype(float)); g.SetMarkerStyle(mk); g.SetMarkerSize(1.6); g.SetMarkerColor(cl); g.SetLineColor(cl); g.Draw('pl same'); keep.append(g); lg.AddEntry(g, lab, 'pl')
lg.Draw(); keep.append(lg)
c.SaveAs('toa_baseline.png'); c.SaveAs('toa_baseline.pdf'); say('wrote toa_baseline.png')
np.savez('toa_baseline.npz', pad=pad, even=even, **{k + e: T[e][k] for e in ('S', 'N') for k in ('t', 'ok', 'tc', 'A', 'base')})
say('wrote toa_baseline.npz')
