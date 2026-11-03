"""SFIL2 cell 3: x from North - South (calibrated TOA, and the ADC 50 % crossing) against the AstroPix x of the matched cluster.
x from N - S uses each estimator's own paddle fit (as in toa_waveforms.py); AstroPix in the stack frame x = astro_x - 88.75 mm.
    python3 toa_astropix.py [min_purity] -> toa_astropix.txt, slides/astropix.png"""
import numpy as np, uproot, sys, ROOT
sys.path.insert(0, '/sphenix/user/tmengel/eic/bic_testbeam/common'); from plotstyle import setup, legend, colz
setup(); ROOT.gROOT.SetBatch(True); out = open('toa_astropix.txt', 'w')
def say(*a): print(*a, flush=True); print(*a, file=out)
PUR = float(sys.argv[1]) if len(sys.argv) > 1 else 0.3
z = np.load('toa_baseline.npz'); NT = len(z['pad']); pad = z['pad']; odd = ~z['even']
PB = {4: 2, 3: 1, 6: 8, 5: 4}; PX = {4: -80., 3: -30., 6: 30., 5: 80.}; RED, BLUE, BLK = ROOT.kRed + 1, ROOT.kBlue + 1, ROOT.kBlack
def core(v, w=2.5):
    v = v[np.isfinite(v)]
    if len(v) < 30: return np.nan, np.nan
    m = np.median(v); s = 1.4826 * np.median(np.abs(v - m))
    for _ in range(6): g = np.abs(v - m) < w * s; m, s = v[g].mean(), v[g].std()
    return m, s / 0.9
ts0 = uproot.open('/sphenix/user/tmengel/eic/bic_testbeam/data/kcu0_3062.root')['header']['ts0'].array(library='np', entry_stop=NT)
ev = uproot.open('/sphenix/user/tmengel/eic/bic_testbeam/data/events_3062.root')['events'].arrays(['trig_ts', 'astro_x', 'astro_purity'], library='np')
ets = ev['trig_ts']; order = np.argsort(ets); ets_s = ets[order]; pos = np.clip(np.searchsorted(ets_s, ts0), 0, len(ets_s) - 1); hit = ets_s[pos] == ts0; row = np.where(hit, order[pos], -1)
ax = np.full(NT, np.nan); ax[hit] = ev['astro_x'][row[hit]] - 88.75; pur = np.zeros(NT); pur[hit] = ev['astro_purity'][row[hit]]
both = z['okS'] & z['okN'] & odd
xs = np.array([PX[p] for p in (4, 3, 6, 5)]); X = {}
for lab, d in (('TOA', z['tN'] - z['tS']), ('ADC', z['tcN'] - z['tcS'])):
    m_ = [core(d[both & (pad == PB[p])])[0] for p in (4, 3, 6, 5)]; a, b = np.polyfit(xs, m_, 1); X[lab] = (d - b) / a
sel = both & np.isfinite(ax) & (pur >= PUR)
say('odd triggers, both TOA ends good: %d; with an AstroPix cluster of purity >= %.1f: %d; AstroPix x 5-95 %%: %.0f to %.0f mm' % (both.sum(), PUR, sel.sum(), *np.percentile(ax[sel], [5, 95])))
RES = {}
for lab in ('TOA', 'ADC'):
    r = X[lab][sel] - ax[sel]; m_, s_ = core(r); g = np.abs(r - m_) < 600
    sl = np.polyfit(ax[sel][g], X[lab][sel][g], 1)[0]; cc = np.corrcoef(ax[sel][g], X[lab][sel][g])[0, 1]
    # within a paddle: does x follow AstroPix once the paddle mean is removed?
    rin = []
    for p in (4, 3, 6, 5):
        mm = sel & (pad == PB[p])
        if mm.sum() > 100:
            dx = ax[mm] - np.median(ax[mm]); dX = X[lab][mm] - np.median(X[lab][mm]); gg = np.abs(dX) < 600; rin.append((p, mm.sum(), np.corrcoef(dx[gg], dX[gg])[0, 1], ax[mm].std()))
    RES[lab] = (m_, s_)
    say('%s: x - x_AstroPix core mean %.0f mm, sigma %.0f mm; slope of x vs x_AstroPix %.2f, corr %.3f; within paddle corr %s' % (lab, m_, s_, sl, cc, ', '.join('P%d (n %d) corr %.3f, AstroPix rms %.0f mm' % q for q in rin)))
keep = []
c = ROOT.TCanvas('c', '', 2400, 800); c.Divide(3, 1)
c.cd(1); colz(right=0.15); h2 = ROOT.TH2F('h2', ';x_{AstroPix} [mm];x from TOA t_{N} - t_{S} [mm]', 40, -200, 200, 60, -600, 600); keep.append(h2)
h2.FillN(int(sel.sum()), np.ascontiguousarray(ax[sel]), np.ascontiguousarray(X['TOA'][sel]), np.ones(int(sel.sum()))); h2.GetXaxis().SetNdivisions(505); h2.Draw('colz')
l = ROOT.TLine(-200, -200, 200, 200); l.SetLineStyle(2); l.SetLineWidth(3); l.Draw(); keep.append(l)
c.cd(2); fr = ROOT.gPad.DrawFrame(-200, -200, 200, 330, ';median x_{AstroPix} [mm];mean x from t_{N} - t_{S} [mm]'); fr.GetXaxis().SetNdivisions(505); keep.append(fr); lg = legend(0.2, 0.7, 0.9, 0.88)
for lab, mk, cl in (('TOA', 20, RED), ('ADC', 21, BLUE)):
    px, py, ex, ey = [], [], [], []
    for p in (4, 3, 6, 5):
        mm = sel & (pad == PB[p])
        if mm.sum() < 50: continue
        m_, s_ = core(X[lab][mm]); px.append(np.median(ax[mm])); ex.append(ax[mm].std() / np.sqrt(mm.sum())); py.append(m_); ey.append(s_ / np.sqrt(mm.sum()))
        say('  %s P%d: n %d, AstroPix x median %.0f mm, x from N - S %.0f +- %.0f mm' % (lab, p, mm.sum(), px[-1], m_, ey[-1]))
    g = ROOT.TGraphErrors(len(px), np.array(px), np.array(py), np.array(ex), np.array(ey)); g.SetMarkerStyle(mk); g.SetMarkerSize(2); g.SetMarkerColor(cl); g.SetLineColor(cl); g.SetLineWidth(2); g.Draw('p same'); keep.append(g)
    lg.AddEntry(g, ('TOA' if lab == 'TOA' else 'ADC 50 % crossing') + ' (per paddle)', 'pe')
l2 = ROOT.TLine(-200, -200, 200, 200); l2.SetLineStyle(2); l2.SetLineWidth(2); l2.Draw(); keep.append(l2); lg.AddEntry(l2, 'x = x_{AstroPix}', 'l'); lg.Draw(); keep.append(lg)
c.cd(3); hs = []
for lab, cl in (('TOA', RED), ('ADC', BLUE)):
    h = ROOT.TH1F('hr' + lab, '', 80, -800, 800); r = X[lab][sel] - ax[sel]; h.FillN(len(r), np.ascontiguousarray(r), np.ones(len(r))); h.Scale(1. / max(h.Integral(), 1)); h.SetLineColor(cl); h.SetLineWidth(3); hs.append((lab, h)); keep.append(h)
fr = ROOT.gPad.DrawFrame(-800, 0, 800, max(h.GetMaximum() for _, h in hs) * 1.6, ';x - x_{AstroPix} [mm];fraction'); fr.GetXaxis().SetNdivisions(505); keep.append(fr); lg = legend(0.2, 0.74, 0.9, 0.88, size=0.042)
for lab, h in hs: h.Draw('hist same'); lg.AddEntry(h, '%s: #sigma = %.0f mm' % ('TOA' if lab == 'TOA' else 'ADC crossing', RES[lab][1]), 'l')
lg.Draw(); keep.append(lg)
c.SaveAs('slides/astropix.png'); say('wrote slides/astropix.png')
