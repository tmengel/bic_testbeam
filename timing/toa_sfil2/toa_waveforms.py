"""SFIL2 cell 3 TOA after the full calibration (toa_baseline.npz): example waveforms with the sub-sample TOA time, pulse shape
aligned on the 25 ns tick vs on the TOA time, and the North - South horizontal position resolution in mm.
    python3 toa_waveforms.py -> toa_waveforms.txt, slides/wf_events.png, slides/wf_persist.png, slides/xres.png"""
import numpy as np, uproot, sys, ROOT
sys.path.insert(0, '/sphenix/user/tmengel/eic/bic_testbeam/common'); from plotstyle import setup, legend, colz
setup(); ROOT.gROOT.SetBatch(True); out = open('toa_waveforms.txt', 'w')
def say(*a): print(*a, flush=True); print(*a, file=out)
z = np.load('toa_baseline.npz'); NT = len(z['pad']); pad = z['pad']; odd = ~z['even']
PB = {4: 2, 3: 1, 6: 8, 5: 4}; PX = {4: -80., 3: -30., 6: 30., 5: 80.}; PCOL = {4: ROOT.kBlue + 1, 3: ROOT.kGreen + 2, 6: ROOT.kOrange + 7, 5: ROOT.kRed + 1}
RED, BLUE, BLK = ROOT.kRed + 1, ROOT.kBlue + 1, ROOT.kBlack
def core(v, w=2.5):
    v = v[np.isfinite(v)]
    if len(v) < 30: return np.nan, np.nan
    m = np.median(v); s = 1.4826 * np.median(np.abs(v - m))
    for _ in range(6): g = np.abs(v - m) < w * s; m, s = v[g].mean(), v[g].std()
    return m, s / 0.9
keep = []
f = uproot.open('/sphenix/user/tmengel/eic/bic_testbeam/data/kcu0_3062.root')
NW = 300000   # waveforms read for the displays
Y = {}
for e, br in (('S', 'Sc3'), ('N', 'Nc3')):
    W = f['sfil2'].arrays([br + '_adc'], entry_stop=NW, library='np')[br + '_adc'].astype(float); Y[e] = W - W[:, 1:4].mean(axis=1)[:, None]
both = z['okS'] & z['okN'] & odd

# 1. example events: one per paddle, typical amplitude, N - S near the paddle's mean
dT = z['tN'] - z['tS']; mu = {p: core(dT[both & (pad == PB[p])])[0] for p in PX}
c = ROOT.TCanvas('cev', '', 2400, 800); c.Divide(3, 1); ip = 0
for p in (4, 6, 5):
    cand = np.where(both[:NW] & (pad[:NW] == PB[p]) & (z['AS'][:NW] > 400) & (z['AS'][:NW] < 700) & (z['AN'][:NW] > 400) & (z['AN'][:NW] < 700) & (np.abs(dT[:NW] - mu[p]) < 0.3))[0]
    i = cand[0]; c.cd(ip + 1); ip += 1
    ymax = max(Y['S'][i].max(), Y['N'][i].max()); tlo = 25. * 3; thi = 25. * 14
    fr = ROOT.gPad.DrawFrame(tlo, -0.1 * ymax, thi, 2.1 * ymax, ';time [ns];ADC - pedestal'); fr.GetXaxis().SetNdivisions(506); keep.append(fr)
    lg = legend(0.2, 0.6, 0.92, 0.88, size=0.042)
    lg.SetHeader('P%d (x = %+.0f mm): t_{N} - t_{S} = %.2f ns' % (p, PX[p], dT[i]))
    for e, col, mk in (('S', BLUE, 20), ('N', RED, 21)):
        ts = 25. * np.arange(20); g = ROOT.TGraph(20, ts, np.ascontiguousarray(Y[e][i])); g.SetMarkerStyle(mk); g.SetMarkerSize(1.8); g.SetMarkerColor(col); g.SetLineColor(col); g.SetLineWidth(2); g.Draw('pl same'); keep.append(g)
        t = z['t' + e][i]; l = ROOT.TLine(t, -0.1 * ymax, t, 1.1 * ymax); l.SetLineColor(col); l.SetLineStyle(2); l.SetLineWidth(3); l.Draw(); keep.append(l)
        lg.AddEntry(g, '%s: t_{TOA} = %.1f ns' % ('South' if e == 'S' else 'North', t), 'pl')
    lg.Draw(); keep.append(lg)
    say('event %d (P%d): A_S %.0f A_N %.0f, t_S %.2f t_N %.2f, N - S %.2f ns (paddle mean %.2f)' % (i, p, z['AS'][i], z['AN'][i], z['tS'][i], z['tN'][i], dT[i], mu[p]))
c.SaveAs('slides/wf_events.png')

# 2. pulse shape: samples placed by the 25 ns tick (relative to the peak sample) vs by the TOA time
c = ROOT.TCanvas('cps', '', 2000, 1300); c.Divide(2, 2); ip = 0
for e, lab in (('S', 'South'), ('N', 'North')):
    sel = z['ok' + e][:NW] & odd[:NW] & (z['A' + e][:NW] > 300) & (z['A' + e][:NW] < 700); y = Y[e][sel] / z['A' + e][:NW][sel][:, None]
    ipk = Y[e][sel][:, 4:16].argmax(axis=1) + 4; ts = 25. * np.arange(20)[None, :]
    for how, xx in (('tick', ts - 25. * ipk[:, None]), ('TOA', ts - z['t' + e][:NW][sel][:, None])):
        ip += 1; c.cd(ip); colz(right=0.15); ROOT.gPad.SetLogz()
        h = ROOT.TH2F('h%s%s' % (e, how), ';%s;ADC / amplitude' % ('t - t_{peak sample} [ns]' if how == 'tick' else 't - t_{TOA} [ns]'), 200, -75, 125, 100, -0.2, 1.3); keep.append(h)
        h.FillN(xx.size, np.ascontiguousarray(xx.ravel()), np.ascontiguousarray(y.ravel()), np.ones(xx.size)); h.GetXaxis().SetNdivisions(505); h.Draw('colz')
        tx = ROOT.TLatex(); tx.SetNDC(); tx.SetTextSize(0.05); tx.DrawLatex(0.2, 0.8, '%s, aligned on the %s' % (lab, '25 ns tick' if how == 'tick' else 'TOA time')); keep.append(tx)
c.SaveAs('slides/wf_persist.png')

# 3. horizontal position from N - S: x = (dt - b) / a with each estimator's own paddle fit; resolution in mm
xs = np.array([PX[p] for p in (4, 3, 6, 5)]); X = {}; slope = {}
for lab, d in (('TOA', z['tN'] - z['tS']), ('ADC', z['tcN'] - z['tcS'])):
    m_ = [core(d[both & (pad == PB[p])])[0] for p in (4, 3, 6, 5)]; a, b = np.polyfit(xs, m_, 1); X[lab] = (d - b) / a; slope[lab] = a
    say('%s: slope %.4f ns/mm (v_eff %.0f mm/ns), offset %.2f ns' % (lab, a, 2 / a, b))
c = ROOT.TCanvas('cx', '', 2400, 800); c.Divide(3, 1)
for ip, lab in enumerate(('TOA', 'ADC')):
    c.cd(ip + 1); hs = []
    for p in (4, 3, 6, 5):
        h = ROOT.TH1F('hx%s%d' % (lab, p), '', 80, -600, 600); v = X[lab][both & (pad == PB[p])]; h.FillN(len(v), np.ascontiguousarray(v), np.ones(len(v))); h.Scale(1. / max(h.Integral(), 1)); hs.append((p, h, core(v)))
    ymax = max(h.GetMaximum() for _, h, _ in hs)
    fr = ROOT.gPad.DrawFrame(-600, 0, 600, ymax * 1.75, ';x from t_{N} - t_{S} (%s) [mm];fraction' % ('TOA' if lab == 'TOA' else 'ADC 50 % crossing')); fr.GetXaxis().SetNdivisions(505); keep.append(fr)
    lg = legend(0.2, 0.64, 0.92, 0.88, size=0.04)
    for p, h, (m_, s_) in hs:
        h.SetLineColor(PCOL[p]); h.SetLineWidth(2); h.Draw('hist same'); keep.append(h); lg.AddEntry(h, 'P%d (%+.0f mm): #sigma_{x} = %.0f mm' % (p, PX[p], s_), 'l')
        say('%s P%d: mean x %.0f mm, core sigma %.0f mm' % (lab, p, m_, s_))
    lg.Draw(); keep.append(lg)
# per bin of min(A): sigma_x for TOA, ADC, inverse-variance combination (three-cornered split in mm)
AE = [150, 250, 350, 450, 600, 950]; Amin = np.minimum(z['AS'], z['AN']); rows = []
rx = {}
for lab in ('TOA', 'ADC'):
    r = np.full(NT, np.nan)
    for p in (4, 3, 6, 5):
        mm = both & (pad == PB[p]); r[mm] = X[lab][mm] - core(X[lab][mm])[0]
    rx[lab] = r
say('\nper bin of min(A_S, A_N), odd triggers: core sigma_x [mm] TOA / ADC / combined')
for lo, hi in zip(AE[:-1], AE[1:]):
    m = both & (Amin >= lo) & (Amin < hi) & np.isfinite(rx['TOA']) & np.isfinite(rx['ADC'])
    sT = core(rx['TOA'][m])[1]; sA = core(rx['ADC'][m])[1]; sD = core((rx['TOA'] - rx['ADC'])[m])[1]
    vT = max((sT ** 2 + sD ** 2 - sA ** 2) / 2, 1.); vA = max((sA ** 2 + sD ** 2 - sT ** 2) / 2, 1.)
    sC = core(((rx['TOA'] / vT + rx['ADC'] / vA) / (1 / vT + 1 / vA))[m])[1]; rows.append((sT, sA, sC))
    say('  A %3d-%3d n %6d: TOA %.0f  ADC %.0f  combined %.0f mm' % (lo, hi, m.sum(), sT, sA, sC))
rr = np.array(rows); xm = np.array([200, 300, 400, 525, 775.])
c.cd(3); fr = ROOT.gPad.DrawFrame(100, 0, 1000, 450, ';min(A_{S}, A_{N}) [ADC];core #sigma_{x} [mm]'); fr.GetXaxis().SetNdivisions(505); keep.append(fr); lg = legend(0.25, 0.7, 0.92, 0.88, size=0.042)
for j, mk, cl, lab in ((0, 20, RED, 'TOA'), (1, 21, BLUE, 'ADC 50 % crossing'), (2, 22, BLK, 'TOA + ADC combined')):
    g = ROOT.TGraph(5, xm, np.ascontiguousarray(rr[:, j])); g.SetMarkerStyle(mk); g.SetMarkerSize(1.8); g.SetMarkerColor(cl); g.SetLineColor(cl); g.SetLineWidth(2); g.Draw('pl same'); keep.append(g); lg.AddEntry(g, lab, 'pl')
lg.Draw(); keep.append(lg)
c.SaveAs('slides/xres.png'); say('wrote slides/wf_events.png, wf_persist.png, xres.png')
