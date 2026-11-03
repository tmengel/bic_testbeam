"""Trigger rates of run 3062 versus time from the stream clocks (common/cache/index.npz made by common/bic_streams.py).
    python3 beam_rate.py  -> ../docs/beam_rate.png and numbers on stdout
"""
import numpy as np, ROOT, sys, os
sys.path.insert(0, '/sphenix/user/tmengel/eic/bic_testbeam/common'); from plotstyle import setup, label
d = np.load('../common/cache/index.npz')
k0 = d['k0_ts'] * 6.25e-9; k1 = (d['k1_ts'] + 1143665946) * 6.25e-9; ca = d['c_gtt'] * 8.533e-9 + 12.302395
T = np.arange(0, 1490, 10.)
def rate(t): h, _ = np.histogram(t, bins=T); return h / 10.
r0, r1, rc = rate(k0), rate(k1), rate(ca)
print('KCU0 recorded triggers: %d over %.0f s = %.2f kHz mean; 10 s bins: median %.2f kHz, p10 %.2f, p90 %.2f, bins at 0: %d of %d' % (len(k0), k0[-1], len(k0) / k0[-1] / 1e3, np.median(r0) / 1e3, np.percentile(r0, 10) / 1e3, np.percentile(r0, 90) / 1e3, (r0 == 0).sum(), len(r0)))
print('KCU1: %d, mean %.2f kHz, bins at 0: %d' % (len(k1), len(k1) / 1479 / 1e3, (r1 == 0).sum()))
print('CAEN: %d, mean %.2f kHz, median per 10 s %.2f kHz, bins at 0: %d' % (len(ca), len(ca) / 1471 / 1e3, np.median(rc) / 1e3, (rc == 0).sum()))
dt = np.diff(k0) * 1e6; print('KCU0 spacing between recorded triggers: min %.1f us, p1 %.1f, median %.1f us' % (dt.min(), np.percentile(dt, 1), np.median(dt)))
setup(); c = ROOT.TCanvas('c', '', 1400, 500); keep = []
fr = ROOT.gPad.DrawFrame(0, 0, 1490, 5, ';time in run [s];recorded triggers [kHz]'); keep.append(fr)
x = 0.5 * (T[1:] + T[:-1])
for r, col, nm in ((r0, 1, 'KCU0'), (rc, 2, 'CAEN'), (r1, 4, 'KCU1')):
    g = ROOT.TGraph(len(x), x, r / 1e3); g.SetLineColor(col); g.SetLineWidth(2); g.Draw('L'); keep.append(g)
lg = ROOT.TLegend(0.6, 0.6, 0.92, 0.84); lg.SetBorderSize(0); lg.SetFillStyle(0)
for g, nm in zip(keep[1:], ['KCU0 recorded (4.10 kHz issued)', 'CAEN recorded', 'KCU1 recorded']): lg.AddEntry(g, nm, 'l')
lg.Draw(); keep.append(lg); label('run 3062, 10 s bins', x=0.17, y=0.93, size=0.04)
c.SaveAs('../docs/beam_rate.png'); print('wrote docs/beam_rate.png')
