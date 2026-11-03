"""The TOA North-vs-South scatter of fig_tNtS.png without the median subtraction: times from the start of the KCU0 readout
window, t = 25 s + T(code). Dashed lines mark the 25 ns sample boundaries.
    python3 toa_tNtS_abs.py -> fig_tNtS_abs.png"""
import numpy as np, sys, ROOT
sys.path.insert(0, '/sphenix/user/tmengel/eic/bic_testbeam/common'); from plotstyle import setup, finish, colz
setup(); keep = []
z = np.load('toa_clean.npz'); b = z['both'] & ~z['even']; tS = z['tS'][b]; tN = z['tN'][b]
lo, hi = 125, 225
c = ROOT.TCanvas('c', '', 1000, 850); colz(); ROOT.gPad.SetLogz()
h = ROOT.TH2D('h', 'SFIL2 cell 3, TOA only, %d showers;t_{S} [ns];t_{N} [ns]' % len(tS), 100, lo, hi, 100, lo, hi)
h.FillN(len(tS), np.ascontiguousarray(tS), np.ascontiguousarray(tN), ROOT.nullptr); h.SetMinimum(0.5)
for ax in (h.GetXaxis(), h.GetYaxis()): ax.SetNdivisions(4, 5, 0, False)            # labels on the sample boundaries
h.Draw('colz')
for v in range(150, hi, 25):
    for l in (ROOT.TLine(v, lo, v, hi), ROOT.TLine(lo, v, hi, v)):
        l.SetLineStyle(2); l.SetLineColor(ROOT.kGray + 1); l.SetLineWidth(1); l.Draw(); keep.append(l)
print('showers %d; outside the range: %d; median t_S %.2f, t_N %.2f' % (len(tS), ((tS < lo) | (tS > hi) | (tN < lo) | (tN > hi)).sum(), np.median(tS), np.median(tN)))
finish(c); c.SaveAs('fig_tNtS_abs.png')
