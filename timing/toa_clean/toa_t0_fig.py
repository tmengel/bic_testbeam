"""Why the TOA North-vs-South scatter is correlated, and where the North-South anticorrelation is (toa_t0.py).
    python3 toa_t0_fig.py -> fig_t0.png
(a) t_N vs t_S as before; (b) (t_N + t_S)/2 against the CAEN event time t0; (c) t_N vs t_S with t0 removed;
(d) per paddle, median t_N - t0 and t_S - t0 against the paddle position."""
import numpy as np, sys, ROOT
sys.path.insert(0, '/sphenix/user/tmengel/eic/bic_testbeam/common'); from plotstyle import setup, legend, headroom, finish, colz
setup(); keep = []
z = np.load('toa_clean.npz'); t0 = np.load('toa_t0.npz')['t0']
PB = {4: 2, 3: 1, 6: 8, 5: 4}; PX = {4: -80., 3: -30., 6: 30., 5: 80.}
b = z['both'] & ~z['even'] & np.isfinite(t0); tN = z['tN'][b]; tS = z['tS'][b]; tt = t0[b]; pad = z['pad'][b]
mN = np.median(tN); mS = np.median(tS); mNc = np.median(tN - tt); mSc = np.median(tS - tt); ms = np.median((tN + tS) / 2)
def h2(name, title, x, y, lo=-30, hi=30, yhi=30):
    h = ROOT.TH2D(name, title, 100, lo, hi, 100, lo, yhi); h.FillN(len(x), np.ascontiguousarray(x, dtype=float), np.ascontiguousarray(y, dtype=float), ROOT.nullptr)
    h.SetMinimum(0.5); h.Draw('colz'); keep.append(h); return h
c = ROOT.TCanvas('c', '', 1800, 1600); c.Divide(2, 2)
c.cd(1); colz(); ROOT.gPad.SetLogz()
h2('ha', 'SFIL2 cell 3, TOA, %d showers;t_{S} - median [ns];t_{N} - median [ns]' % len(tN), tS - mS, tN - mN)
c.cd(2); colz(); ROOT.gPad.SetLogz()
h2('hb', 'common time vs CAEN event time;t_{0} (CAEN) - median [ns];(t_{N} + t_{S})/2 - median [ns]', tt - np.median(tt), (tN + tS) / 2 - ms, yhi=45)
ln = ROOT.TLine(-25, -25, 25, 25); ln.SetLineColor(ROOT.kRed + 1); ln.SetLineWidth(3); ln.SetLineStyle(2); ln.Draw(); keep.append(ln)
lg = legend(0.20, 0.80, 0.70, 0.86, size=0.040); lg.AddEntry(ln, 'slope 1 (fit: %.2f)' % np.polyfit(tt, (tN + tS) / 2, 1)[0], 'l'); lg.Draw(); keep.append(lg)
c.cd(3); colz(); ROOT.gPad.SetLogz()
h2('hc', 'event time t_{0} removed;t_{S} - t_{0} - median [ns];t_{N} - t_{0} - median [ns]', tS - tt - mSc, tN - tt - mNc)
c.cd(4)
gs = {}
for e, v, col, mk in (('N', tN - tt - mNc, ROOT.kBlue + 1, 20), ('S', tS - tt - mSc, ROOT.kRed + 1, 21), ('C', (tN + tS) / 2 - tt - np.median((tN + tS) / 2 - tt), ROOT.kGray + 2, 24)):
    g = ROOT.TGraphErrors()
    for p in (4, 3, 6, 5):
        x = v[pad == PB[p]]; w = (np.percentile(x, 84) - np.percentile(x, 16)) / 2
        g.SetPoint(g.GetN(), PX[p], np.median(x)); g.SetPointError(g.GetN() - 1, 0, 1.2533 * w / np.sqrt(len(x)))
    g.SetMarkerColor(col); g.SetLineColor(col); g.SetMarkerStyle(mk); g.SetMarkerSize(1.8); g.SetLineWidth(2); gs[e] = g; keep.append(g)
fr = ROOT.gPad.DrawFrame(-110, -3, 110, 4, 'per paddle, t_{0} removed;paddle position x [mm];median time - median [ns]'); keep.append(fr)
fr.GetXaxis().SetNdivisions(505)
gs['C'].SetLineStyle(2)
for e in ('C', 'N', 'S'): gs[e].Draw('PL same')
lg = legend(0.20, 0.74, 0.90, 0.86, size=0.040, ncol=2)
lg.AddEntry(gs['N'], 't_{N} - t_{0}', 'pl'); lg.AddEntry(gs['S'], 't_{S} - t_{0}', 'pl'); lg.AddEntry(gs['C'], '(t_{N} + t_{S})/2 - t_{0}', 'pl'); lg.Draw(); keep.append(lg)
finish(c); c.SaveAs('fig_t0.png')
