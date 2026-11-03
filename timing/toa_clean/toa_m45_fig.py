"""The -45 deg view of the SFIL2 cell 3 TOA: t_N against t_S per paddle, with the event time t0 removed (toa_t0.py).
A shift of the shower along the fibre moves (t_S, t_N) along -45 deg; the paddles are the only position handle that is
larger than the per-end resolution, so the plot shows the per-paddle distributions and their centroids.
    python3 toa_m45_fig.py -> fig_m45.png, toa_m45.txt"""
import numpy as np, sys, ROOT
sys.path.insert(0, '/sphenix/user/tmengel/eic/bic_testbeam/common'); from plotstyle import setup, legend, finish
setup(); keep = []
out = open('toa_m45.txt', 'w')
def say(*a): print(*a, flush=True); print(*a, file=out)
z = np.load('toa_clean.npz'); t0 = np.load('toa_t0.npz')['t0']
PB = {4: 2, 3: 1, 6: 8, 5: 4}; PX = {4: -80., 3: -30., 6: 30., 5: 80.}; PADS = (4, 3, 6, 5)
COL = {4: ROOT.kRed + 1, 3: ROOT.kOrange + 7, 6: ROOT.kGreen + 2, 5: ROOT.kBlue + 1}; MK = {4: 20, 3: 21, 6: 22, 5: 33}
b = z['both'] & ~z['even'] & np.isfinite(t0); pad = z['pad'][b]
xN = z['tN'][b] - t0[b]; xS = z['tS'][b] - t0[b]; xN -= np.median(xN); xS -= np.median(xS)
C = {}
for p in PADS:
    m = pad == PB[p]; mN = np.median(xN[m]); mS = np.median(xS[m])
    eN = 1.2533 * (np.percentile(xN[m], 84) - np.percentile(xN[m], 16)) / 2 / np.sqrt(m.sum()); eS = 1.2533 * (np.percentile(xS[m], 84) - np.percentile(xS[m], 16)) / 2 / np.sqrt(m.sum())
    C[p] = (mS, mN, eS, eN, m.sum())
    say('P%d x %+4.0f mm: %5d showers; centroid t_S %+.2f +- %.2f, t_N %+.2f +- %.2f ns; along -45 deg (t_N - t_S)/sqrt2 %+.2f, along +45 deg (t_N + t_S)/sqrt2 %+.2f'
        % (p, PX[p], m.sum(), mS, eS, mN, eN, (mN - mS) / np.sqrt(2), (mN + mS) / np.sqrt(2)))
d = np.array(C[5][:2]) - np.array(C[4][:2])
say('P4 -> P5 (160 mm): dt_S %+.2f, dt_N %+.2f ns; -45 deg component %.2f ns, +45 deg component %.2f ns' % (d[0], d[1], (d[1] - d[0]) / np.sqrt(2), (d[1] + d[0]) / np.sqrt(2)))
c = ROOT.TCanvas('c', '', 1800, 850); c.Divide(2, 1)
# (a) per-paddle 1 sigma contours (39 % of the showers inside)
c.cd(1); R = 12
fr = ROOT.gPad.DrawFrame(-R, -R, R, R, 'per-paddle contours, 39 % of showers inside;t_{S} - t_{0} - median [ns];t_{N} - t_{0} - median [ns]'); keep.append(fr)
ln = ROOT.TLine(-6, 6, 3.5, -3.5); ln.SetLineStyle(2); ln.SetLineWidth(2); ln.SetLineColor(ROOT.kGray + 2); ln.Draw(); keep.append(ln)
lg = legend(0.56, 0.17, 0.93, 0.40, size=0.038)
for p in PADS:
    m = pad == PB[p]; h = ROOT.TH2D('h%d' % p, '', 48, -R, R, 48, -R, R)
    h.FillN(int(m.sum()), np.ascontiguousarray(xS[m], dtype=float), np.ascontiguousarray(xN[m], dtype=float), ROOT.nullptr); h.Smooth(1)
    v = np.sort(np.array([h.GetBinContent(i, j) for i in range(1, 49) for j in range(1, 49)]))[::-1]
    lev = v[np.searchsorted(np.cumsum(v) / v.sum(), 0.39)]
    h.SetContour(1, np.array([lev], dtype=float)); h.SetLineColor(COL[p]); h.SetLineWidth(3); h.Draw('cont3 same'); keep.append(h)
    g = ROOT.TGraph(1, np.array([C[p][0]]), np.array([C[p][1]])); g.SetMarkerStyle(MK[p]); g.SetMarkerColor(COL[p]); g.SetLineColor(COL[p]); g.SetLineWidth(3); g.SetMarkerSize(2.0); g.Draw('P same'); keep.append(g)
    lg.AddEntry(g, 'P%d, x = %+.0f mm' % (p, PX[p]), 'pl')
lg.AddEntry(ln, '-45#circ: position', 'l'); lg.Draw(); keep.append(lg)
# (b) the centroids, zoomed
c.cd(2)
fr = ROOT.gPad.DrawFrame(-3, -3, 3, 3, 'per-paddle centroids (medians);t_{S} - t_{0} - median [ns];t_{N} - t_{0} - median [ns]'); keep.append(fr)
lg2 = legend(0.58, 0.17, 0.93, 0.29, size=0.036)
for a, st, lab in ((-1, 2, '-45#circ: position'), (1, 3, '+45#circ: common time')):
    l = ROOT.TLine(-1.2, -1.2 * a, 1.2, 1.2 * a); l.SetLineStyle(st); l.SetLineWidth(2); l.SetLineColor(ROOT.kGray + 1); l.Draw(); keep.append(l); lg2.AddEntry(l, lab, 'l')
lg2.Draw(); keep.append(lg2)
for p in PADS:
    mS, mN, eS, eN, n = C[p]
    g = ROOT.TGraphErrors(1, np.array([mS]), np.array([mN]), np.array([eS]), np.array([eN])); g.SetMarkerStyle(MK[p]); g.SetMarkerColor(COL[p]); g.SetLineColor(COL[p]); g.SetMarkerSize(2.4); g.SetLineWidth(2); g.Draw('P same'); keep.append(g)
    t = ROOT.TLatex(); t.SetTextFont(42); t.SetTextSize(0.040); t.SetTextColor(COL[p]); t.SetTextAlign(12); t.DrawLatex(mS + 0.22, mN, 'P%d (%+.0f mm)' % (p, PX[p])); keep.append(t)
finish(c); c.SaveAs('fig_m45.png')
