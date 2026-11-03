"""Figures for the ADC-free SFIL2 cell 3 TOA study, from toa_clean.npz.
    python3 toa_clean_figs.py -> fig_codes.png, fig_tNtS.png, fig_ns_paddle.png, fig_xres.png"""
import numpy as np, sys, ROOT
sys.path.insert(0, '/sphenix/user/tmengel/eic/bic_testbeam/common'); from plotstyle import setup, legend, headroom, finish, colz
setup(); z = np.load('toa_clean.npz'); keep = []
PB = {4: 2, 3: 1, 6: 8, 5: 4}; PX = {4: -80., 3: -30., 6: 30., 5: 80.}; PADS = (4, 3, 6, 5)
COL = {4: ROOT.kRed + 1, 3: ROOT.kOrange + 7, 6: ROOT.kGreen + 2, 5: ROOT.kBlue + 1}
def th1(vals, nb, lo, hi, col, width=3, style=1, norm=True):
    h = ROOT.TH1D('h%d' % len(keep), '', nb, lo, hi); v = np.ascontiguousarray(vals, dtype=float); h.FillN(len(v), v, ROOT.nullptr)
    if norm: h.Scale(1. / max(h.Integral(), 1))
    h.SetLineColor(col); h.SetLineWidth(width); h.SetLineStyle(style); keep.append(h); return h
# 1. code density of the pulse word against the words seen without a pulse
c = ROOT.TCanvas('c1', '', 1800, 800); c.Divide(2, 1)
for i, (e, nm) in enumerate((('S', 'South'), ('N', 'North'))):
    c.cd(i + 1); c0 = int(z['c0' + e]); code = z['code' + e][z['has' + e]]; bg = z['bg' + e]
    hp = th1((code - c0) % 1024, 64, 0, 1024, ROOT.kBlack); hb = th1((bg - c0) % 1024, 64, 0, 1024, ROOT.kRed + 1, style=2)
    fr = ROOT.gPad.DrawFrame(0, 0, 1024, 1, 'SFIL2 cell 3 %s;unwrapped TOA code;fraction of words per 16 codes' % nm); headroom(fr, 0, max(hp.GetMaximum(), hb.GetMaximum()), 0.35); keep.append(fr)
    hp.Draw('hist same'); hb.Draw('hist same')
    lg = legend(0.20, 0.72, 0.85, 0.86, size=0.040); lg.AddEntry(hp, 'word on the pulse (sample k or k+1)', 'l'); lg.AddEntry(hb, 'words in triggers with no pulse', 'l'); lg.Draw(); keep.append(lg)
finish(c); c.SaveAs('fig_codes.png')
# 2. TOA time North against South, per shower
b = z['both'] & ~z['even']; tS = z['tS'][b]; tN = z['tN'][b]; mS = np.median(tS); mN = np.median(tN)
c = ROOT.TCanvas('c2', '', 1000, 850); colz(); ROOT.gPad.SetLogz()
h2 = ROOT.TH2D('h2', 'SFIL2 cell 3, TOA only;t_{S} - median [ns];t_{N} - median [ns]', 100, -30, 30, 100, -30, 30)
h2.FillN(len(tS), np.ascontiguousarray(tS - mS), np.ascontiguousarray(tN - mN), ROOT.nullptr); h2.Draw('colz'); keep.append(h2)
finish(c); c.SaveAs('fig_tNtS.png')
# 3. North - South per paddle: TOA (amplitude/baseline corrected) and ADC
c = ROOT.TCanvas('c3', '', 1800, 800); c.Divide(2, 1)
for i, (lab, key) in enumerate((('TOA', 'dT'), ('ADC 50 % crossing', 'dA'))):
    c.cd(i + 1); hs = []
    for p in PADS:
        v = z[key][b & (z['pad'] == PB[p])]; hs.append(th1(v - np.median(z[key][b]), 80, -15, 15, COL[p]))
    fr = ROOT.gPad.DrawFrame(-15, 0, 15, 1, 'SFIL2 cell 3, %s;t_{N} - t_{S} - median [ns];fraction of showers' % lab); headroom(fr, 0, max(h.GetMaximum() for h in hs), 0.35); keep.append(fr)
    lg = legend(0.20, 0.72, 0.85, 0.86, size=0.040, ncol=2)
    for h, p in zip(hs, PADS): h.Draw('hist same'); lg.AddEntry(h, 'P%d  x = %+.0f mm' % (p, PX[p]), 'l')
    lg.Draw(); keep.append(lg)
finish(c); c.SaveAs('fig_ns_paddle.png')
# 4. mean N - S against paddle position, and sigma_x per paddle
c = ROOT.TCanvas('c4', '', 1800, 800); c.Divide(2, 1); x = np.array([PX[p] for p in PADS], float)
SETS = (('TOA_raw', 'TOA, codes only', 24, ROOT.kGray + 2), ('TOA', 'TOA + amplitude/baseline', 20, ROOT.kRed + 1), ('ADC', 'ADC 50 % crossing', 21, ROOT.kBlue + 1))
c.cd(1); fr = ROOT.gPad.DrawFrame(-100, -2, 100, 2, 'mean North - South per paddle;paddle x along the fibres [mm];#LTt_{N} - t_{S}#GT - fit at x = 0 [ns]'); headroom(fr, -2, 1.5, 0.30); keep.append(fr); fr.GetXaxis().SetNdivisions(505)
lg = legend(0.20, 0.72, 0.85, 0.86, size=0.038)
for key, lab, mk, cl in SETS:
    g = ROOT.TGraph(4, x, (z[key + '_mu'] - z[key + '_c0']).astype(float)); g.SetMarkerStyle(mk); g.SetMarkerSize(1.8); g.SetMarkerColor(cl); g.SetLineColor(cl); g.Draw('pl same'); keep.append(g)
    lg.AddEntry(g, '%s: %.0f mm/ns' % (lab, 2 / float(z[key + '_c1'])), 'pl')
lg.Draw(); keep.append(lg)
c.cd(2); fr = ROOT.gPad.DrawFrame(-100, 0, 100, 400, 'position resolution per shower;paddle x along the fibres [mm];#sigma_{x} [mm]'); headroom(fr, 0, 300, 0.30); keep.append(fr); fr.GetXaxis().SetNdivisions(505)
lg = legend(0.20, 0.72, 0.85, 0.86, size=0.038)
for key, lab, mk, cl in SETS:
    g = ROOT.TGraph(4, x, (z[key + '_sig'] / float(z[key + '_c1'])).astype(float)); g.SetMarkerStyle(mk); g.SetMarkerSize(1.8); g.SetMarkerColor(cl); g.SetLineColor(cl); g.Draw('pl same'); keep.append(g)
    lg.AddEntry(g, lab, 'pl')
lg.Draw(); keep.append(lg)
finish(c); c.SaveAs('fig_xres.png')
