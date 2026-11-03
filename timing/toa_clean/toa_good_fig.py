"""Figure for the good-TOA selection: sigma_x per cut stage (TOA and ADC on the same showers) and N - S per paddle after all cuts.
    python3 toa_good_fig.py -> fig_good.png"""
import numpy as np, sys, ROOT
sys.path.insert(0, '/sphenix/user/tmengel/eic/bic_testbeam/common'); from plotstyle import setup, legend, headroom, finish
setup(); z = np.load('toa_clean.npz'); g = np.load('toa_good.npz'); R = g['rows']; keep = []
PB = {4: 2, 3: 1, 6: 8, 5: 4}; PX = {4: -80., 3: -30., 6: 30., 5: 80.}; PADS = (4, 3, 6, 5)
COL = {4: ROOT.kRed + 1, 3: ROOT.kOrange + 7, 6: ROOT.kGreen + 2, 5: ROOT.kBlue + 1}
c = ROOT.TCanvas('c', '', 1800, 800); c.Divide(2, 1)
c.cd(1); ROOT.gPad.SetBottomMargin(0.24); n = len(R)
h = ROOT.TH1F('fr', 'position resolution per shower;;#sigma_{x} [mm]', n, 0, n); h.SetMinimum(0); h.SetMaximum(380); h.SetStats(0)
for i, lab in enumerate(['all', '+quiet', '+1 word', '+pure']): h.GetXaxis().SetBinLabel(i + 1, '%s %.0f%%' % (lab, 100 * R[i, 0]))
h.GetXaxis().SetLabelSize(0.055); h.GetXaxis().LabelsOption('h'); ROOT.gPad.SetBottomMargin(0.16); h.Draw('axis'); keep.append(h)
x = np.arange(n) + 0.5
for lo, hi, mk, cl, lab in ((7, 8, 20, ROOT.kRed + 1, 'TOA'), (14, 15, 21, ROOT.kBlue + 1, 'ADC 50 % crossing, same showers')):
    mid = 0.5 * (R[:, lo] + R[:, hi]); err = 0.5 * (R[:, hi] - R[:, lo])
    gr = ROOT.TGraphErrors(n, x.astype(float), mid.astype(float), np.zeros(n), err.astype(float)); gr.SetMarkerStyle(mk); gr.SetMarkerSize(1.8); gr.SetMarkerColor(cl); gr.SetLineColor(cl); gr.SetLineWidth(3); gr.Draw('p same'); keep.append(gr)
lg = legend(0.20, 0.74, 0.90, 0.86, size=0.040); lg.AddEntry(keep[-2], 'TOA (range over paddles)', 'pe'); lg.AddEntry(keep[-1], 'ADC 50 % crossing, same showers', 'pe'); lg.Draw(); keep.append(lg)
c.cd(2); m = g['mgood'] & ~z['even']; d = g['dTgood']; hs = []
for p in PADS:
    v = d[m & (z['pad'] == PB[p])] - np.median(d[m]); hh = ROOT.TH1D('h%d' % p, '', 60, -12, 12); hh.FillN(len(v), np.ascontiguousarray(v), ROOT.nullptr); hh.Scale(1. / hh.Integral())
    hh.SetLineColor(COL[p]); hh.SetLineWidth(3); hs.append(hh); keep.append(hh)
fr = ROOT.gPad.DrawFrame(-12, 0, 12, 1, 'TOA North - South, good TOAs only;t_{N} - t_{S} - median [ns];fraction of showers'); headroom(fr, 0, max(x_.GetMaximum() for x_ in hs), 0.35); keep.append(fr)
lg = legend(0.20, 0.72, 0.85, 0.86, size=0.040, ncol=2)
for hh, p in zip(hs, PADS): hh.Draw('hist same'); lg.AddEntry(hh, 'P%d  x = %+.0f mm' % (p, PX[p]), 'l')
lg.Draw(); keep.append(lg)
finish(c); c.SaveAs('fig_good.png')
