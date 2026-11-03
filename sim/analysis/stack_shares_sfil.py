"""Presentation figure: share of the visible energy in the four KCU0 SFIL layers (SFIL1-4 only, renormalised to
their sum), 3-6 GeV e+ simulation (original pure-Pb geometry, sim/v1_purePb) against run 3062 (all paddles,
data_shares.npz).  Baby BCAL, SFIL5 and the per-paddle series are left out.
    python3 stack_shares_sfil.py [simdir] [tag]   -> stack_shares_sfil<tag>.pdf (+ .png)
"""
import glob, sys, numpy as np, uproot, ROOT
sys.path.insert(0, '/sphenix/user/tmengel/eic/bic_testbeam/common'); from plotstyle import setup, legend, label
SIMDIR = sys.argv[1] if len(sys.argv) > 1 else '/sphenix/tg/tg01/jets/tmengel/eic/sim/v1_purePb'
TAG = sys.argv[2] if len(sys.argv) > 2 else ''
NL = 4
share = {}
for f in sorted(f for f in glob.glob(SIMDIR + '/bic_sim_e*_s*.root') if '_x' not in f):
    a = uproot.open(f)['sim'].arrays(['ebeam', 'sfil_evis'], library='np'); E = a['ebeam'].mean()
    if not 2500 <= E <= 6500: continue
    l = a['sfil_evis'][:, :NL, :].sum(axis=2); tot = l.sum(axis=1); ok = tot > 0
    share[E] = (l[ok] / tot[ok][:, None]).mean(axis=0)
D = np.load('data_shares.npz'); dat = np.array(D['all'][:NL]); dat /= dat.sum()
for e in sorted(share): print('sim %.0f GeV' % (e / 1e3), share[e].round(3))
print('data', dat.round(3), 'n =', int(D['n']))

setup()
ROOT.gStyle.SetPadLeftMargin(0.14); ROOT.gStyle.SetPadRightMargin(0.04); ROOT.gStyle.SetPadTopMargin(0.08); ROOT.gStyle.SetPadBottomMargin(0.14); ROOT.gROOT.ForceStyle()
c = ROOT.TCanvas('c', '', 1000, 750); keep = []
YMAX = 0.75
fr = ROOT.TH1F('fr', '', NL, 0.5, NL + 0.5); fr.SetMinimum(0); fr.SetMaximum(YMAX); fr.SetStats(0); fr.Draw('AXIS'); keep.append(fr)
ax, ay = fr.GetXaxis(), fr.GetYaxis()
for b in range(NL): ax.SetBinLabel(b + 1, 'SFIL%d' % (b + 1))          # bin labels sit at the bin centres x = 1..4
ax.SetLabelSize(0.06); ax.SetLabelOffset(0.015); ax.SetTickLength(0)   # ROOT would put the ticks at the bin EDGES: draw our own below
ay.SetTitle('share of the visible energy in SFIL1-4'); ay.SetTitleSize(0.05); ay.SetTitleOffset(1.25); ay.SetLabelSize(0.045); ay.SetNdivisions(505)
tick = 0.03 * YMAX                                                     # ticks at the layer positions, bottom and top, aligned with the labels
for x in range(1, NL + 1):
    for y0, y1 in ((0, tick), (YMAX, YMAX - tick)):
        ln = ROOT.TLine(x, y0, x, y1); ln.SetLineWidth(2); ln.Draw(); keep.append(ln)

COL = {3000: ROOT.kRed + 1, 4000: ROOT.kAzure + 2, 5000: ROOT.kGreen + 2, 6000: ROOT.kViolet - 1}; LS = {3000: 1, 4000: 2, 5000: 7, 6000: 9}
xs = np.arange(1, NL + 1, dtype=float)
lsim = legend(0.17, 0.66, 0.42, 0.89, size=0.04); lsim.SetHeader('simulation, e^{+}')
for e in sorted(share):
    g = ROOT.TGraph(NL, xs, np.array(share[e])); g.SetLineColor(COL[e]); g.SetLineWidth(4); g.SetLineStyle(LS[e]); g.Draw('L'); keep.append(g)
    lsim.AddEntry(g, '%.0f GeV' % (e / 1000.), 'l')
ldat = legend(0.4, 0.78, 0.85, 0.90, size=0.04); ldat.SetHeader('run 3062, e^{+}, %dk showers' % round(int(D['n']) / 1000.))
gd = ROOT.TGraph(NL, xs, dat); gd.SetMarkerStyle(20); gd.SetMarkerSize(2.2); gd.SetLineWidth(3); gd.Draw('LP'); keep.append(gd)
ldat.AddEntry(gd, '#sqrt{S x N} ADC amplitude (all paddles)', 'lp')
lsim.Draw(); ldat.Draw(); keep += [lsim, ldat]
keep.append(label('BIC test beam, JLab Hall D, August 2026', size=0.045))
t = ROOT.TLatex(); t.SetNDC(); t.SetTextFont(42); t.SetTextSize(0.032); t.SetTextColor(ROOT.kGray + 2); t.SetTextAlign(31)
t.DrawLatex(0.9, 0.72, 'data uncalibrated: ADC shares, not energy shares'); keep.append(t)
ROOT.gPad.RedrawAxis()
c.SaveAs('stack_shares_sfil%s.pdf' % TAG); c.SaveAs('stack_shares_sfil%s.png' % TAG); print('wrote stack_shares_sfil%s.pdf' % TAG)
