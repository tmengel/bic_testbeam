"""Stand-alone figure: share of the visible energy per element of the stack (SFIL1-5, Baby BCAL), 3-6 GeV positron
simulation and run 3062 (data_shares.npz: all paddles, and paddles P4 / P5 at the two ends of the module).
    python3 stack_shares.py  -> stack_shares.png
"""
import glob, sys, numpy as np, uproot, ROOT
sys.path.insert(0, '/sphenix/user/tmengel/eic/bic_testbeam/common'); from plotstyle import setup
share = {}
for f in sorted([f for f in glob.glob('/sphenix/tg/tg01/jets/tmengel/eic/sim/v1_purePb/bic_sim_e*_s*.root') if '_x' not in f]):
    a = uproot.open(f)['sim'].arrays(['ebeam', 'sfil_evis', 'bb_evis'], library='np'); E = a['ebeam'].mean()
    if not 2500 <= E <= 6500: continue
    l5 = a['sfil_evis'][:, :5, :].sum(axis=2); bb = a['bb_evis'].sum(axis=(1, 2)); tot = l5.sum(axis=1) + bb; ok = tot > 0
    share[E] = np.concatenate([l5[ok] / tot[ok][:, None], (bb[ok] / tot[ok])[:, None]], axis=1).mean(axis=0)
D = np.load('data_shares.npz')
setup(); ROOT.gStyle.SetPadLeftMargin(0.13); ROOT.gStyle.SetPadRightMargin(0.04); ROOT.gStyle.SetPadTopMargin(0.05); ROOT.gStyle.SetPadBottomMargin(0.13); ROOT.gROOT.ForceStyle()
c = ROOT.TCanvas('c', '', 1000, 750); keep = []
fr = ROOT.gPad.DrawFrame(0.5, 0, 6.5, 0.65, ';;share of the visible energy'); fr.GetXaxis().SetNdivisions(6); keep.append(fr)
for b, lab in enumerate(['SFIL1', 'SFIL2', 'SFIL3', 'SFIL4', 'SFIL5', 'Baby BCAL']): fr.GetXaxis().SetBinLabel(fr.GetXaxis().FindBin(b + 1), lab)
fr.GetXaxis().SetLabelSize(0.045); fr.GetXaxis().LabelsOption('h'); fr.GetYaxis().SetTitleOffset(1.3); fr.GetYaxis().SetTitleSize(0.045); fr.GetYaxis().SetLabelSize(0.04)
COL = {3000: 2, 4000: 4, 5000: 8, 6000: 6}; LS = {3000: 1, 4000: 2, 5000: 7, 6000: 9}; x6 = np.arange(1, 7, dtype=float)
lsim = ROOT.TLegend(0.16, 0.7, 0.4, 0.93); lsim.SetBorderSize(0); lsim.SetFillStyle(0); lsim.SetTextSize(0.034); lsim.SetHeader('simulation, e^{+}')
for e in sorted(share):
    g = ROOT.TGraph(6, x6, np.array(share[e])); g.SetLineColor(COL[e]); g.SetLineWidth(3); g.SetLineStyle(LS[e]); g.Draw('L'); keep.append(g); lsim.AddEntry(g, '%.0f GeV' % (e / 1000.), 'l')
ldat = ROOT.TLegend(0.47, 0.7, 0.94, 0.93); ldat.SetBorderSize(0); ldat.SetFillStyle(0); ldat.SetTextSize(0.034); ldat.SetHeader('run 3062, %d showers seen by all boards' % int(D['n']))
gd = ROOT.TGraph(6, x6, np.array(D['all'])); gd.SetMarkerStyle(20); gd.SetMarkerSize(2.0); gd.SetLineWidth(3); gd.Draw('LP'); keep.append(gd); ldat.AddEntry(gd, 'all paddles', 'lp')
g4 = ROOT.TGraph(6, x6 - 0.13, np.array(D['P4'])); g4.SetMarkerStyle(24); g4.SetMarkerSize(1.7); g4.Draw('P'); keep.append(g4); ldat.AddEntry(g4, 'paddle P4, North end', 'p')
g5 = ROOT.TGraph(6, x6 + 0.13, np.array(D['P5'])); g5.SetMarkerStyle(25); g5.SetMarkerSize(1.7); g5.Draw('P'); keep.append(g5); ldat.AddEntry(g5, 'paddle P5, South end', 'p')
lsim.Draw(); ldat.Draw(); keep += [lsim, ldat]
t = ROOT.TLatex(); t.SetNDC(); t.SetTextFont(42); t.SetTextSize(0.03); t.SetTextColor(ROOT.kGray + 2)
t.DrawLatex(0.16, 0.64, 'SFIL1-4: #sqrt{S#timesN} per cell, uncalibrated, beam row saturates'); t.DrawLatex(0.16, 0.60, 'SFIL5 (KCU1) and Baby BCAL (CAEN) scales set to the 4 GeV simulation,'); t.DrawLatex(0.16, 0.56, 'so their all-paddle points are a normalisation, not a measurement'); keep.append(t)
c.SaveAs('stack_shares.pdf'); print('wrote stack_shares.png')
