"""Stand-alone figure, paddle P5 only: South-end share of the visible energy per element (SFIL1-5, Baby BCAL).
Simulation: beam under P5 (x = +80 mm, 70 mm from the South fibre end), South-end light (sfil_ls, bb_ls), 5 and 6 GeV e+
(P5 is the high-energy end of the module: 0.68 GeV in May at 0.19 T, about 5.9 GeV at 1.662 T if the geometry held).
Data: run 3062 P5 triple-matched showers, South-end amplitudes only (data_p5_south.npz); the KCU1 (SFIL5) and CAEN
(Baby BCAL) scales are set so that their means relative to SFIL1-4 equal the 4 GeV simulation.
    python3 stack_shares_p5.py [simdir]  -> stack_shares_p5.png   (data_p5_south.npz regenerated 2026-10-23 with the correct paddle bitmask: 33468 P5 showers)
"""
import glob, sys, numpy as np, uproot, ROOT
sys.path.insert(0, '/sphenix/user/tmengel/eic/bic_testbeam/common'); from plotstyle import setup
share, means = {}, {}
SIMDIR = sys.argv[1] if len(sys.argv) > 1 else '/sphenix/tg/tg01/jets/tmengel/eic/sim/v1_purePb'
for f in sorted(glob.glob(SIMDIR + '/bic_sim_e*_x80_s*.root')):
    a = uproot.open(f)['sim'].arrays(['ebeam', 'sfil_ls', 'bb_ls'], library='np'); E = a['ebeam'].mean()
    if E < 4500: continue   # show the 5 and 6 GeV curves only
    l5 = a['sfil_ls'][:, :5, :].sum(axis=2); bb = a['bb_ls'].sum(axis=(1, 2)); tot = l5.sum(axis=1) + bb; ok = tot > 0
    share[E] = np.concatenate([l5[ok] / tot[ok][:, None], (bb[ok] / tot[ok])[:, None]], axis=1).mean(axis=0); means[E] = np.concatenate([l5.mean(axis=0), [bb.mean()]])
if not share: sys.exit('no x80 simulation files yet')
D = np.load('data_p5_south.npz'); raw = D['raw']; ref = np.mean([means[e] for e in share], axis=0)   # normalisation: mean of the shown energies
m = raw.mean(axis=0); k5 = (ref[4] / ref[:4].sum()) * m[:4].sum() / m[4]; kbb = (ref[5] / ref[:4].sum()) * m[:4].sum() / m[5]
cal = raw.copy(); cal[:, 4] *= k5; cal[:, 5] *= kbb; dshare = (cal / cal.sum(axis=1)[:, None]).mean(axis=0); derr = (cal / cal.sum(axis=1)[:, None]).std(axis=0) / np.sqrt(len(cal))
print('P5 data shares (South end):', np.round(dshare, 3), ' k5 %.4f kbb %.5f, normalised to the mean of the %s GeV simulations' % (k5, kbb, '/'.join('%.0f' % (e / 1000.) for e in sorted(share))))
for e in sorted(share): print('  sim %.0f GeV South-end shares:' % e, np.round(share[e], 3))
setup(); ROOT.gStyle.SetPadLeftMargin(0.13); ROOT.gStyle.SetPadRightMargin(0.04); ROOT.gStyle.SetPadTopMargin(0.05); ROOT.gStyle.SetPadBottomMargin(0.13); ROOT.gROOT.ForceStyle()
c = ROOT.TCanvas('c', '', 1000, 750); keep = []
fr = ROOT.gPad.DrawFrame(0.5, 0, 6.5, 0.65, ';;South-end share of the visible energy of the stack'); fr.GetXaxis().SetNdivisions(6); keep.append(fr)
for b, lab in enumerate(['SFIL1', 'SFIL2', 'SFIL3', 'SFIL4', 'SFIL5', 'Baby BCAL']): fr.GetXaxis().SetBinLabel(fr.GetXaxis().FindBin(b + 1), lab)
fr.GetXaxis().SetLabelSize(0.045); fr.GetXaxis().LabelsOption('h'); fr.GetYaxis().SetTitleOffset(1.3); fr.GetYaxis().SetTitleSize(0.045); fr.GetYaxis().SetLabelSize(0.04)
COL = {3000: 2, 4000: 4, 5000: 8, 6000: 6}; LS = {3000: 1, 4000: 2, 5000: 7, 6000: 9}; x6 = np.arange(1, 7, dtype=float)
lsim = ROOT.TLegend(0.16, 0.7, 0.45, 0.93); lsim.SetBorderSize(0); lsim.SetFillStyle(0); lsim.SetTextSize(0.034); lsim.SetHeader('simulation, e^{+} under P5')
for e in sorted(share):
    g = ROOT.TGraph(6, x6, np.array(share[e])); g.SetLineColor(COL.get(int(e), 1)); g.SetLineWidth(3); g.SetLineStyle(LS.get(int(e), 1)); g.Draw('L'); keep.append(g); lsim.AddEntry(g, '%.0f GeV' % (e / 1000.), 'l')
ldat = ROOT.TLegend(0.52, 0.78, 0.94, 0.93); ldat.SetBorderSize(0); ldat.SetFillStyle(0); ldat.SetTextSize(0.034); ldat.SetHeader('run 3062, paddle P5 (%d showers)' % int(D['n']))
gd = ROOT.TGraphErrors(6, x6, dshare, np.zeros(6), derr); gd.SetMarkerStyle(20); gd.SetMarkerSize(2.0); gd.SetLineWidth(3); gd.Draw('LP'); keep.append(gd); ldat.AddEntry(gd, 'South-end amplitudes', 'lp')
lsim.Draw(); ldat.Draw(); keep += [lsim, ldat]
t = ROOT.TLatex(); t.SetNDC(); t.SetTextFont(42); t.SetTextSize(0.03); t.SetTextColor(ROOT.kGray + 2)
t.DrawLatex(0.16, 0.64, 'SFIL1-4: South-end peak ADC per cell, uncalibrated, beam row saturates'); t.DrawLatex(0.16, 0.60, 'SFIL5 (KCU1) and Baby BCAL (CAEN) scales set to the 5-6 GeV simulation mean,'); t.DrawLatex(0.16, 0.56, 'so those two data points are a normalisation, not a measurement'); keep.append(t)
c.SaveAs('stack_shares_p5.png'); print('wrote stack_shares_p5.png')
