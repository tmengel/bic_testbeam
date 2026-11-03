"""SFIL1-4 layer shares with the improved data points, the preliminary per-layer gains and an equal-gain energy estimate.

Data: data_fit_layers.npz (extract_fit_layers.py): per trigger the South-end amplitude of each layer summed over its four
cells, template-fit amplitude where a fit was accepted (saturation recovered), raw pedestal-subtracted peak minus the
noise-peak bias elsewhere; showers = South total > TOTCUT ADC.  Simulation: South light per layer (sfil_ls, all cells).
    python3 stack_shares_calib.py [simdir] [tag]
      -> stack_shares_sfil_calib<tag>.pdf/.png   shares, sim 3-6 GeV lines, data all paddles + P3-P6 spread band
         sfil_gains<tag>.pdf/.png, sfil_gain_calibration<tag>.txt   per-layer gain = data / sim South light, per energy
"""
import glob, sys, numpy as np, uproot, ROOT
sys.path.insert(0, '/sphenix/user/tmengel/eic/bic_testbeam/common'); from plotstyle import setup, legend, label
SIMDIR = sys.argv[1] if len(sys.argv) > 1 else '/sphenix/tg/tg01/jets/tmengel/eic/sim/v1_purePb'
TAG = sys.argv[2] if len(sys.argv) > 2 else ''
NL, TOTCUT = 4, 600.
out = open('sfil_gain_calibration%s.txt' % TAG, 'w')
def say(*a): print(*a); print(*a, file=out)
say('simulation:', SIMDIR)

# ---- simulation: mean South light per layer (MeV), by beam position (6: x=30 under P6, 5: x=80 under P5) and energy
sim = {}
for f in sorted(glob.glob(SIMDIR + '/bic_sim_e*_s*.root')):
    g = uproot.open(f)['sim'].arrays(['ebeam', 'sfil_ls'], library='np'); E = int(round(g['ebeam'].mean()))
    if not 2500 <= E <= 6500: continue
    sim[(5 if '_x80' in f else 6, E)] = g['sfil_ls'][:, :NL, :].sum(axis=2).mean(axis=0)
E6 = sorted(E for p, E in sim if p == 6)
say('simulation points (paddle, E):', sorted(sim))

# ---- data
D = np.load('data_fit_layers.npz'); S = D['S']; pad = D['pad']; tot = S.sum(axis=1); shower = tot > TOTCUT
def mean_shares(m): v = S[m].mean(axis=0).astype(float); return v, v / v.sum()
vall, dall = mean_shares(shower); nall = int(shower.sum())
vp, dp = {}, {}
for p in (3, 4, 5, 6): vp[p], dp[p] = mean_shares(shower & (pad == p))
old = np.array(np.load('data_shares.npz')['all'][:NL], dtype=float); old = old / old.sum()
say('\ndata: %d showers (South total > %.0f ADC) of %d fitted triggers' % (nall, TOTCUT, len(tot)))
say('mean South amplitude per layer (ADC), all paddles:', vall.round(1)); say('shares all paddles:', dall.round(3))
for p in (3, 4, 5, 6): say('  P%d: %7d showers, mean' % (p, int((shower & (pad == p)).sum())), vp[p].round(1), 'shares', dp[p].round(3))
say('previous points (peak amplitude sqrt(SxN), beam-row cells, renormalised to SFIL1-4):', old.round(3))

# ---- per-layer gains: data / sim South light, P6 showers vs x=30, P5 showers vs x=80, per energy
say('\nPer-layer gain = <data South amplitude> / <sim South light>  [ADC per MeV]')
say('%-22s' % 'sample, energy' + ''.join('%9s' % ('SFIL%d' % (l + 1)) for l in range(NL)) + '   SFIL2-4 mean (rms/mean)')
gains = {}
for p in (6, 5):
    for E in E6:
        if (p, E) not in sim: continue
        g = vp[p] / sim[(p, E)]; gains[(p, E)] = g
        say('%-22s' % ('P%d, %d GeV' % (p, E // 1000)) + ''.join('%9.2f' % x for x in g) + '   %6.2f  (%.3f)' % (g[1:].mean(), g[1:].std() / g[1:].mean()))
say('SFIL1 / mean(SFIL2-4) gain ratio at 5 GeV, P6: %.2f' % (gains[(6, 5000)][0] / gains[(6, 5000)][1:].mean()) if (6, 5000) in gains else '')

# ---- equal-gain energy estimate from the SFIL2-4 profile (SFIL1 excluded): chi2 of the renormalised shares vs E, interpolated
def shares24(v): w = np.array(v[1:4], dtype=float); return w / w.sum()
Es = np.array(E6, dtype=float); grid = np.arange(Es.min(), Es.max() + 1, 50.)
say('\nEqual-gain energy from the SFIL2-4 share profile (model shares interpolated between the scan points, 2%% share uncertainty)')
energy = {}
for name, v in [('all paddles', vall)] + [('P%d' % p, vp[p]) for p in (3, 4, 5, 6)]:
    ref = 5 if name == 'P5' and all((5, E) in sim for E in E6) else 6; d = shares24(v)   # P5 curves only when the x=80 scan covers every energy
    curves = np.array([shares24(sim[(ref, E)]) for E in E6]); chi = []
    for e in grid:
        mdl = np.array([np.interp(e, Es, curves[:, i]) for i in range(3)]); chi.append(np.sum(((d - mdl) / 0.02) ** 2))
    chi = np.array(chi); i = chi.argmin(); inside = chi < chi.min() + 1
    energy[name] = (grid[i], grid[inside].min(), grid[inside].max(), chi.min())
    say('  %-12s SFIL2-4 shares %s  ->  E = %.1f GeV  (chi2 < min+1: %.1f-%.1f, chi2min %.1f; at the scan edge if %.0f or %.0f)' % (name, d.round(3), grid[i] / 1e3, grid[inside].min() / 1e3, grid[inside].max() / 1e3, chi.min(), Es.min() / 1e3, Es.max() / 1e3))
out.close()

# ---- figure 1: shares
setup(); ROOT.gStyle.SetPadLeftMargin(0.14); ROOT.gStyle.SetPadRightMargin(0.04); ROOT.gStyle.SetPadTopMargin(0.08); ROOT.gStyle.SetPadBottomMargin(0.14); ROOT.gROOT.ForceStyle()
c = ROOT.TCanvas('c', '', 1000, 750); keep = []; YMAX = 0.75
fr = ROOT.TH1F('fr', '', NL, 0.5, NL + 0.5); fr.SetMinimum(0); fr.SetMaximum(YMAX); fr.SetStats(0); fr.Draw('AXIS'); keep.append(fr)
ax, ay = fr.GetXaxis(), fr.GetYaxis()
for b in range(NL): ax.SetBinLabel(b + 1, 'SFIL%d' % (b + 1))
ax.SetLabelSize(0.06); ax.SetLabelOffset(0.015); ax.SetTickLength(0)
ay.SetTitle('share of the visible energy in SFIL1-4'); ay.SetTitleSize(0.05); ay.SetTitleOffset(1.25); ay.SetLabelSize(0.045); ay.SetNdivisions(505)
tick = 0.03 * YMAX
for x in range(1, NL + 1):
    for y0, y1 in ((0, tick), (YMAX, YMAX - tick)): ln = ROOT.TLine(x, y0, x, y1); ln.SetLineWidth(2); ln.Draw(); keep.append(ln)
xs = np.arange(1, NL + 1, dtype=float)
# P3-P6 spread band
lo = np.min([dp[p] for p in dp], axis=0).astype(float); hi = np.max([dp[p] for p in dp], axis=0).astype(float)
band = ROOT.TGraphAsymmErrors(NL, xs, dall, np.zeros(NL), np.zeros(NL), dall - lo, hi - dall); band.SetFillColorAlpha(ROOT.kGray + 1, 0.35); band.SetLineWidth(0); band.Draw('3'); keep.append(band)
COL = {3000: ROOT.kRed + 1, 4000: ROOT.kAzure + 2, 5000: ROOT.kGreen + 2, 6000: ROOT.kViolet - 1}; LS = {3000: 1, 4000: 2, 5000: 7, 6000: 9}
lsim = legend(0.17, 0.66, 0.42, 0.90, size=0.04); lsim.SetHeader('simulation, e^{+}')
for E in E6:
    sh = sim[(6, E)] / sim[(6, E)].sum(); g = ROOT.TGraph(NL, xs, sh); g.SetLineColor(COL.get(E, 1)); g.SetLineWidth(4); g.SetLineStyle(LS.get(E, 1)); g.Draw('L'); keep.append(g); lsim.AddEntry(g, '%.0f GeV' % (E / 1000.), 'l')
ldat = legend(0.44, 0.69, 0.94, 0.90, size=0.04); ldat.SetHeader('run 3062, e^{+}, %dk showers' % round(nall / 1000.))
gd = ROOT.TGraph(NL, xs, dall); gd.SetMarkerStyle(20); gd.SetMarkerSize(2.2); gd.SetLineWidth(3); gd.Draw('LP'); keep.append(gd)
ldat.AddEntry(gd, 'fitted amplitude, South end, all cells', 'lp'); ldat.AddEntry(band, 'spread of paddles P3-P6', 'f')
go = ROOT.TGraph(NL, xs, old); go.SetMarkerStyle(24); go.SetMarkerSize(1.8); go.SetMarkerColor(ROOT.kGray + 2); go.SetLineColor(ROOT.kGray + 2); go.SetLineStyle(3); go.SetLineWidth(2); go.Draw('LP'); keep.append(go)
ldat.AddEntry(go, 'previous: #sqrt{S#timesN} peak, beam row only', 'lp')
lsim.Draw(); ldat.Draw(); keep += [lsim, ldat]
keep.append(label('BIC test beam, JLab Hall D, August 2026', size=0.045))
t = ROOT.TLatex(); t.SetNDC(); t.SetTextFont(42); t.SetTextSize(0.03); t.SetTextColor(ROOT.kGray + 2); t.SetTextAlign(31)
t.DrawLatex(0.94, 0.64, 'data: equal gains assumed, ADC shares; sim: South light'); keep.append(t)
ROOT.gPad.RedrawAxis(); c.SaveAs('stack_shares_sfil_calib%s.pdf' % TAG); c.SaveAs('stack_shares_sfil_calib%s.png' % TAG)

# ---- figure 2: per-layer gains versus energy hypothesis
c2 = ROOT.TCanvas('c2', '', 1000, 750); keep2 = []; GMAX = max(g.max() for g in gains.values()) * 1.9
fr2 = ROOT.TH1F('fr2', '', NL, 0.5, NL + 0.5); fr2.SetMinimum(0); fr2.SetMaximum(GMAX); fr2.SetStats(0); fr2.Draw('AXIS'); keep2.append(fr2)
ax, ay = fr2.GetXaxis(), fr2.GetYaxis()
for b in range(NL): ax.SetBinLabel(b + 1, 'SFIL%d' % (b + 1))
ax.SetLabelSize(0.06); ax.SetLabelOffset(0.015); ax.SetTickLength(0)
ay.SetTitle('gain  [ADC per MeV of South light]'); ay.SetTitleSize(0.05); ay.SetTitleOffset(1.25); ay.SetLabelSize(0.045); ay.SetNdivisions(505)
tick = 0.03 * GMAX
for x in range(1, NL + 1):
    for y0, y1 in ((0, tick), (GMAX, GMAX - tick)): ln = ROOT.TLine(x, y0, x, y1); ln.SetLineWidth(2); ln.Draw(); keep2.append(ln)
lg = legend(0.17, 0.66, 0.50, 0.90, size=0.038); lg.SetHeader('if the beam energy is')
MK = {6: 20, 5: 25}
for p in (6, 5):
    for E in E6:
        if (p, E) not in gains: continue
        g = ROOT.TGraph(NL, xs + (0.08 if p == 5 else -0.08), gains[(p, E)]); g.SetMarkerStyle(MK[p]); g.SetMarkerSize(2.0); g.SetMarkerColor(COL.get(E, 1)); g.SetLineColor(COL.get(E, 1)); g.SetLineWidth(2); g.SetLineStyle(LS.get(E, 1)); g.Draw('LP' if p == 6 else 'P'); keep2.append(g)
        if p == 6: lg.AddEntry(g, '%.0f GeV' % (E / 1000.), 'l')
lm = legend(0.58, 0.76, 0.94, 0.90, size=0.036)
m6 = ROOT.TMarker(0, 0, 20); m6.SetMarkerSize(2.0); m5 = ROOT.TMarker(0, 0, 25); m5.SetMarkerSize(2.0); keep2 += [m6, m5]
lm.AddEntry(m6, 'P6 showers, beam under P6', 'p'); lm.AddEntry(m5, 'P5 showers, beam under P5', 'p')
lg.Draw(); lm.Draw(); keep2 += [lg, lm]; keep2.append(label('BIC test beam, JLab Hall D, August 2026', size=0.045))
t2 = ROOT.TLatex(); t2.SetNDC(); t2.SetTextFont(42); t2.SetTextSize(0.03); t2.SetTextColor(ROOT.kGray + 2); t2.SetTextAlign(31)
t2.DrawLatex(0.94, 0.66, 'South-end amplitude summed over the four cells of a layer'); t2.DrawLatex(0.94, 0.62, 'over the simulated South light of that layer'); keep2.append(t2)
ROOT.gPad.RedrawAxis(); c2.SaveAs('sfil_gains%s.pdf' % TAG); c2.SaveAs('sfil_gains%s.png' % TAG); print('wrote figures')
