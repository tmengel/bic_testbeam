"""SFIL1-4 shares with the route-3 (photostatistics) gain equalisation: SFIL2 and SFIL3 South amplitudes are converted
to photoelectrons with the measured ADC/pe of their beam-row channels (photostat_ns.py: SFIL2 S c3 14.09 +- 0.16,
SFIL3 S c3 9.47 +- 0.06); SFIL1 and SFIL4 have no photostatistics measurement (SFIL1 North masked, SFIL4 North
anomalous) and are ASSUMED to have the SFIL3 gain.  Equal photoelectron yield per MeV in all layers is still assumed.
    python3 stack_shares_pe.py [simdir] [tag]  -> stack_shares_sfil_pe<tag>.pdf/.png/.txt
"""
import glob, sys, numpy as np, uproot, ROOT
sys.path.insert(0, '/sphenix/user/tmengel/eic/bic_testbeam/common'); from plotstyle import setup, legend, label
SIMDIR = sys.argv[1] if len(sys.argv) > 1 else '/sphenix/tg/tg01/jets/tmengel/eic/sim'
TAG = sys.argv[2] if len(sys.argv) > 2 else ''
MODE = sys.argv[3] if len(sys.argv) > 3 else 'fit'      # fit: all showers, fitted amplitudes; nosat: only triggers without any saturated channel; nosat_raw: the same with raw peaks only
NL = 4; A_PE = {1: (14.09, 0.16), 2: (9.47, 0.06)}                     # ADC per photoelectron (x excess noise), South beam-row channels, from photostat_ns.txt
gfac = np.array([A_PE[2][0] / A_PE[2][0], A_PE[2][0] / A_PE[1][0], 1.0, 1.0])   # relative gain correction per layer (SFIL1, SFIL4 assumed = SFIL3)
sim, simcells = {}, {}
for f in sorted(x for x in glob.glob(SIMDIR + '/bic_sim_e*_s*.root') if '_x' not in x):
    g = uproot.open(f)['sim'].arrays(['ebeam', 'sfil_ls'], library='np'); E = int(round(g['ebeam'].mean()))
    if 2500 <= E <= 6500: simcells[E] = g['sfil_ls'][:, :NL, :4]; sim[E] = simcells[E].sum(axis=2).mean(axis=0)
E6 = sorted(sim); Es = np.array(E6, dtype=float)
D = np.load('data_fit_layers.npz'); S = D['S'].astype(float); pad = D['pad']; tot = S.sum(axis=1); shower = tot > 600.
if MODE.startswith('nosat'):
    nosat = (D['nsatfail'] == 0) & (D['Ssat'].sum(axis=1) == 0) & (D['Nsat'].sum(axis=1) == 0)
    print('unsaturated triggers: %d of %d showers (%.1f %%); by paddle P3-P6: %s of %s' % ((shower & nosat).sum(), shower.sum(), 100. * (shower & nosat).sum() / shower.sum(), [int((shower & nosat & (pad == p)).sum()) for p in (3, 4, 5, 6)], [int((shower & (pad == p)).sum()) for p in (3, 4, 5, 6)]))
    print('mean South total: all showers %.0f ADC, unsaturated %.0f ADC' % (tot[shower].mean(), tot[shower & nosat].mean()))
    shower = shower & nosat
    if MODE == 'nosat_raw': S = D['Sraw'].astype(float)
def shares(m, fac=gfac): v = S[m].mean(axis=0) * fac; return v / v.sum()
dall = shares(shower); dp = {p: shares(shower & (pad == p)) for p in (3, 4, 5, 6)}; draw = shares(shower, np.ones(NL))
if MODE.startswith('nosat'):
    # the same cut on the model by QUANTILE MATCHING: the data lose, per layer, the fraction f_L of showers whose South beam-row
    # channel sits at the ceiling; the model drops the same fraction of its events with the largest light in that layer (cells 1-3),
    # layer by layer, so the selection is equivalent whatever the gains (an amplitude threshold needs a gain and the data's
    # amplitude distribution is 1.3 x broader than the model's, so a threshold removes 3 x fewer model events)
    satS = D['satS'] > 0; fL = (satS[(tot > 600.)]).mean(axis=0)
    print('data: fraction of showers with the South beam-row channel saturated, per layer:', fL.round(3), '-> the model drops the same top fractions')
    for E in E6:
        L3 = simcells[E][:, :, :3].sum(axis=2); keep = np.ones(len(L3), bool)
        for L in range(NL):
            if fL[L] > 0: keep &= L3[:, L] < np.quantile(L3[:, L], 1. - fL[L])
        print('sim %d GeV: kept fraction %.2f (data %.2f)' % (E // 1000, keep.mean(), (satS[(tot > 600.)].sum(axis=1) == 0).mean()))
        sim[E] = simcells[E][keep][:, :, :3].sum(axis=2).mean(axis=0)      # cells 1-3: cell 4 is empty in the data
out = open('stack_shares_sfil_pe%s.txt' % TAG, 'w')
def say(*a): print(*a); print(*a, file=out)
say('simulation:', SIMDIR, ' mode:', MODE, ' gain factors applied to the South layer sums (SFIL1-4):', gfac.round(3))
say('equal-gain shares:', draw.round(3)); say('pe-equalised shares, all paddles:', dall.round(3))
for p in (3, 4, 5, 6): say('  P%d:' % p, dp[p].round(3))
for E in E6: say('  sim %d GeV:' % (E // 1000), (sim[E] / sim[E].sum()).round(3))
# uncertainty: photostat errors on the SFIL2/SFIL3 ratio (statistical) and +-10 % for the excess-noise assumption, plus +-17.5 % on the unmeasured SFIL1 and SFIL4 gains
v0 = S[shower].mean(axis=0) * gfac
def dshare(fac): v = v0 * fac; return v / v.sum() - dall
e_ratio = np.sqrt((A_PE[1][1] / A_PE[1][0]) ** 2 + (A_PE[2][1] / A_PE[2][0]) ** 2 + 0.10 ** 2)
inner = np.sqrt(0.5 * (dshare(np.array([1, 1 + e_ratio, 1, 1])) ** 2 + dshare(np.array([1, 1 - e_ratio, 1, 1])) ** 2))
outer = np.sqrt(inner ** 2 + sum(0.5 * (dshare(np.exp(sg * 0.175 * np.eye(NL)[L])) ** 2 + dshare(np.exp(-sg * 0.175 * np.eye(NL)[L])) ** 2) for L in (0, 3) for sg in (1,)))
say('inner bar (SFIL2/SFIL3 gain ratio: stat (+) 10 %):', inner.round(3)); say('outer bar (+ SFIL1, SFIL4 gains unknown, +-17.5 %):', outer.round(3))
# energy from the SFIL2-4 profile
def s24(v): w = np.array(v[1:4], dtype=float); return w / w.sum()
grid = np.arange(Es.min(), Es.max() + 1, 50.); curves = np.array([s24(sim[E]) for E in E6])
for name, v in [('all', S[shower].mean(axis=0) * gfac)] + [('P%d' % p, S[shower & (pad == p)].mean(axis=0) * gfac) for p in (3, 4, 5, 6)]:
    d = s24(v); chi = np.array([np.sum(((d - np.array([np.interp(e, Es, curves[:, i]) for i in range(3)])) / 0.02) ** 2) for e in grid]); i = chi.argmin(); ok = chi < chi.min() + 1
    say('  %-4s SFIL2-4 pe-equalised shares %s -> E = %.1f GeV (%.1f-%.1f), chi2min %.1f' % (name, d.round(3), grid[i] / 1e3, grid[ok].min() / 1e3, grid[ok].max() / 1e3, chi.min()))
out.close()
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
xs = np.arange(1, NL + 1, dtype=float); z = np.zeros(NL)
lo = np.min([dp[p] for p in dp], axis=0); hi = np.max([dp[p] for p in dp], axis=0)
band = ROOT.TGraphAsymmErrors(NL, xs, dall, z, z, dall - lo, hi - dall); band.SetFillColorAlpha(ROOT.kGray + 1, 0.35); band.SetLineWidth(0); band.Draw('3'); keep.append(band)
COL = {3000: ROOT.kRed + 1, 4000: ROOT.kAzure + 2, 5000: ROOT.kGreen + 2, 6000: ROOT.kViolet - 1}; LS = {3000: 1, 4000: 2, 5000: 7, 6000: 9}
lsim = legend(0.17, 0.66, 0.42, 0.90, size=0.04); lsim.SetHeader('simulation, e^{+}' + (', same cut' if MODE.startswith('nosat') else ''))
for E in E6:
    g = ROOT.TGraph(NL, xs, sim[E] / sim[E].sum()); g.SetLineColor(COL.get(E, 1)); g.SetLineWidth(4); g.SetLineStyle(LS.get(E, 1)); g.Draw('L'); keep.append(g); lsim.AddEntry(g, '%.0f GeV' % (E / 1000.), 'l')
go = ROOT.TGraphAsymmErrors(NL, xs, dall, z, z, outer, outer); go.SetLineWidth(2); go.SetLineColor(ROOT.kGray + 3); go.SetMarkerSize(0); go.Draw('E'); keep.append(go)
gi = ROOT.TGraphAsymmErrors(NL, xs, dall, 0.12 + z, 0.12 + z, inner, inner); gi.SetLineWidth(3); gi.SetMarkerStyle(20); gi.SetMarkerSize(2.2); gi.Draw('P'); keep.append(gi)
gl = ROOT.TGraph(NL, xs, dall); gl.SetLineWidth(3); gl.Draw('L'); keep.append(gl)
gr = ROOT.TGraph(NL, xs, draw); gr.SetMarkerStyle(24); gr.SetMarkerSize(1.8); gr.SetMarkerColor(ROOT.kGray + 2); gr.SetLineColor(ROOT.kGray + 2); gr.SetLineStyle(3); gr.SetLineWidth(2); gr.Draw('LP'); keep.append(gr)
ldat = legend(0.42, 0.64, 0.94, 0.90, size=0.029); ldat.SetHeader('run 3062, e^{+}, all paddles, South end, all cells' + (', no saturated channel' if MODE.startswith('nosat') else '') + (', raw peaks' if MODE == 'nosat_raw' else ''))
ldat.AddEntry(gi, 'SFIL2/SFIL3 gain from photostatistics (bar: its error)', 'lp'); ldat.AddEntry(go, 'outer bar: SFIL1, SFIL4 gains assumed = SFIL3 (#pm17.5%)', 'l'); ldat.AddEntry(band, 'spread of paddles P3-P6', 'f'); ldat.AddEntry(gr, 'equal gains (before)', 'lp')
lsim.Draw(); ldat.Draw(); keep += [lsim, ldat]; keep.append(label('BIC test beam, JLab Hall D, August 2026', size=0.045))
ROOT.gPad.RedrawAxis(); c.SaveAs('stack_shares_sfil_pe%s.pdf' % TAG); c.SaveAs('stack_shares_sfil_pe%s.png' % TAG); print('wrote stack_shares_sfil_pe%s.pdf' % TAG)
