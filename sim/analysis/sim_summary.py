"""Summary of the bic_sim energy scan: sampling fraction, resolution, longitudinal profile, N-S timing, AstroPix multiplicity.
    python3 sim_summary.py [tag=_v4] /sphenix/tg/tg01/jets/tmengel/eic/sim/bic_sim_e*_s*.root   -> sim_summary<tag>.png, sim_summary<tag>.txt
"""
import sys, glob, numpy as np, uproot, ROOT
sys.path.insert(0, '/sphenix/user/tmengel/eic/bic_testbeam/common'); from plotstyle import setup, finish, label
TAG = ""
if sys.argv[1:] and sys.argv[1].startswith("tag="): TAG = sys.argv[1][4:]; sys.argv = sys.argv[:1] + sys.argv[2:]   # tag=_v4: output names sim_summary_v4.*
files = sorted(sys.argv[1:]) or sorted([f for f in glob.glob('/sphenix/tg/tg01/jets/tmengel/eic/sim/bic_sim_e*_s*.root') if '_x' not in f])
rows = []; prof = {}; share = {}
for f in files:
    a = uproot.open(f)['sim'].arrays(['ebeam', 'sfil_evis', 'sfil_etot', 'bb_evis', 'bb_etot', 'e_escape', 'sfil_tn', 'sfil_ts', 'sfil_x', 'astro_n', 'sfil_e', 'bb_e'], library='np')
    E = a['ebeam'].mean(); vis = a['sfil_evis'][:, :5, :].sum(axis=(1, 2)); vbb = a['bb_evis'].sum(axis=(1, 2)); tot = a['sfil_etot'] + a['bb_etot']
    v = vis + vbb; mu, sig = v.mean(), v.std()
    core = v[(v > mu - 2 * sig) & (v < mu + 2 * sig)]; mu, sig = core.mean(), core.std()
    d = a['sfil_tn'][:, 2, 2] - a['sfil_ts'][:, 2, 2]; ok = (a['sfil_tn'][:, 2, 2] > 0) & (a['sfil_ts'][:, 2, 2] > 0)
    rows.append((E, len(v), vis.mean() / E, vbb.mean() / E, vis.mean() / a['sfil_etot'].mean(), (a['sfil_etot'] + a['bb_etot']).mean() / E, a['e_escape'].mean() / E, mu, sig / mu, d[ok].std(), a['astro_n'].mean()))
    prof[E] = (a['sfil_e'][:, :5].mean(axis=0) / E, a['bb_e'][:, :12].mean(axis=0) / E)
    # share of the SFIL1-4 visible energy per layer, as measured in the data (sqrt(S x N) per cell summed over the layer)
    # shares over the whole stack: SFIL1-5 and the Baby BCAL (visible fibre energy), plus the SFIL1-4-only normalisation
    l5 = a['sfil_evis'][:, :5, :].sum(axis=2); bbv = a['bb_evis'].sum(axis=(1, 2)); tot6 = l5.sum(axis=1) + bbv; ok = tot6 > 0
    share[E] = np.concatenate([l5[ok] / tot6[ok][:, None], (bbv[ok] / tot6[ok])[:, None]], axis=1).mean(axis=0)   # 6 entries: SFIL1..5, BBCAL
rows.sort(); out = open('sim_summary%s.txt' % TAG, 'w')
hdr = 'E[MeV]  n   SFILvis/E  BBvis/E  SF(SFIL)  (SFIL+BB)/E  escape/E  vis_mu[MeV]  sigma/mu  tN-tS rms[ns]  astro_n'
print(hdr); out.write(hdr + '\n')
for r in rows: line = '%6.0f %5d  %.4f    %.4f   %.4f   %.3f        %.3f     %8.1f   %.4f    %.3f       %.1f' % r; print(line); out.write(line + '\n')
# ---- figure: 2 x 2, legends inside the frames
E = np.array([r[0] for r in rows]) / 1000.; res = np.array([r[8] for r in rows]); mus = np.array([r[7] for r in rows])
SHOW = [e for e in sorted(prof) if 2500 <= e <= 6500]     # curves: 3-6 GeV, the Hall D positron range
COL = {3000: 2, 4000: 4, 5000: 8, 6000: 6}; MK = {3000: 20, 4000: 21, 5000: 22, 6000: 23}
def legend(x1, y1, x2, y2):
    l = ROOT.TLegend(x1, y1, x2, y2); l.SetBorderSize(0); l.SetFillStyle(0); l.SetTextSize(0.042); return l
setup(); ROOT.gStyle.SetPadTopMargin(0.06); ROOT.gStyle.SetPadRightMargin(0.04); ROOT.gStyle.SetPadLeftMargin(0.15); ROOT.gStyle.SetPadBottomMargin(0.14); ROOT.gROOT.ForceStyle()
c = ROOT.TCanvas('c', '', 1400, 1100); c.Divide(2, 2, 0.005, 0.005); keep = []
# (1) resolution
c.cd(1); fr = ROOT.gPad.DrawFrame(0, 0, 7, 0.11, ';beam energy [GeV];#sigma/#mu of the visible energy'); keep.append(fr)
g = ROOT.TGraph(len(E), E, res); g.SetMarkerStyle(20); g.SetMarkerSize(1.4); g.Draw('P'); keep.append(g)
fn = ROOT.TF1('fn', 'sqrt([0]*[0]/x + [1]*[1])', 0.3, 7); fn.SetParameters(0.1, 0.02); fn.SetLineColor(2); g.Fit(fn, 'Q'); fn.Draw('same'); keep.append(fn)
l1 = legend(0.45, 0.68, 0.92, 0.88); l1.AddEntry(g, 'e^{+}, visible fibre energy', 'p'); l1.AddEntry(fn, '%.1f%% / #sqrt{E} #oplus %.1f%%' % (100 * abs(fn.GetParameter(0)), 100 * abs(fn.GetParameter(1))), 'l'); l1.Draw(); keep.append(l1)
# (2) linearity
c.cd(2); fr = ROOT.gPad.DrawFrame(0, 0, 7, 500, ';beam energy [GeV];mean visible energy [MeV]'); keep.append(fr)
lin = ROOT.TGraph(len(E), E, mus); lin.SetMarkerStyle(20); lin.SetMarkerSize(1.4); lin.Draw('P'); keep.append(lin)
fl = ROOT.TF1('fl', '[0]*x', 0, 7); fl.SetLineColor(2); lin.Fit(fl, 'Q'); fl.Draw('same'); keep.append(fl)
l2 = legend(0.2, 0.7, 0.7, 0.88); l2.AddEntry(lin, 'e^{+}, mean visible energy', 'p'); l2.AddEntry(fl, 'linear: %.1f MeV per GeV' % fl.GetParameter(0), 'l'); l2.Draw(); keep.append(l2)
# (3) longitudinal profile
c.cd(3); fr = ROOT.gPad.DrawFrame(0.5, 0, 17.5, 0.24, ';layer (1-5 SFIL, 6-17 Baby BCAL);fraction of the beam energy per layer'); keep.append(fr)
l3 = legend(0.55, 0.62, 0.92, 0.88); l3.SetHeader('longitudinal profile')
for e in SHOW:
    ps, pb = prof[e]; y = np.concatenate([ps, pb]); x = np.arange(len(y), dtype=float) + 1
    gg = ROOT.TGraph(len(y), x, y); gg.SetMarkerStyle(MK[e]); gg.SetMarkerSize(1.2); gg.SetMarkerColor(COL[e]); gg.SetLineColor(COL[e]); gg.Draw('LP'); keep.append(gg); l3.AddEntry(gg, 'e^{+} %.0f GeV' % (e / 1000.), 'lp')
l3.Draw(); keep.append(l3)
# (4) layer shares over the whole stack (SFIL1-5 + Baby BCAL), simulation and run 3062
# data_shares.npz (made in bic_testbeam, see README): per-event shares of sqrt(S x N) SFIL1-4, KCU1 SFIL5 arrays B+C+D and the
# CAEN Baby BCAL charge, for triple-matched showering triggers; the KCU1 and CAEN scales are fixed so that the all-paddle
# means of SFIL5 and BBCAL relative to SFIL1-4 equal the 4 GeV simulation (those two points are a normalisation, not a measurement)
D = np.load('data_shares.npz')
c.cd(4); fr = ROOT.gPad.DrawFrame(0.5, 0, 6.5, 0.6, ';;share of the visible energy of the stack'); fr.GetXaxis().SetNdivisions(6); keep.append(fr)
for b, lab in enumerate(['SFIL1', 'SFIL2', 'SFIL3', 'SFIL4', 'SFIL5', 'BBCAL']): fr.GetXaxis().SetBinLabel(fr.GetXaxis().FindBin(b + 1), lab)
fr.GetXaxis().SetLabelSize(0.055); fr.GetXaxis().LabelsOption('h')
l4 = legend(0.42, 0.47, 0.93, 0.9); l4.SetTextSize(0.033); l4.SetHeader('stack = SFIL1-5 + Baby BCAL')
for e in SHOW:
    gg = ROOT.TGraph(6, np.arange(1, 7, dtype=float), np.array(share[e])); gg.SetLineColor(COL[e]); gg.SetLineWidth(3); gg.SetLineStyle({3000: 1, 4000: 2, 5000: 7, 6000: 9}[e]); gg.Draw('L'); keep.append(gg); l4.AddEntry(gg, 'sim e^{+} %.0f GeV' % (e / 1000.), 'l')
x6 = np.arange(1, 7, dtype=float)
gd = ROOT.TGraph(6, x6, np.array(D['all'])); gd.SetMarkerStyle(20); gd.SetMarkerSize(1.8); gd.SetLineWidth(3); gd.Draw('LP'); keep.append(gd); l4.AddEntry(gd, 'run 3062, all paddles (%d triggers)' % int(D['n']), 'lp')
g4 = ROOT.TGraph(6, x6 - 0.12, np.array(D['P4'])); g4.SetMarkerStyle(24); g4.SetMarkerSize(1.5); g4.Draw('P'); keep.append(g4); l4.AddEntry(g4, 'run 3062, paddle P4 (North end)', 'p')
g5 = ROOT.TGraph(6, x6 + 0.12, np.array(D['P5'])); g5.SetMarkerStyle(25); g5.SetMarkerSize(1.5); g5.Draw('P'); keep.append(g5); l4.AddEntry(g5, 'run 3062, paddle P5 (South end)', 'p')
l4.AddEntry(gd, 'SFIL5, BBCAL scale: set by the 4 GeV sim', ''); l4.AddEntry(gd, 'SFIL1-4: uncalibrated, saturating', '')
l4.Draw(); keep.append(l4)
c.SaveAs('sim_summary%s.png' % TAG); print('wrote sim_summary%s.png/.txt' % TAG)
