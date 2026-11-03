"""Per-layer gains and the beam energy per paddle from two observables per layer: the mean South layer amplitude and the
fraction of showers whose South beam-row channel saturates (the 10-bit ceiling is an absolute ADC reference).

Data (data_fit_layers.npz): showers = South total > TOTCUT; per paddle p and layer L: <A_L> (fitted layer sum), f_pL = P(beam-row
channel at the ceiling), phi_pL = beam-row share of the layer (from the unsaturated part of the sample: S3 / S).
Model (per energy E, beam under P6): per-event South light of the layer L_e (cells 1-3; cell 4 is empty in the data), with the
beam-row channel taken as phi_pL x L_e. Two gain estimates per (p, L, E):
    g_mean = <A_L> / <L_e>      and      g_sat such that P(g_sat x phi x L_e >= ceiling_L) = f_pL,  ceiling_L = 1023 - <pedestal_L>.
They agree only at the right energy if the model's amplitude distribution is right; the crossing gives E per paddle and layer.
    python3 sat_gain.py [simdir] [totcut]  -> sat_gain.txt, sat_gain.png
"""
import sys, glob, numpy as np, uproot, ROOT
sys.path.insert(0, '/sphenix/user/tmengel/eic/bic_testbeam/common'); from plotstyle import setup, legend, label
SIMDIR = sys.argv[1] if len(sys.argv) > 1 else '/sphenix/tg/tg01/jets/tmengel/eic/sim'
TOTCUT = float(sys.argv[2]) if len(sys.argv) > 2 else 600.
PED = {1: 87., 2: 95., 3: 70.}                      # mean pedestal of the South beam-row channels SFIL2-4 (fit tree)
D = np.load('data_fit_layers.npz'); S = D['S'].astype(float); S3 = D['S3'].astype(float); sat = D['satS'] > 0; pad = D['pad']; tot = S.sum(axis=1); shower = tot > TOTCUT
sim = {}
for f in sorted(x for x in glob.glob(SIMDIR + '/bic_sim_e*_s*.root') if '_x' not in x):
    g = uproot.open(f)['sim'].arrays(['ebeam', 'sfil_ls'], library='np'); E = int(round(g['ebeam'].mean()))
    if 2500 <= E <= 6500: sim[E] = g['sfil_ls'][:, :4, :3].sum(axis=2)       # cells 1-3, South light, MeV
E6 = sorted(sim); out = open('sat_gain.txt', 'w')
def say(*a): print(*a); print(*a, file=out)
say('simulation', SIMDIR, ' showers: South total >', TOTCUT, ' model layers use cells 1-3 (cell 4 empty in the data)')
res = {}
for p in (4, 3, 6, 5):
    m = shower & (pad == p); n = m.sum()
    say('\n=== P%d, %d showers' % (p, n))
    say('%-6s %7s %6s %6s | ' % ('layer', '<A>', 'f_sat', 'phi') + ' '.join('%16s' % ('E=%d: g_mean/g_sat' % (E // 1000)) for E in E6))
    for L in (1, 2, 3):
        A = S[m, L].mean(); f = sat[m, L].mean(); uns = m & ~sat[:, L]; phi = (S3[uns, L].sum() / S[uns, L].sum()); ceiling = 1023. - PED[L]
        line = 'SFIL%d %7.1f %6.3f %6.3f | ' % (L + 1, A, f, phi); res[(p, L)] = {}
        for E in E6:
            Le = sim[E][:, L]; g_mean = A / Le.mean()
            thr = np.quantile(Le, 1. - f) if 0 < f < 1 else np.nan; g_sat = ceiling / (phi * thr) if thr > 0 else np.nan   # P(phi g L >= ceiling) = f  ->  g = ceiling / (phi x L_(1-f))
            res[(p, L)][E] = (g_mean, g_sat); line += '%7.2f /%7.2f  ' % (g_mean, g_sat)
        say(line)
    # crossing energy per layer: where ln(g_sat/g_mean) changes sign (linear interpolation), and the gains there
    say('crossing: layer, E* (GeV), g at E*')
    for L in (1, 2, 3):
        Es = np.array(E6, float); r = np.array([np.log(res[(p, L)][E][1] / res[(p, L)][E][0]) for E in E6])
        if np.all(np.isnan(r)): continue
        idx = np.where(np.diff(np.sign(r)) != 0)[0]
        if len(idx): i = idx[0]; Ex = Es[i] + (Es[i + 1] - Es[i]) * (-r[i]) / (r[i + 1] - r[i]); gm = np.interp(Ex, Es, [res[(p, L)][E][0] for E in E6]); say('  SFIL%d  E* = %.2f GeV  g = %.1f ADC/MeV' % (L + 1, Ex / 1e3, gm))
        else: say('  SFIL%d  no crossing in 3-6 GeV (g_sat/g_mean = %s)' % (L + 1, np.exp(r).round(2)))
out.close()
# figure: g_mean and g_sat versus E per paddle and layer
setup(); ROOT.gStyle.SetPadLeftMargin(0.14); ROOT.gStyle.SetPadRightMargin(0.04); ROOT.gStyle.SetPadTopMargin(0.1); ROOT.gStyle.SetPadBottomMargin(0.14); ROOT.gROOT.ForceStyle()
c = ROOT.TCanvas('c', '', 1400, 1000); c.Divide(2, 2, 0.002, 0.002); keep = []
COL = {1: ROOT.kAzure + 2, 2: ROOT.kGreen + 2, 3: ROOT.kViolet - 1}
for i, p in enumerate((4, 3, 6, 5)):
    c.cd(i + 1); fr = ROOT.gPad.DrawFrame(2.5, 0, 6.5, 45, ';beam energy hypothesis [GeV];gain [ADC per MeV of South light]'); fr.GetXaxis().SetNdivisions(505); fr.GetYaxis().SetNdivisions(505); fr.GetYaxis().SetTitleOffset(1.3); keep.append(fr)
    lg = legend(0.17, 0.62, 0.6, 0.88, size=0.04); lg.SetHeader('paddle P%d' % p)
    for L in (1, 2, 3):
        Es = np.array(E6, float) / 1e3; gm = np.array([res[(p, L)][E][0] for E in E6]); gs = np.array([res[(p, L)][E][1] for E in E6])
        a = ROOT.TGraph(len(Es), Es, gm); a.SetLineColor(COL[L]); a.SetLineWidth(3); a.Draw('L'); keep.append(a)
        b = ROOT.TGraph(len(Es), Es, gs); b.SetLineColor(COL[L]); b.SetLineWidth(3); b.SetLineStyle(2); b.Draw('L'); keep.append(b)
        lg.AddEntry(a, 'SFIL%d from the mean amplitude' % (L + 1), 'l'); lg.AddEntry(b, 'SFIL%d from the saturation fraction' % (L + 1), 'l')
    lg.Draw(); keep.append(lg)
c.cd(); keep.append(label('run 3062: the energy at which the mean amplitude and the saturated fraction of a layer imply the same gain', x=0.01, y=0.965, size=0.03))
c.SaveAs('sat_gain.png'); c.SaveAs('sat_gain.pdf'); print('wrote sat_gain.png')
