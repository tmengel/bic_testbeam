"""Five-element shares, SFIL1-5 South, with SFIL5 brought onto the KCU0 scale by the measured SFIL5/SFIL4 ratio.

Sample: showers (South total > 600 ADC) whose trigger also has the KCU1 stream (join through the KCU0 timestamp).
SFIL1-4: template-fit South layer sums, SFIL2 scaled by the photostatistics gain ratio (9.47/14.09), SFIL1 and SFIL4 at the SFIL3 gain.
SFIL5: KCU1 South arrays B+C+D divided by g5/g4 = (<A5>/<A4>) / (model L5/L4 at the anchor energy), so SFIL5 is in SFIL4-equivalent
ADC; by construction <SFIL5>/<SFIL4> then equals the model ratio at the anchor, and the bar spans the 3-6 GeV range of that ratio.
Model: South light per layer, cells 1-3 (cell 4 is empty in the data).
    python3 stack_shares_sfil5.py [simdir] [tag] [anchor_GeV]  -> stack_shares_sfil5_pe<tag>.pdf/.png/.txt
"""
import sys, glob, numpy as np, uproot, ROOT
sys.path.insert(0, '/sphenix/user/tmengel/eic/bic_testbeam/common'); from plotstyle import setup, legend, label
SIMDIR = sys.argv[1] if len(sys.argv) > 1 else '/sphenix/tg/tg01/jets/tmengel/eic/sim/v3_singleexp'
TAG = sys.argv[2] if len(sys.argv) > 2 else ''
ANCHOR = int(float(sys.argv[3]) * 1000) if len(sys.argv) > 3 else 4000
NL = 5; GF4 = np.array([1.0, 9.47 / 14.09, 1.0, 1.0])
out = open('stack_shares_sfil5_pe%s.txt' % TAG, 'w')
def say(*a): print(*a); print(*a, file=out)
# ---- data: join the fit cache to the event table for SFIL5
D = np.load('data_fit_layers.npz'); S = D['S'].astype(float); trig = D['trig']; pad = D['pad']; tot = S.sum(axis=1)
ts0 = uproot.open('/sphenix/user/tmengel/eic/bic_testbeam/data/kcu0_3062.root')['header']['ts0'].array(library='np')
ev = uproot.open('/sphenix/user/tmengel/eic/bic_testbeam/data/events_3062.root')['events'].arrays(['trig_ts', 'has_kcu1', 'sfil5_s'], library='np')
order = np.argsort(ev['trig_ts']); ets = ev['trig_ts'][order]; pos = np.clip(np.searchsorted(ets, ts0[trig]), 0, len(ets) - 1); row = order[pos]; has5 = (ets[pos] == ts0[trig]) & (ev['has_kcu1'][row] == 1)
s5 = np.clip(ev['sfil5_s'][row].astype(float), 0, None); A5 = s5[:, 18:35].sum(axis=1) + s5[:, 36:53].sum(axis=1) + s5[:, 54:71].sum(axis=1)
REC = 'rec' in sys.argv[1:]
if REC:   # recovered SFIL5 sums (kcu1_sat_recovery.py: TOT / neighbour predictor for saturated single SiPMs), keyed by event row
    Rz = np.load('sfil5_recovered.npz'); arr = np.full(len(ets), np.nan); arr[Rz['row']] = Rz['A5']; A5 = arr[row]; has5 = has5 & np.isfinite(A5)
shower = (tot > 600.) & has5
# ---- model
sim, simx80 = {}, {}
for f in sorted(glob.glob(SIMDIR + '/bic_sim_e*_s*.root')):
    g = uproot.open(f)['sim'].arrays(['ebeam', 'sfil_ls'], library='np'); E = int(round(g['ebeam'].mean()))
    if not 2500 <= E <= 6500: continue
    (simx80 if '_x80' in f else sim)[E] = g['sfil_ls'][:, :NL, :3].sum(axis=2).mean(axis=0)
E6 = sorted(sim); Es = np.array(E6, float)
def l54(E, p=None): m = simx80[E] if (p == 5 and E in simx80) else sim[E]; return m[4] / m[3]
say('simulation', SIMDIR, ' anchor %d MeV; matched showers %d' % (ANCHOR, shower.sum()))
# ---- g5/g4 per sample from the sample's own means, at the anchor; SFIL5 in SFIL4-equivalent ADC
def shares(m, E=ANCHOR, p=None):
    a14 = S[m].mean(axis=0) * GF4; a5 = A5[m].mean(); g54 = (a5 / S[m, 3].mean()) / l54(E, p)
    v = np.concatenate([a14, [a5 / g54]]); return v / v.sum(), g54
dall, g54 = shares(shower); say('all paddles: g5/g4 at the anchor = %.2f; shares SFIL1-5 =' % g54, dall.round(3))
dp = {}
for p in (3, 4, 5, 6): dp[p], gp = shares(shower & (pad == p), p=p); say('  P%d: g5/g4 %.2f, shares' % (p, gp), dp[p].round(3))
for E in E6: say('  sim %d GeV:' % (E // 1000), (sim[E] / sim[E].sum()).round(3))
# ---- bars: inner = photostat ratio error (+) 10 % on SFIL2, and the 3-6 GeV range of the SFIL5 anchor; outer adds +-17.5 % on SFIL1 and SFIL4
v0 = np.concatenate([S[shower].mean(axis=0) * GF4, [A5[shower].mean() / g54]])
def dsh(fac): v = v0 * fac; return v / v.sum() - dall
e2 = np.sqrt((0.16 / 14.09) ** 2 + (0.06 / 9.47) ** 2 + 0.10 ** 2)
inner = np.sqrt(0.5 * (dsh(np.array([1, 1 + e2, 1, 1, 1])) ** 2 + dsh(np.array([1, 1 - e2, 1, 1, 1])) ** 2))
r = np.array([l54(ANCHOR) / l54(E) for E in E6]); lo5, hi5 = dsh(np.array([1, 1, 1, 1, r.min()])), dsh(np.array([1, 1, 1, 1, r.max()]))
inner_lo = np.sqrt(inner ** 2 + np.minimum(lo5, hi5).clip(max=0) ** 2); inner_hi = np.sqrt(inner ** 2 + np.maximum(lo5, hi5).clip(min=0) ** 2)
outer = np.sqrt(sum(0.5 * (dsh(np.exp(0.175 * np.eye(NL)[L])) ** 2 + dsh(np.exp(-0.175 * np.eye(NL)[L])) ** 2) for L in (0, 3)))
say('inner bar -/+:', inner_lo.round(3), inner_hi.round(3), ' outer add:', outer.round(3)); out.close()
# ---- figure
setup(); ROOT.gStyle.SetPadLeftMargin(0.14); ROOT.gStyle.SetPadRightMargin(0.04); ROOT.gStyle.SetPadTopMargin(0.08); ROOT.gStyle.SetPadBottomMargin(0.14); ROOT.gROOT.ForceStyle()
c = ROOT.TCanvas('c', '', 1000, 750); keep = []; YMAX = 0.62
fr = ROOT.TH1F('fr', '', NL, 0.5, NL + 0.5); fr.SetMinimum(0); fr.SetMaximum(YMAX); fr.SetStats(0); fr.Draw('AXIS'); keep.append(fr)
ax, ay = fr.GetXaxis(), fr.GetYaxis()
for b in range(NL): ax.SetBinLabel(b + 1, 'SFIL%d' % (b + 1))
ax.SetLabelSize(0.06); ax.SetLabelOffset(0.015); ax.SetTickLength(0)
ay.SetTitle('share of the visible energy in SFIL1-5'); ay.SetTitleSize(0.05); ay.SetTitleOffset(1.25); ay.SetLabelSize(0.045); ay.SetNdivisions(505)
tick = 0.03 * YMAX
for x in range(1, NL + 1):
    for y0, y1 in ((0, tick), (YMAX, YMAX - tick)): ln = ROOT.TLine(x, y0, x, y1); ln.SetLineWidth(2); ln.Draw(); keep.append(ln)
xs = np.arange(1, NL + 1, dtype=float); z = np.zeros(NL)
lo = np.min([dp[p] for p in dp], axis=0); hi = np.max([dp[p] for p in dp], axis=0)
band = ROOT.TGraphAsymmErrors(NL, xs, dall, z, z, dall - lo, hi - dall); band.SetFillColorAlpha(ROOT.kGray + 1, 0.35); band.SetLineWidth(0); band.Draw('3'); keep.append(band)
COL = {3000: ROOT.kRed + 1, 4000: ROOT.kAzure + 2, 5000: ROOT.kGreen + 2, 6000: ROOT.kViolet - 1}; LS = {3000: 1, 4000: 2, 5000: 7, 6000: 9}
lsim = legend(0.17, 0.66, 0.42, 0.90, size=0.04); lsim.SetHeader('simulation, e^{+}')
for E in E6:
    g = ROOT.TGraph(NL, xs, sim[E] / sim[E].sum()); g.SetLineColor(COL.get(E, 1)); g.SetLineWidth(4); g.SetLineStyle(LS.get(E, 1)); g.Draw('L'); keep.append(g); lsim.AddEntry(g, '%.0f GeV' % (E / 1000.), 'l')
go = ROOT.TGraphAsymmErrors(NL, xs, dall, z, z, np.sqrt(inner_lo ** 2 + outer ** 2), np.sqrt(inner_hi ** 2 + outer ** 2)); go.SetLineWidth(2); go.SetLineColor(ROOT.kGray + 3); go.SetMarkerSize(0); go.Draw('E'); keep.append(go)
gi = ROOT.TGraphAsymmErrors(NL, xs, dall, 0.12 + z, 0.12 + z, inner_lo, inner_hi); gi.SetLineWidth(3); gi.SetMarkerStyle(20); gi.SetMarkerSize(2.2); gi.Draw('P'); keep.append(gi)
gl = ROOT.TGraph(NL, xs, dall); gl.SetLineWidth(3); gl.Draw('L'); keep.append(gl)
ldat = legend(0.44, 0.64, 0.94, 0.90, size=0.031); ldat.SetHeader('run 3062, e^{+}, %dk showers with KCU1, South' % round(shower.sum() / 1000.))
ldat.AddEntry(gi, 'SFIL1-4 fitted; SFIL5 / (g_{5}/g_{4} = %.1f)' % g54, 'lp'); ldat.AddEntry(go, 'inner: ratio errors; outer: + SFIL1/4 gains', 'l'); ldat.AddEntry(band, 'spread of paddles P3-P6', 'f')
lsim.Draw(); ldat.Draw(); keep += [lsim, ldat]; keep.append(label('BIC test beam, JLab Hall D, August 2026', size=0.045))
t = ROOT.TLatex(); t.SetNDC(); t.SetTextFont(42); t.SetTextSize(0.028); t.SetTextColor(ROOT.kGray + 2); t.SetTextAlign(31)
t.DrawLatex(0.94, 0.60, 'SFIL5 = KCU1 arrays B+C+D / (g_{5}/g_{4}), g_{5}/g_{4} from the model light ratio L_{5}/L_{4} at %.0f GeV' % (ANCHOR / 1000.)); keep.append(t)
ROOT.gPad.RedrawAxis(); c.SaveAs('stack_shares_sfil5_pe%s.pdf' % TAG); c.SaveAs('stack_shares_sfil5_pe%s.png' % TAG); print('wrote stack_shares_sfil5_pe%s.pdf' % TAG)
