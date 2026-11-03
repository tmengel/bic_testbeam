"""Six-element share figure (SFIL1-5 + Baby BCAL) with preliminary gains and systematic uncertainties.

SFIL1-4: template-fitted South-end amplitudes summed over the cells (data_fit_layers.npz), one common gain
(equal gains assumed) = sum of the P6 layer sums 2-4 over the simulated South light of layers 2-4 at the anchor energy.
SFIL5 (KCU1 South arrays B+C+D) and Baby BCAL (CAEN South, 12 channels): means of the triple-matched single-paddle
showers (data_sfil5_bbcal.npz), each with its own gain anchored to the simulation at the anchor energy with P6 showers
versus the beam under P6.  By construction these two points sit on the model at the anchor energy; their bar is the
range of anchors 3-6 GeV.
Systematics on SFIL1-4 (quadrature, inner bar): shower threshold (300 / 1200 ADC instead of 600), saturation recovery
(fitted amplitude of saturated channels +- 10 %), noise-peak bias (+- 1 ADC per raw-peak channel), statistics.
Gain uniformity (outer bar): the equal-gain assumption, taken as a per-layer gain uncertainty equal to the rms of
ln(North/South) over SFIL2-4 (the two ends of the same layer differ by that much), propagated to the shares.
Grey band: spread of paddles P3-P6 (physics: the beam energy changes across the paddles, not an error).
    python3 stack_shares_calib6.py [simdir] [tag] [anchor_GeV] [pe]   (pe: SFIL2 scaled by the measured SFIL3/SFIL2 gain ratio 9.47/14.09)  -> stack_shares_calib6<tag>.pdf/.png, stack_shares_calib6<tag>.txt
"""
import glob, sys, numpy as np, uproot, ROOT
sys.path.insert(0, '/sphenix/user/tmengel/eic/bic_testbeam/common'); from plotstyle import setup, legend, label
SIMDIR = sys.argv[1] if len(sys.argv) > 1 else '/sphenix/tg/tg01/jets/tmengel/eic/sim/v1_purePb'
TAG = sys.argv[2] if len(sys.argv) > 2 else ''
ANCHOR = int(float(sys.argv[3]) * 1000) if len(sys.argv) > 3 else 4000
PE = len(sys.argv) > 4 and sys.argv[4] == 'pe'          # apply the photostatistics gain equalisation of SFIL2 vs SFIL3 (method 4.13) before the common gain
GFAC = np.array([1.0, 9.47 / 14.09, 1.0, 1.0]) if PE else np.ones(4)   # SFIL1 and SFIL4 assumed at the SFIL3 gain
NE, TOTCUT = 6, 600.; NAMES = ['SFIL1', 'SFIL2', 'SFIL3', 'SFIL4', 'SFIL5', 'Baby BCAL']
out = open('stack_shares_calib6%s.txt' % TAG, 'w')
def say(*a): print(*a); print(*a, file=out)
say('simulation:', SIMDIR, ' anchor energy for the SFIL5 / Baby BCAL gains and the common SFIL1-4 gain: %d MeV' % ANCHOR, ' photostatistics gain factors on SFIL1-4:', GFAC.round(3))

# ---- simulation: mean South light (MeV) of the six elements, by beam position (6: x=30, 5: x=80) and energy
sim = {}
for f in sorted(glob.glob(SIMDIR + '/bic_sim_e*_s*.root')):
    g = uproot.open(f)['sim'].arrays(['ebeam', 'sfil_ls', 'bb_ls'], library='np'); E = int(round(g['ebeam'].mean()))
    if not 2500 <= E <= 6500: continue
    sim[(5 if '_x80' in f else 6, E)] = np.concatenate([g['sfil_ls'][:, :5, :].sum(axis=2).mean(axis=0), [g['bb_ls'].sum(axis=(1, 2)).mean()]])
E6 = sorted(E for p, E in sim if p == 6); Es = np.array(E6, dtype=float)
def simshare(p, E): v = sim[(p, E)]; return v / v.sum()

# ---- data: SFIL1-4 per-event layer sums (full fitted sample), SFIL5 / BBCAL means (triple-matched sample)
D = np.load('data_fit_layers.npz'); S = D['S'].astype(float) * GFAC; Ssat = D['Ssat'].astype(float) * GFAC; Snfit = D['Snfit'].astype(float); pad = D['pad']; tot = D['S'].astype(float).sum(axis=1)
T = np.load('data_sfil5_bbcal.npz')
def layer_means(cut=TOTCUT, satscale=1.0, biasshift=0.0, p=None):
    m = tot > cut
    if p is not None: m &= pad == p
    v = S[m] + (satscale - 1.0) * Ssat[m] + biasshift * Snfit[m]
    return v.mean(axis=0), int(m.sum()), v.std(axis=0) / np.sqrt(m.sum())
def shares6(l14, e5, bb, anchor=ANCHOR, l14_ref=None, e5_ref=None, bb_ref=None):
    """MeV per element and shares; the gains come from the reference (P6) means against the beam under P6 at the anchor"""
    ref = sim[(6, anchor)]
    g24 = l14_ref[1:4].sum() / ref[1:4].sum(); g5 = e5_ref / ref[4]; gbb = bb_ref / ref[5]
    mev = np.concatenate([l14 / g24, [e5 / g5, bb / gbb]]); return mev / mev.sum(), mev, (g24, g5, gbb)
l_all, n_all, l_all_err = layer_means(); l_p6, n_p6, _ = layer_means(p=6)
lp = {p: layer_means(p=p)[0] for p in (3, 4, 5, 6)}
e5 = {k: float(T['e5_' + k]) for k in ('all', 'P3', 'P4', 'P5', 'P6')}; bb = {k: float(T['bb_' + k]) for k in e5}
e5err = {k: float(T['e5_err_' + k]) for k in e5}; bberr = {k: float(T['bb_err_' + k]) for k in e5}
central, mev, gains = shares6(l_all, e5['all'], bb['all'], ANCHOR, l_p6, e5['P6'], bb['P6'])
say('\nSFIL1-4: %d showers (South total > %.0f ADC); SFIL5 / Baby BCAL: %d triple-matched single-paddle showers' % (n_all, TOTCUT, int(T['n_all'])))
say('gains at the anchor: SFIL1-4 common %.2f ADC/MeV (SFIL2-4 vs sim), SFIL5 %.1f KCU1 ADC/MeV, Baby BCAL %.1f CAEN units/MeV' % gains)
say('data MeV per element at the anchor:', mev.round(1), ' sum %.0f MeV (sim at the anchor: %.0f)' % (mev.sum(), sim[(6, ANCHOR)].sum()))
say('shares, all paddles:', central.round(3))
for E in E6: say('  sim %d GeV:' % (E // 1000), simshare(6, E).round(3))
# per paddle (own means for all six elements, gains fixed)
pshare = {}
for p in (3, 4, 5, 6):
    pshare[p] = shares6(lp[p], e5['P%d' % p], bb['P%d' % p], ANCHOR, l_p6, e5['P6'], bb['P6'])[0]; say('  P%d:' % p, pshare[p].round(3))
lo = np.min([pshare[p] for p in pshare], axis=0); hi = np.max([pshare[p] for p in pshare], axis=0)

# ---- systematics
sysd = {}
def var(name, l14=None, e5v=None, bbv=None, anchor=ANCHOR):
    s = shares6(l_all if l14 is None else l14, e5['all'] if e5v is None else e5v, bb['all'] if bbv is None else bbv, anchor, l_p6, e5['P6'], bb['P6'])[0]
    sysd.setdefault(name, []).append(s - central)
var('threshold', l14=layer_means(cut=300.)[0]); var('threshold', l14=layer_means(cut=1200.)[0])
var('saturation', l14=layer_means(satscale=1.1)[0]); var('saturation', l14=layer_means(satscale=0.9)[0])
var('noise bias', l14=layer_means(biasshift=1.0)[0]); var('noise bias', l14=layer_means(biasshift=-1.0)[0])
var('statistics', l14=l_all + l_all_err, e5v=e5['all'] + e5err['all'], bbv=bb['all'] + bberr['all'])
for E in E6:
    if E != ANCHOR: var('anchor energy', anchor=E)
def halfrange(d): d = np.array(d); return np.maximum(d.max(axis=0), 0) if len(d) == 1 else 0.5 * (d.max(axis=0) - d.min(axis=0))
meas = {k: halfrange(v) for k, v in sysd.items() if k != 'anchor energy'}
anchor_lo = -np.minimum(np.array(sysd['anchor energy']).min(axis=0), 0); anchor_hi = np.maximum(np.array(sysd['anchor energy']).max(axis=0), 0)
# gain uniformity: rms of ln(N/S) over SFIL2-4 as the per-layer gain uncertainty, propagated to the six shares (SFIL1-4 only)
N = D['N'].astype(float); m = tot > TOTCUT; ns = N[m].mean(axis=0)[1:4] / S[m].mean(axis=0)[1:4]; sig_g = np.std(np.log(ns))
say('\nNorth/South of the layer sums, SFIL2-4:', ns.round(3), ' -> per-layer gain uncertainty (rms of ln N/S) %.3f' % sig_g)
gain_unc = np.zeros(NE)
for L in range(4):
    for sgn in (+1, -1):
        sg = (0.10 if L == 1 else 0.0 if L == 2 else sig_g) if PE else sig_g
        v = mev.copy(); v[L] *= np.exp(sgn * sg); gain_unc += (v / v.sum() - central) ** 2 / 2   # average of the +/- variations, summed in quadrature over the four layers
gain_unc = np.sqrt(gain_unc)
inner = np.sqrt(sum(meas[k] ** 2 for k in meas)); inner_lo = np.sqrt(inner ** 2 + anchor_lo ** 2); inner_hi = np.sqrt(inner ** 2 + anchor_hi ** 2)
say('\nSystematic uncertainties on the shares (absolute):')
say('%-16s' % 'source' + ''.join('%11s' % n for n in NAMES))
for k in meas: say('%-16s' % k + ''.join('%11.4f' % x for x in meas[k]))
say('%-16s' % 'anchor 3-6 GeV' + ''.join('%11s' % ('-%.3f/+%.3f' % (a, b)) for a, b in zip(anchor_lo, anchor_hi)))
say('%-16s' % 'measurement sum' + ''.join('%11s' % ('-%.3f/+%.3f' % (a, b)) for a, b in zip(inner_lo, inner_hi)))
say('%-16s' % 'gain uniformity' + ''.join('%11.4f' % x for x in gain_unc))
say('%-16s' % 'paddle spread' + ''.join('%11s' % ('-%.3f/+%.3f' % (a, b)) for a, b in zip(central - lo, hi - central)))
out.close()

# ---- figure
setup(); ROOT.gStyle.SetPadLeftMargin(0.13); ROOT.gStyle.SetPadRightMargin(0.04); ROOT.gStyle.SetPadTopMargin(0.08); ROOT.gStyle.SetPadBottomMargin(0.14); ROOT.gROOT.ForceStyle()
c = ROOT.TCanvas('c', '', 1100, 750); keep = []; YMAX = 0.62
fr = ROOT.TH1F('fr', '', NE, 0.5, NE + 0.5); fr.SetMinimum(0); fr.SetMaximum(YMAX); fr.SetStats(0); fr.Draw('AXIS'); keep.append(fr)
ax, ay = fr.GetXaxis(), fr.GetYaxis()
for b in range(NE): ax.SetBinLabel(b + 1, NAMES[b])
ax.SetLabelSize(0.055); ax.SetLabelOffset(0.015); ax.SetTickLength(0)
ay.SetTitle('share of the visible energy'); ay.SetTitleSize(0.05); ay.SetTitleOffset(1.2); ay.SetLabelSize(0.045); ay.SetNdivisions(505)
tick = 0.03 * YMAX
for x in range(1, NE + 1):
    for y0, y1 in ((0, tick), (YMAX, YMAX - tick)): ln = ROOT.TLine(x, y0, x, y1); ln.SetLineWidth(2); ln.Draw(); keep.append(ln)
xs = np.arange(1, NE + 1, dtype=float); z = np.zeros(NE)
band = ROOT.TGraphAsymmErrors(NE, xs, central, z, z, central - lo, hi - central); band.SetFillColorAlpha(ROOT.kGray + 1, 0.35); band.SetLineWidth(0); band.Draw('3'); keep.append(band)
COL = {3000: ROOT.kRed + 1, 4000: ROOT.kAzure + 2, 5000: ROOT.kGreen + 2, 6000: ROOT.kViolet - 1}; LS = {3000: 1, 4000: 2, 5000: 7, 6000: 9}
lsim = legend(0.16, 0.66, 0.40, 0.90, size=0.038); lsim.SetHeader('simulation, e^{+}')
for E in E6:
    g = ROOT.TGraph(NE, xs, simshare(6, E)); g.SetLineColor(COL.get(E, 1)); g.SetLineWidth(4); g.SetLineStyle(LS.get(E, 1)); g.Draw('L'); keep.append(g); lsim.AddEntry(g, '%.0f GeV' % (E / 1000.), 'l')
gouter = ROOT.TGraphAsymmErrors(NE, xs, central, z, z, np.sqrt(inner_lo ** 2 + gain_unc ** 2), np.sqrt(inner_hi ** 2 + gain_unc ** 2)); gouter.SetLineWidth(2); gouter.SetLineColor(ROOT.kGray + 3); gouter.SetMarkerSize(0); gouter.Draw('E'); keep.append(gouter)
ginner = ROOT.TGraphAsymmErrors(NE, xs, central, 0.12 + z, 0.12 + z, inner_lo, inner_hi); ginner.SetLineWidth(3); ginner.SetMarkerStyle(20); ginner.SetMarkerSize(2.2); ginner.Draw('P'); keep.append(ginner)
gline = ROOT.TGraph(NE, xs, central); gline.SetLineWidth(3); gline.Draw('L'); keep.append(gline)
ldat = legend(0.40, 0.66, 0.94, 0.90, size=0.030); ldat.SetHeader('run 3062, e^{+}, all paddles, anchor %.0f GeV' % (ANCHOR / 1000.))
ldat.AddEntry(ginner, 'fitted amplitudes, %s; bar: meas. syst.' % ('SFIL2/3 gain from photostat.' if PE else 'equal gains'), 'lp')
ldat.AddEntry(gouter, 'outer bar: + %s, anchor 3-6 GeV' % ('SFIL1/4 gains unmeasured' if PE else 'gain uniformity'), 'l')
ldat.AddEntry(band, 'spread of paddles P3-P6', 'f')
lsim.Draw(); ldat.Draw(); keep += [lsim, ldat]
keep.append(label('BIC test beam, JLab Hall D, August 2026', size=0.045))
t = ROOT.TLatex(); t.SetNDC(); t.SetTextFont(42); t.SetTextSize(0.028); t.SetTextColor(ROOT.kGray + 2); t.SetTextAlign(31)
t.DrawLatex(0.94, 0.60, 'SFIL5 (KCU1) and Baby BCAL (CAEN): gains anchored to the simulation at %.0f GeV (normalisation)' % (ANCHOR / 1000.)); keep.append(t)
ROOT.gPad.RedrawAxis(); c.SaveAs('stack_shares_calib6%s.pdf' % TAG); c.SaveAs('stack_shares_calib6%s.png' % TAG); print('wrote stack_shares_calib6%s.pdf' % TAG)
