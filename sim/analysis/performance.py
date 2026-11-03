"""Physics performance of run 3062 at the present calibration level: energy resolution and linearity per paddle, paddle
energies from the longitudinal profile, containment, against the v3 model.
Energy per shower in visible MeV with the 4 GeV model-anchored gains: SFIL1-4 South fitted layer sums / (8.19, 23.82, 18.53,
15.01 ADC/MeV), SFIL5 South rebuilt sum / (14.0 x 15.01), Baby BCAL South = sum over live channels of q_c / g_c (bbcal_gains_yfine).
Model: South light of SFIL1-5 (cells 1-3) + Baby BCAL per event, 3-6 GeV (x = 30 mm).  The paddle energy is fitted from the
six-element shares (model shares interpolated in energy); the resolution is the Gaussian core of the total.
    python3 performance.py -> performance.txt/.png, performance.npz
"""
import numpy as np, uproot, glob, sys, ROOT
sys.path.insert(0, '/sphenix/user/tmengel/eic/bic_testbeam/common'); from plotstyle import setup, legend, label
setup(); ROOT.gROOT.SetBatch(True); out = open('performance.txt', 'w')
def say(*a): print(*a, flush=True); print(*a, file=out)
# ---- data
D = np.load('data_fit_layers.npz'); S = D['S'].astype(float); pad = D['pad']; trig = D['trig']; satS = D['satS']
ts0 = uproot.open('/sphenix/user/tmengel/eic/bic_testbeam/data/kcu0_3062.root')['header']['ts0'].array(library='np')[trig]
ev = uproot.open('/sphenix/user/tmengel/eic/bic_testbeam/data/events_3062.root')['events'].arrays(['trig_ts', 'has_kcu1', 'has_caen'], library='np')
ets = ev['trig_ts']; order = np.argsort(ets); ets_s = ets[order]; pos = np.clip(np.searchsorted(ets_s, ts0), 0, len(ets_s) - 1); hit = ets_s[pos] == ts0; row = np.where(hit, order[pos], -1)
N = len(ets)
R5 = np.load('sfil5_recovered.npz'); A5row = np.full(N, np.nan); A5row[R5['row']] = R5['A5']; A5 = A5row[np.maximum(row, 0)]; A5[row < 0] = np.nan
NSEG = [687701, 694330, 125979]; OFF = np.cumsum([0] + NSEG); qn = np.full((OFF[-1], 16), np.nan, np.float32)
for s_ in (0, 1, 2): qn[OFF[s_]:OFF[s_ + 1]] = np.load('/sphenix/tg/tg01/jets/tmengel/eic/caen/calib/caen_calib_seg%d.npz' % s_)['qn']
ic = np.load('/sphenix/user/tmengel/eic/bic_testbeam/sync/out/triggers_3062.npz')['i_caen']
labels = ['A0', 'A1', 'A3', 'A4', 'A2', 'A5', 'A6', 'A7', 'A8', 'A9', 'A14', 'A11']
G = np.load('bbcal_gains_yfine.npz', allow_pickle=True); gain = {str(l): g for l, g in zip(G['labels'], G['gain'])}
w = np.array([1. / gain[l] if l in gain else 0. for l in labels])        # MeV per CAEN unit per channel, 0 for dead
Ebb_row = np.full(N, np.nan); hc = ic >= 0; Ebb_row[hc] = (np.clip(qn[ic[hc]][:, :12], 0, None) * w[None, :]).sum(axis=1)
Ebb = Ebb_row[np.maximum(row, 0)]; Ebb[row < 0] = np.nan
g14 = np.array([8.19, 23.82, 18.53, 15.01]); g5 = 14.0 * 15.01
PXm = {4: -80., 3: -30., 6: 30., 5: 80.}; xp = np.array([PXm.get(int(p), 0.) for p in pad])
LAMB = np.array([1500., 1329., 1646., 1500.]); LAM5, LAMBB = 1500., 1500.    # measured SFIL2/3; others assumed
att14 = np.exp(-xp[:, None] / LAMB[None, :]); att5 = np.exp(-xp / LAM5); attbb = np.exp(-xp / LAMBB)   # South light grows with x: divide it out
E14 = S / g14[None, :] * att14; E5 = A5 / g5 * att5; Ebb = Ebb * attbb
Etot = E14.sum(axis=1) + E5 + Ebb
say('attenuation correction to x = 0 with the paddle position: lambda_eff %s mm (SFIL1-4), %.0f (SFIL5), %.0f (Baby BCAL); factors at P4/P3/P6/P5: %s' % (LAMB.tolist(), LAM5, LAMBB, [round(float(np.exp(-x / 1500.)), 3) for x in (-80, -30, 30, 80)]))
shower = (S.sum(axis=1) > 600) & np.isfinite(Etot)
say('showers with SFIL1-5 and Baby BCAL: %d of %d fit-cache showers' % (shower.sum(), (S.sum(axis=1) > 600).sum()))
say('gains used (4 GeV anchor): SFIL1-4 %s ADC/MeV, SFIL5 %.0f KCU1 ADC/MeV, Baby BCAL per channel from bbcal_gains_yfine' % (g14.tolist(), g5))
# ---- model
sim = {}
for fn in sorted(glob.glob('/sphenix/tg/tg01/jets/tmengel/eic/sim/bic_sim_e*_s*.root')):
    if '_x80' in fn or '_bby' in fn or '_yb' in fn: continue
    g = uproot.open(fn)['sim'].arrays(['ebeam', 'sfil_ls', 'bb_ls'], library='np'); E = int(round(g['ebeam'].mean()))
    comp = np.concatenate([g['sfil_ls'][:, :5, :3].sum(axis=2), g['bb_ls'][:, :10, :4].sum(axis=(1, 2))[:, None]], axis=1)   # [n, 6] South light
    sim[E] = comp
Es = sorted(sim)
def core(v):   # Gaussian core: iterate mean/sigma within +-2 sigma
    m, s = np.median(v), 1.4826 * np.median(np.abs(v - np.median(v)))
    for _ in range(5):
        sel = np.abs(v - m) < 2 * s; m, s = v[sel].mean(), v[sel].std()
    return m, s
say('\nmodel per energy: mean visible (MeV), core sigma/mu, shares SFIL1-5 + BB:')
msh = {}; mres = {}
for E in Es:
    tot = sim[E].sum(axis=1); m, s = core(tot); msh[E] = sim[E].mean(axis=0) / sim[E].mean(axis=0).sum(); mres[E] = (m, s / m)
    say('  %d MeV: <E_vis> %.1f  sigma/mu %.4f  shares %s' % (E, m, s / m, msh[E].round(3)))
def model_share(E):
    Es_ = np.array(Es); return np.array([np.interp(E, Es_, [msh[e][k] for e in Es]) for k in range(6)])
# ---- per paddle: energy from the shares, resolution of the total, linearity
PX = {4: -80, 3: -30, 6: 30, 5: 80}; res = {}
EGEO = {6: 4000., 5: 4000. * 1.125, 3: 4000. * 1.00}        # P6 is the gain anchor; P5 and P3 from the geometric-mean E(x) curve
def model_share4(E):
    Es_ = np.array(Es); sh = np.array([np.interp(E, Es_, [msh[e][k] for e in Es]) for k in range(1, 5)]); return sh / sh.sum()
say('\nper paddle: n, energy (P6 anchor + geometric-mean curve; P4 from the SFIL2-5 shape), shape-fit energy from the SFIL2-5 shares (10 % systematic), data and model SFIL2-5 shares, <E_vis> (attenuation-corrected), core sigma/mu, model sigma/mu')
Egrid = np.arange(2500, 6501, 25); shapeE = {}
for p in (4, 3, 6, 5):
    m = shower & (pad == p); comp = np.concatenate([E14[m], E5[m, None], Ebb[m, None]], axis=1); dsh6 = comp.mean(axis=0) / comp.mean(axis=0).sum()
    d4 = comp.mean(axis=0)[1:5]; d4 = d4 / d4.sum()
    chi = np.array([(((d4 - model_share4(E)) / (0.10 * model_share4(E) + 0.005)) ** 2).sum() for E in Egrid]); shapeE[p] = Egrid[chi.argmin()]
for p in (4, 3, 6, 5):
    m = shower & (pad == p); comp = np.concatenate([E14[m], E5[m, None], Ebb[m, None]], axis=1); dsh = comp.mean(axis=0) / comp.mean(axis=0).sum()
    Efit = EGEO[p] if p in EGEO else EGEO[3] * shapeE[4] / shapeE[3]
    tot = Etot[m]; mu, sg = core(tot); mr = np.interp(Efit, Es, [mres[e][1] for e in Es]); mmu = np.interp(Efit, Es, [mres[e][0] for e in Es])
    d4 = comp.mean(axis=0)[1:5]; d4 = d4 / d4.sum()
    res[p] = (m.sum(), Efit, shapeE[p], dsh, mu, sg / mu, mr, mmu)
    say('  P%d (x %+d): n %6d  E %.2f GeV (shape fit %.2f)  SFIL2-5 shares data %s model %s  <E_vis> %.1f MeV (model %.1f)  sigma/mu %.4f (model %.4f)' % (p, PX[p], m.sum(), Efit / 1000, shapeE[p] / 1000, d4.round(3), model_share4(Efit).round(3), mu, mmu, sg / mu, mr))
# resolution with SFIL1-4 only (template-recovered KCU0, the best understood part) for comparison
say('core sigma/mu of SFIL1-4 alone per paddle: ' + '  '.join('P%d %.4f' % (p, core(E14[shower & (pad == p)].sum(axis=1))[1] / core(E14[shower & (pad == p)].sum(axis=1))[0]) for p in (4, 3, 6, 5)) + '; model SFIL1-4 alone at 4 GeV %.4f' % (core(sim[4000][:, :4].sum(axis=1))[1] / core(sim[4000][:, :4].sum(axis=1))[0]))
say('\nenergy spread implied per paddle if the model resolution holds: sigma_E/E = sqrt(data^2 - model^2): ' + '  '.join('P%d %.1f %%' % (p, 100 * np.sqrt(max(res[p][5] ** 2 - res[p][6] ** 2, 0))) for p in (4, 3, 6, 5)))
say('linearity: <E_vis>/E_fit per paddle: ' + '  '.join('P%d %.1f MeV/GeV' % (p, res[p][4] / (res[p][1] / 1000)) for p in (4, 3, 6, 5)) + '   (model %.1f MeV/GeV at 4 GeV)' % (mres[4000][0] / 4))
say('shape-fit energy ratios (SFIL2-5 profile): P5/P6 %.3f, P3/P6 %.3f, P4/P6 %.3f  (geometric-mean curve: P5/P6 1.12, P3/P6 1.00)' % (shapeE[5] / shapeE[6], shapeE[3] / shapeE[6], shapeE[4] / shapeE[6]))
say('containment: Baby BCAL share per paddle data %s vs model at E_fit %s; model energy behind the Baby BCAL (escape) 15-16 %%' % ([round(float(res[p][3][5]), 3) for p in (4, 3, 6, 5)], [round(float(model_share(res[p][1])[5]), 3) for p in (4, 3, 6, 5)]))
np.savez('performance.npz', Etot=Etot, pad=pad, shower=shower, E14=E14, E5=E5, Ebb=Ebb)
# ---- figure
c = ROOT.TCanvas('c', '', 2400, 1600); c.Divide(3, 2); keep = []
cols = {4: ROOT.kBlue + 1, 3: ROOT.kGreen + 2, 6: ROOT.kOrange + 7, 5: ROOT.kRed + 1}
c.cd(1); fr = ROOT.gPad.DrawFrame(0, 0, 700, 0.08, ';E_{vis} = SFIL1-5 + Baby BCAL [MeV, 4 GeV anchor];fraction per 10 MeV'); keep.append(fr); lg = legend(0.5, 0.6, 0.88, 0.88)
for p in (4, 3, 6, 5):
    m = shower & (pad == p); h = ROOT.TH1F('h%d' % p, '', 70, 0, 700); [h.Fill(v) for v in Etot[m]]; h.Scale(1. / max(1, h.GetEntries())); h.SetLineColor(cols[p]); h.SetLineWidth(2); h.Draw('hist same'); keep.append(h); lg.AddEntry(h, 'P%d: #sigma/#mu %.1f %%' % (p, 100 * res[p][5]), 'l')
fr.SetMaximum(max(k.GetMaximum() for k in keep[-4:]) * 1.3); lg.Draw(); keep.append(lg)
c.cd(2); fr2 = ROOT.gPad.DrawFrame(2.5, 0, 6.5, 0.12, ';energy [GeV];#sigma/#mu of E_{vis}'); keep.append(fr2)
gm = ROOT.TGraph(len(Es), np.array(Es) / 1000., np.array([mres[e][1] for e in Es])); gm.SetLineColor(ROOT.kRed + 1); gm.SetLineWidth(2); gm.Draw('l same'); keep.append(gm)
gd = ROOT.TGraph(4, np.array([res[p][1] / 1000. for p in (4, 3, 6, 5)]), np.array([res[p][5] for p in (4, 3, 6, 5)])); gd.SetMarkerStyle(20); gd.SetMarkerSize(1.4); gd.Draw('p same'); keep.append(gd)
lg2 = legend(0.4, 0.7, 0.88, 0.88); lg2.AddEntry(gm, 'model, monoenergetic', 'l'); lg2.AddEntry(gd, 'data per paddle at the fitted energy', 'p'); lg2.Draw(); keep.append(lg2)
c.cd(3); fr3 = ROOT.gPad.DrawFrame(2.5, 150, 6.5, 500, ';fitted paddle energy [GeV];<E_{vis}> [MeV]'); keep.append(fr3)
gl = ROOT.TGraph(len(Es), np.array(Es) / 1000., np.array([mres[e][0] for e in Es])); gl.SetLineColor(ROOT.kRed + 1); gl.SetLineWidth(2); gl.Draw('l same'); keep.append(gl)
gld = ROOT.TGraph(4, np.array([res[p][1] / 1000. for p in (4, 3, 6, 5)]), np.array([res[p][4] for p in (4, 3, 6, 5)])); gld.SetMarkerStyle(20); gld.SetMarkerSize(1.4); gld.Draw('p same'); keep.append(gld)
lg3 = legend(0.15, 0.7, 0.6, 0.88); lg3.AddEntry(gl, 'model', 'l'); lg3.AddEntry(gld, 'data (4 GeV-anchored gains)', 'p'); lg3.Draw(); keep.append(lg3)
c.cd(4); fr4 = ROOT.gPad.DrawFrame(0.5, 0, 6.5, 0.5, ';element (SFIL1-5, Baby BCAL);share of E_{vis}'); keep.append(fr4); lg4 = legend(0.15, 0.62, 0.55, 0.88)
for p in (4, 3, 6, 5):
    g_ = ROOT.TGraph(6, np.arange(1, 7, dtype=float), res[p][3]); g_.SetMarkerStyle(20); g_.SetMarkerColor(cols[p]); g_.SetLineColor(cols[p]); g_.Draw('lp same'); keep.append(g_); lg4.AddEntry(g_, 'P%d data' % p, 'lp')
    gm_ = ROOT.TGraph(6, np.arange(1, 7, dtype=float) + 0.1, model_share(res[p][1])); gm_.SetMarkerStyle(24); gm_.SetMarkerColor(cols[p]); gm_.SetLineColor(cols[p]); gm_.SetLineStyle(2); gm_.Draw('lp same'); keep.append(gm_)
lg4.Draw(); keep.append(lg4); keep.append(label('open: model at the fitted energy', x=0.55, y=0.85, size=0.035))
c.cd(5); fr5 = ROOT.gPad.DrawFrame(-100, 3.0, 100, 4.5, ';paddle position x [mm];fitted energy [GeV]'); keep.append(fr5)
ge = ROOT.TGraph(4, np.array([-80., -30., 30., 80.]), np.array([res[p][1] / 1000. for p in (4, 3, 6, 5)])); ge.SetMarkerStyle(20); ge.SetMarkerSize(1.4); ge.Draw('p same'); keep.append(ge)
ge3 = ROOT.TGraph(4, np.array([-80., -30., 30., 80.]), np.array([res[p][2] / 1000. for p in (4, 3, 6, 5)])); ge3.SetMarkerStyle(24); ge3.SetMarkerSize(1.4); ge3.Draw('p same'); keep.append(ge3)
lg5 = legend(0.15, 0.7, 0.6, 0.88); lg5.AddEntry(ge, 'P6 anchor + geometric mean', 'p'); lg5.AddEntry(ge3, 'SFIL2-5 shape fit', 'p'); lg5.Draw(); keep.append(lg5)
c.cd(6); fr6 = ROOT.gPad.DrawFrame(0, 0, 700, 0.1, ';E_{vis} [MeV];fraction'); keep.append(fr6); m = shower & (pad == 6)
h = ROOT.TH1F('hd6', '', 70, 0, 700); [h.Fill(v) for v in Etot[m]]; h.Scale(1. / h.GetEntries()); h.SetLineWidth(2); h.Draw('hist same'); keep.append(h)
Ef = res[6][1]; lo, hi = max(e for e in Es if e <= Ef), min(e for e in Es if e >= Ef); wgt = 0 if hi == lo else (Ef - lo) / (hi - lo)
hm = ROOT.TH1F('hm6', '', 70, 0, 700); [hm.Fill(v, 1 - wgt) for v in sim[lo].sum(axis=1)]; [hm.Fill(v, wgt) for v in sim[hi].sum(axis=1)]; hm.Scale(1. / hm.Integral()); hm.SetLineColor(ROOT.kRed + 1); hm.SetLineWidth(2); hm.Draw('hist same'); keep.append(hm)
fr6.SetMaximum(max(h.GetMaximum(), hm.GetMaximum()) * 1.3); lg6 = legend(0.5, 0.7, 0.88, 0.88); lg6.AddEntry(h, 'P6 data', 'l'); lg6.AddEntry(hm, 'model at %.2f GeV' % (Ef / 1000), 'l'); lg6.Draw(); keep.append(lg6)
c.SaveAs('performance.png'); c.SaveAs('performance.pdf'); say('wrote performance.png')
