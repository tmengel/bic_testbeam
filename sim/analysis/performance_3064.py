"""Run 3064 (12 h, 18 x 3062) with the KCU0 part only: SFIL1-4 South energy resolution and linearity per paddle, their stability
over the run, and the equal-gain / pe-equalised shape energy, with the same definitions as the run-3062 performance chain
(performance.py: Gaussian core of the total above the shower threshold, iterated +-2 sigma; stack_shares_pe.py: SFIL2-4 shares
with the photostatistics SFIL2/SFIL3 gain ratio against the model share profile).

    python3 performance_3064.py [glob of data_fit_layers_3064_part*.npz] [simdir] [tag]  -> performance_3064<tag>.txt/.png

Per shower E_vis = sum_L S_L / g_L x exp(-x_pad / lambda_L) with the v4 4 GeV / P6 anchored gains (8.6 / 24.2 / 18.4 / 14.7 ADC
per MeV of South light, sfil_gain_calibration_v4.txt) and the measured attenuation (SFIL2 1329, SFIL3 1646, else 1500 mm), i.e. the
same recipe as performance.py with the v4 gains; the SFIL5 and Baby BCAL terms wait for the KCU1 rebuild and the CAEN join.
The 14 parts are slices of ~3.85 M triggers, 0.65 h each at the beam rate; the TAC period (3-6 Hz after 9.08 h, ~50k triggers) is only a tail
inside part 13 and is not separated here (the cache has no time stamp); the two beam-off gaps (1.5-2.1 h, 5.8-6.3 h) fall inside parts.
Run 3062 for comparison (performance.txt, v3 gains): core sigma/mu of SFIL1-4 alone P4 0.299, P3 0.265, P6 0.204, P5 0.188.
"""
import sys, glob, numpy as np, uproot, ROOT
sys.path.insert(0, '/sphenix/user/tmengel/eic/bic_testbeam/common'); from plotstyle import setup, finish, legend, headroom
setup(); ROOT.gROOT.SetBatch(True)
PAT = sys.argv[1] if len(sys.argv) > 1 else 'data_fit_layers_3064_part*.npz'
SIMDIR = sys.argv[2] if len(sys.argv) > 2 else '/sphenix/tg/tg01/jets/tmengel/eic/sim/v4'
TAG = sys.argv[3] if len(sys.argv) > 3 else ''
out = open('performance_3064%s.txt' % TAG, 'w')
def say(*a): s = ' '.join(str(x) for x in a); print(s, flush=True); out.write(s + '\n')
G14 = np.array([8.6, 24.2, 18.4, 14.7])                 # ADC per MeV of South light, v4 model, 4 GeV under P6
LAMB = np.array([1500., 1329., 1646., 1500.]); PXm = {4: -80., 3: -30., 6: 30., 5: 80.}
A_PE = {1: 14.09, 2: 9.47}; gfac = np.array([1., A_PE[2] / A_PE[1], 1., 1.])   # photostatistics gain ratio SFIL2/SFIL3 (stack_shares_pe.py)
THR = 600.                                                 # shower: South total > 600 ADC, as in 3062
def core(v):
    if len(v) < 50: return np.nan, np.nan
    m, s = np.median(v), 1.4826 * np.median(np.abs(v - np.median(v)))
    for _ in range(5):
        sel = np.abs(v - m) < 2 * s
        if sel.sum() < 20: break
        m, s = v[sel].mean(), v[sel].std()
    return m, s
# ---------------------------------------------------------------- data
files = sorted(glob.glob(PAT)); S, pad, part = [], [], []
for i, f in enumerate(files):
    D = np.load(f); S.append(D['S'].astype(float)); pad.append(D['pad']); part.append(np.full(len(D['pad']), i))
    say('%s: %d triggers, %d showers' % (f.split('/')[-1], len(D['pad']), (D['S'].sum(axis=1) > THR).sum()))
S = np.concatenate(S); pad = np.concatenate(pad); part = np.concatenate(part); tot = S.sum(axis=1); shower = tot > THR
NPART = len(files); TACPART = -1      # no part is a pure TAC sample (see the docstring)
xp = np.array([PXm.get(int(p), 0.) for p in pad]); att = np.exp(-xp[:, None] / LAMB[None, :])
E14 = S / G14[None, :] * att; Evis = E14.sum(axis=1)
beam = shower & (part != TACPART) if TACPART >= 0 else shower
say('\n%d parts, %d triggers, %d showers (South total > %.0f ADC), %d in the beam parts' % (NPART, len(pad), shower.sum(), THR, beam.sum()))
# ---------------------------------------------------------------- model (SFIL1-4 South light, cells 1-3 as the data; v4 scan)
sim = {}
for fn in sorted(glob.glob(SIMDIR + '/bic_sim_e*_s*.root')):
    if '_x' in fn: continue
    g = uproot.open(fn)['sim'].arrays(['ebeam', 'sfil_ls'], library='np'); E = int(round(g['ebeam'].mean()))
    if 2500 <= E <= 6500: sim[E] = g['sfil_ls'][:, :4, :3].sum(axis=2)
Es = sorted(sim); msh = {E: sim[E].mean(axis=0) / sim[E].mean(axis=0).sum() for E in Es}
mres = {}
for E in Es:
    m, s = core(sim[E].sum(axis=1)); mres[E] = (m, s / m)
say('model (%s) SFIL1-4 South light: ' % SIMDIR + '  '.join('%d GeV <L> %.1f MeV sigma/mu %.4f' % (E // 1000, mres[E][0], mres[E][1]) for E in Es))
def share24(E):
    sh = np.array([np.interp(E, Es, [msh[e][k] for e in Es]) for k in (1, 2, 3)]); return sh / sh.sum()
Egrid = np.arange(3000, 6001, 25)
def shape_energy(d):
    d = d / d.sum(); chi = np.array([(((d - share24(E)) / (0.02 + 0.0 * share24(E))) ** 2).sum() for E in Egrid])
    j = chi.argmin(); ok = chi < chi[j] + 1; return Egrid[j], Egrid[ok].min(), Egrid[ok].max(), chi[j]
# ---------------------------------------------------------------- per paddle, whole beam period
say('\nper paddle (beam parts): n, SFIL1-4 South shares (equal gains), pe-equalised SFIL2-4 shares -> shape energy, <E_vis> [MeV] and core sigma/mu with the v4 4 GeV gains; run 3062 SFIL1-4 alone: P4 0.299 P3 0.265 P6 0.204 P5 0.188')
res = {}
for p in (4, 3, 6, 5):
    m = beam & (pad == p); sh = S[m].mean(axis=0); sh = sh / sh.sum(); d = (S[m].mean(axis=0) * gfac)[1:]; Ef, lo, hi, c2 = shape_energy(d)
    mu, sg = core(Evis[m]); res[p] = (m.sum(), sh, Ef, lo, hi, mu, sg / mu)
    say('  P%d: n %7d  shares %s  pe-eq SFIL2-4 %s -> E %.2f GeV (%.2f-%.2f, chi2 %.1f)  <E_vis> %.1f MeV  sigma/mu %.4f  (model at that E: %.4f)' % (
        p, m.sum(), sh.round(3), (d / d.sum()).round(3), Ef / 1000, lo / 1000, hi / 1000, c2, mu, sg / mu, np.interp(Ef, Es, [mres[e][1] for e in Es])))
m = beam; d = (S[m].mean(axis=0) * gfac)[1:]; Ef, lo, hi, c2 = shape_energy(d); mu, sg = core(Evis[m])
say('  all: n %7d  pe-eq SFIL2-4 %s -> E %.2f GeV (%.2f-%.2f)  <E_vis> %.1f MeV  sigma/mu %.4f' % (m.sum(), (d / d.sum()).round(3), Ef / 1000, lo / 1000, hi / 1000, mu, sg / mu))
say('linearity (E_vis per GeV of the shape energy): ' + '  '.join('P%d %.1f MeV/GeV' % (p, res[p][5] / (res[p][2] / 1000)) for p in (4, 3, 6, 5)) + '  (model %.1f MeV of SFIL1-4 South light per GeV at 4 GeV)' % (mres[4000][0] / 4 if 4000 in mres else float('nan')))
if TACPART >= 0:
    m = shower & (part == TACPART)
    say('\nTAC period (part %02d): %d showers; per paddle n %s; <E_vis> %.1f MeV, core sigma/mu %.4f; pe-eq SFIL2-4 shares %s' % (
        TACPART, m.sum(), [int((m & (pad == p)).sum()) for p in (4, 3, 6, 5)], core(Evis[m])[0], core(Evis[m])[1] / core(Evis[m])[0] if m.sum() > 50 else float('nan'), ((S[m].mean(axis=0) * gfac)[1:] / (S[m].mean(axis=0) * gfac)[1:].sum()).round(3) if m.sum() else None))
# ---------------------------------------------------------------- per part
say('\nper part (~%d triggers each): per paddle <E_vis> [MeV] and core sigma/mu' % (len(pad) // max(NPART, 1)))
say('%-5s' % 'part' + ''.join('  P%d: <E>  s/mu ' % p for p in (4, 3, 6, 5)) + '   n showers')
PP = {p: [] for p in (4, 3, 6, 5)}
for i in range(NPART):
    line = '%-5s' % ('%02d%s' % (i, 'T' if i == TACPART else ' '))
    for p in (4, 3, 6, 5):
        m = shower & (part == i) & (pad == p); mu, sg = core(Evis[m]); PP[p].append((mu, sg / mu if mu == mu and mu > 0 else np.nan, m.sum())); line += '  %6.1f %6.4f' % (mu, sg / mu if mu == mu and mu > 0 else np.nan)
    say(line + '   %7d' % (shower & (part == i)).sum())
for p in (4, 3, 6, 5):
    a = np.array([v for i, v in enumerate(PP[p]) if i != TACPART and v[2] > 500])
    if len(a): say('  P%d over the beam parts: <E_vis> mean %.1f rms %.1f (%.2f %%), sigma/mu mean %.4f rms %.4f' % (p, a[:, 0].mean(), a[:, 0].std(), 100 * a[:, 0].std() / a[:, 0].mean(), a[:, 1].mean(), a[:, 1].std()))
# ---------------------------------------------------------------- figure
keep = []; cv = ROOT.TCanvas('c', '', 2100, 700); cv.Divide(3, 1, 0.002, 0.002); keep.append(cv)
COL = {4: ROOT.kBlue + 1, 3: ROOT.kGreen + 2, 6: ROOT.kOrange + 7, 5: ROOT.kRed + 1}
cv.cd(1); hs = []; lg = legend(0.50, 0.62, 0.92, 0.88); keep.append(lg)
for p in (4, 3, 6, 5):
    m = beam & (pad == p); h = ROOT.TH1D('h%d' % p, 'run 3064, beam period, SFIL1-4 South;E_{vis} SFIL1-4 [MeV, v4 4 GeV gains];fraction per 5 MeV', 100, 0, 500); keep.append(h)
    v = Evis[m]; h.FillN(len(v), np.ascontiguousarray(v, float), np.ones(len(v)))
    if h.GetEntries(): h.Scale(1. / h.GetEntries())
    h.SetLineColor(COL[p]); h.SetLineWidth(2); hs.append(h); lg.AddEntry(h, 'P%d: #sigma/#mu %.1f %%' % (p, 100 * res[p][6]), 'l')
mx = max(h.GetMaximum() for h in hs); headroom(hs[0], 0, mx, 0.35); hs[0].Draw('hist')
for h in hs[1:]: h.Draw('hist same')
lg.Draw()
cv.cd(2); fr = ROOT.TH2D('fr', 'core #sigma/#mu per part;part (0.65 h each);#sigma/#mu', 10, -0.5, NPART - 0.5, 10, 0.15, 0.35); fr.SetStats(0); keep.append(fr); fr.Draw()
lg2 = legend(0.62, 0.66, 0.92, 0.88); keep.append(lg2)
for p in (4, 3, 6, 5):
    a = np.array(PP[p]); x = np.arange(NPART, dtype=float); ok = np.isfinite(a[:, 1]) & (a[:, 2] > 200)
    g = ROOT.TGraph(int(ok.sum()), np.ascontiguousarray(x[ok]), np.ascontiguousarray(a[ok, 1])); g.SetMarkerStyle(20); g.SetMarkerColor(COL[p]); g.SetLineColor(COL[p]); g.SetMarkerSize(1.1); g.Draw('P'); keep.append(g); lg2.AddEntry(g, 'P%d' % p, 'p')
    if TACPART >= 0 and ok[TACPART]:
        gt = ROOT.TGraph(1, np.array([float(TACPART)]), np.array([a[TACPART, 1]])); gt.SetMarkerStyle(24); gt.SetMarkerColor(COL[p]); gt.SetMarkerSize(1.4); gt.Draw('P'); keep.append(gt)
lg2.Draw()
cv.cd(3); fr3 = ROOT.TH2D('fr3', '<E_{vis}> per part;part (0.65 h each);<E_{vis}> SFIL1-4 [MeV]', 10, -0.5, NPART - 0.5, 10, 105, 150); fr3.SetStats(0); keep.append(fr3); fr3.Draw()
lg3 = legend(0.62, 0.66, 0.92, 0.88); keep.append(lg3)
for p in (4, 3, 6, 5):
    a = np.array(PP[p]); x = np.arange(NPART, dtype=float); ok = np.isfinite(a[:, 0]) & (a[:, 2] > 200)
    g = ROOT.TGraph(int(ok.sum()), np.ascontiguousarray(x[ok]), np.ascontiguousarray(a[ok, 0])); g.SetMarkerStyle(20); g.SetMarkerColor(COL[p]); g.SetLineColor(COL[p]); g.SetMarkerSize(1.1); g.Draw('P'); keep.append(g); lg3.AddEntry(g, 'P%d' % p, 'p')
lg3.Draw(); finish(cv); cv.SaveAs('performance_3064%s.png' % TAG); say('wrote performance_3064%s.png' % TAG); out.close()
