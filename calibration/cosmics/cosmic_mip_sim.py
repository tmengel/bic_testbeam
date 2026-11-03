"""ADC per MeV from the cosmic run 3088: the two-end muon tag applied to the data and to the cosmic-mode simulation.

    python3 cosmic_mip_sim.py [data cosmic.root] [sim glob] [tag] [pad]  -> cosmic_mip_sim_<tag>.txt, cosmic_mip_sim_<tag>.png  ('pad' = simulated muons with a paddle hit only)

Data (cosmic_mip.C tree): per cell the events with both ends above THR ADC, peaking at samples 7-9 within one sample of each other
and with the pulse shape at one end at least; spectrum of sqrt(A_S A_N). Simulation (cosmic mode of bic_sim, geometry v4): per cell
the crossings with light at both ends; the readout model per end is  A = Poisson(L x Y) x G + N(0, sigma)  with L the simulated light
at that end (MeV of South/North light, attenuation included), Y the light yield in photoelectrons per MeV of end light (the one free
parameter per cell), G the ADC per photoelectron (the run-3062 photostatistics values where measured, else 10; North = South x the
measured N/S gain ratio of the layer, 0.85 / 1.03 / 0.29 for SFIL2-4) and sigma the channel noise of the cosmic data after the
common-mode subtraction. The same tag (both ends > THR) is applied; with 'pad' only the simulated muons that hit a paddle are used
(all data tags carry a paddle pulse; the sim/cosmic/tgt production sends every muon through a paddle). Y is fitted by the binned
likelihood of the data spectrum against the simulated template, and the gain follows as ADC per MeV of end light = Y x G, directly
comparable with the beam-anchored gains (ADC per MeV of simulated South light: 8.6 / 24.2 / 18.4 / 14.7 for SFIL1-4 at 4 GeV, v4).
The fit is repeated for G x 0.7 and G x 1.4 (the excess-noise uncertainty of the ADC/pe values).

Cells 1 and 2 (2026-10-30): on the v2 SiPM boards (all South, North SFIL2-4) the channels "cell 1" and "cell 2" both read SiPM1 + SiPM2
mixed (r = 0.97-0.99 between them; hookup note), so a muon through physical cell 1 or 2 lights both channels and each per-channel fit
sees only a share of that light. Besides the per-channel rows, the element "c1+2" uses the SUM of the two channels per end (tag: either
channel passes the timing and shape, the sums above THR x sqrt 2) against the simulated light of cells 1 + 2 (noise added in quadrature).
If the passive network conserves the charge, the c1+2 gain is the board gain and should agree with cells 3 and 4.
"""
import sys, glob, numpy as np, uproot, ROOT
sys.path.insert(0, '/sphenix/user/tmengel/eic/bic_testbeam/common'); from plotstyle import setup, finish, legend, headroom
setup()
F = sys.argv[1] if len(sys.argv) > 1 else '/sphenix/tg/tg01/jets/tmengel/eic/kcu0/cosmics/cosmic_3088.root'
SIM = sys.argv[2] if len(sys.argv) > 2 else '/sphenix/tg/tg01/jets/tmengel/eic/sim/cosmic/bic_cosmic_s*.root'
TAG = sys.argv[3] if len(sys.argv) > 3 else '3088'
THR = 8.; TAIL = 0.8; GBEAM = [8.6, 24.2, 18.4, 14.7]
GPE = {(1, 0): 14.1, (2, 0): 9.5}          # ADC per photoelectron (x F) measured on the South beam cells of SFIL2 / SFIL3; default 10 elsewhere
NSR = {1: 0.85, 2: 1.03, 3: 0.29}          # North / South gain ratio of the beam cell from the AstroPix-x attenuation study (SFIL2 / 3 / 4); 1 for SFIL1
PADSIM = (sys.argv[4] == 'pad') if len(sys.argv) > 4 else False   # 'pad' as 4th argument: require a paddle hit (> 0.5 MeV) in the simulated muon
out = open('cosmic_mip_sim_%s.txt' % TAG, 'w')
def say(*a): s = ' '.join(str(x) for x in a); print(s); out.write(s + '\n')
rng = np.random.default_rng(1)
# ---------------------------------------------------------------- data
t = uproot.open(F)['cos']; n = t.num_entries
A = np.zeros((n, 36), np.float32); PK = np.zeros((n, 36), np.int8); shape = np.zeros((n, 36), bool); NZ = np.zeros(36); nz = 0; i0 = 0
import awkward as ak
for arr in t.iterate(['A', 'pk', 'nwf', 'wf', 'noise'], step_size=40000, library='ak'):
    m = len(arr); a = ak.to_numpy(arr['A'])[:, :36]; pkk = ak.to_numpy(arr['pk'])[:, :36].astype(int); nw = ak.to_numpy(arr['nwf']); nzc = ak.to_numpy(arr['noise'])[:, :36]
    A[i0:i0 + m] = a; PK[i0:i0 + m] = pkk; ok = nzc > -1; NZ += np.where(ok, nzc, 0).sum(axis=0); nz += ok.sum(axis=0)
    has = nw > 0
    if has.any():
        wf = ak.to_numpy(arr['wf'][has]).reshape(int(has.sum()), 36, 20) / 10.
        pc = np.clip(pkk[has], 5, 11); tail = np.take_along_axis(wf, np.clip(pc + 4, 0, 19)[:, :, None], axis=2)[:, :, 0]; rise = np.take_along_axis(wf, np.clip(pc - 2, 0, 19)[:, :, None], axis=2)[:, :, 0]
        shape[i0:i0 + m][has] = (tail < TAIL * a[has]) & (rise < 0.5 * a[has])
    i0 += m
noise = NZ / np.maximum(nz, 1)
say('data %s: %d triggers; noise per channel after the common-mode subtraction [ADC]: %s' % (F, n, np.round(noise, 1).tolist()))
good = lambda k: (A[:, k] > THR) & (PK[:, k] >= 7) & (PK[:, k] <= 9)
# ---------------------------------------------------------------- simulation
files = sorted(glob.glob(SIM)); LS, LN, PE = [], [], []
for f in files:
    a = uproot.open(f)['sim'].arrays(['sfil_ls', 'sfil_ln', 'paddle_e'], library='np'); LS.append(a['sfil_ls'][:, :4, :4]); LN.append(a['sfil_ln'][:, :4, :4]); PE.append(a['paddle_e'])
LS = np.concatenate(LS); LN = np.concatenate(LN); PE = np.concatenate(PE); padhit = (PE > 0.5).any(axis=1)
say('simulation: %d files, %d muons, %d with a paddle hit (> 0.5 MeV); paddle-tagged subsample used: %s' % (len(files), len(LS), padhit.sum(), PADSIM))
# ---------------------------------------------------------------- per element
edges = np.arange(0, 121, 4.)
def template(ls, ln, Y, gS, gN, sS, sN, nrep, thr):
    """simulated sqrt(A_S A_N) after the tag, nrep readout realisations per crossing"""
    ls = np.repeat(ls, nrep); ln = np.repeat(ln, nrep); r = np.random.default_rng(12345)   # common random numbers: the same stream for every Y and G, so the likelihood scan is smooth
    aS = r.poisson(ls * Y) * gS + r.normal(0, sS, len(ls)); aN = r.poisson(ln * Y) * gN + r.normal(0, sN, len(ln))
    m = (aS > thr) & (aN > thr); return np.sqrt(aS[m] * aN[m]), m.mean()
def nll(h, tmpl_counts):
    mu = tmpl_counts * h.sum() / max(tmpl_counts.sum(), 1e-9) + 1e-9; return (mu - h * np.log(mu)).sum()
keep = []; cv = ROOT.TCanvas('c', '', 2500, 1500); cv.Divide(5, 4, 0.002, 0.002); keep.append(cv)
say('\nper element: data tags, simulated crossings; fitted light yield Y [pe per MeV of end light] and ADC per MeV of end light Y x G (G = ADC/pe used: South, North),')
say('for G, 0.7 G and 1.4 G; beam-anchored gain for the layer (ADC per MeV of South light, 4 GeV anchor, v4); tag efficiency of the simulated crossings')
say('c1+2 = channel sums of cells 1 and 2 (one SiPM pair read by both channels) against the light of physical cells 1 + 2')
res = {}
def element(L, c, label, m, gm, ls_all, ln_all, sS, sN, thr, ipad):
    cross = (ls_all > 0.05) & (ln_all > 0.05); ncross_all = cross.sum()
    if PADSIM: cross &= padhit
    ls, ln = ls_all[cross], ln_all[cross]; h, _ = np.histogram(gm, edges)
    gS = GPE.get((L, 0), 10.); gN = gS * NSR.get(L, 1.0)
    line = 'SFIL%d %-4s: data %5d, sim crossings %5d (%5d with a paddle, %5d used; light S/N median %.2f / %.2f MeV)' % (L + 1, label, m.sum(), ncross_all, ((ls_all > 0.05) & (ln_all > 0.05) & padhit).sum(), cross.sum(), np.median(ls) if cross.any() else 0, np.median(ln) if cross.any() else 0)
    if m.sum() < 30 or cross.sum() < 50: say(line + '  - too few'); return
    best = {}
    for fac in (1.0, 0.7, 1.4):
        Ys = np.exp(np.linspace(np.log(0.3), np.log(30), 94)); nl = []
        for Y in Ys:
            tm, eff = template(ls, ln, Y, fac * gS, fac * gN, sS, sN, max(2, int(80000 / cross.sum())), thr)
            th, _ = np.histogram(tm, edges); nl.append(nll(h, th + 0.5))
        nl = np.array(nl); j = nl.argmin(); ok = (nl < nl[j] + 0.5); lo, hi = Ys[ok].min() / 1.025, Ys[ok].max() * 1.025   # half a grid step added on each side
        tm, eff = template(ls, ln, Ys[j], fac * gS, fac * gN, sS, sN, max(2, int(100000 / cross.sum())), thr)
        best[fac] = (Ys[j], lo, hi, Ys[j] * fac * gS, eff, tm)
    b = best[1.0]
    line += ' | Y %.2f [%.2f, %.2f] pe/MeV -> %.1f ADC/MeV (G %.1f / %.1f); 0.7G: %.1f, 1.4G: %.1f; beam %.1f; tag eff %.2f' % (b[0], b[1], b[2], b[3], gS, gN, best[0.7][3], best[1.4][3], GBEAM[L], b[4])
    say(line); res[(L, label)] = best
    cv.cd(ipad)
    hd = ROOT.TH1D('hd_%d_%s' % (L, label), 'SFIL%d %s: data (%d) and cosmic model at Y = %.1f pe/MeV;#sqrt{A_{S}A_{N}} [ADC];tags' % (L + 1, label, m.sum(), b[0]), len(edges) - 1, edges); keep.append(hd)
    for v in gm: hd.Fill(float(v))
    hm = ROOT.TH1D('hm_%d_%s' % (L, label), '', len(edges) - 1, edges); keep.append(hm)
    for v in b[5]: hm.Fill(float(v))
    if hm.Integral() > 0: hm.Scale(hd.Integral() / hm.Integral())
    hd.SetLineWidth(2); hm.SetLineColor(ROOT.kRed + 1); hm.SetLineWidth(2); headroom(hd, 0, max(hd.GetMaximum(), hm.GetMaximum()), 0.3); hd.Draw('hist'); hm.Draw('hist same')
    lg = legend(0.45, 0.70, 0.92, 0.86); keep.append(lg); lg.AddEntry(hd, 'data', 'l'); lg.AddEntry(hm, 'model %.1f ADC/MeV' % b[3], 'l'); lg.Draw()
for L in range(4):
    for c in range(4):
        s, nn = 8 * L + c, 8 * L + 4 + c
        m = good(s) & good(nn) & (np.abs(PK[:, s].astype(int) - PK[:, nn].astype(int)) <= 1) & (shape[:, s] | shape[:, nn])
        element(L, c, 'c%d' % (c + 1), m, np.sqrt(A[m, s] * A[m, nn]), LS[:, L, c], LN[:, L, c], max(noise[s], 0.8), max(noise[nn], 0.8), THR, 5 * L + c + 1)
    # cells 1 + 2 as one element: channel sums per end
    s1, s2, n1, n2 = 8 * L, 8 * L + 1, 8 * L + 4, 8 * L + 5
    AS = A[:, s1] + A[:, s2]; AN = A[:, n1] + A[:, n2]; thr2 = THR * np.sqrt(2)
    pkS = np.where(A[:, s1] >= A[:, s2], PK[:, s1], PK[:, s2]).astype(int); pkN = np.where(A[:, n1] >= A[:, n2], PK[:, n1], PK[:, n2]).astype(int)
    m = (good(s1) | good(s2)) & (good(n1) | good(n2)) & (AS > thr2) & (AN > thr2) & (np.abs(pkS - pkN) <= 1) & (shape[:, s1] | shape[:, s2] | shape[:, n1] | shape[:, n2])
    element(L, 0, 'c1+2', m, np.sqrt(AS[m] * AN[m]), LS[:, L, 0] + LS[:, L, 1], LN[:, L, 0] + LN[:, L, 1], np.hypot(max(noise[s1], 0.8), max(noise[s2], 0.8)), np.hypot(max(noise[n1], 0.8), max(noise[n2], 0.8)), thr2, 5 * L + 5)
finish(cv); cv.SaveAs('cosmic_mip_sim_%s.png' % TAG)
out.close()
