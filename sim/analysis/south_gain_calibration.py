"""Preliminary ADC -> MeV gains of the South-end channels from the paddle-tagged showers of run 3062 against the
simulation, and the gain-free double ratios that test the P5 / P6 energy assignment.

Data (bic_testbeam/data/events_3062.root): showers seen by all three boards, tagged by paddle P5 (South end of the
module) or P6 (second from South); per event the South-end peak amplitudes of SFIL1-4 (cells 1+2, 3, 4), the SFIL5 South
array sum (KCU1 B+C+D) and the Baby BCAL South charge (CAEN ch 0-11).
Simulation: South-end light (sfil_ls, bb_ls) with the beam under P6 (x = 30 mm, default scan) and under P5 (x = 80 mm).

    python3 south_gain_calibration.py  -> south_gain_calibration.txt, south_gain_calibration.png
"""
import glob, sys, numpy as np, uproot, ROOT
sys.path.insert(0, '/sphenix/user/tmengel/eic/bic_testbeam/common'); from plotstyle import setup
SIM = (sys.argv[1] if len(sys.argv) > 1 else '/sphenix/tg/tg01/jets/tmengel/eic/sim/v1_purePb') + '/'
out = open('south_gain_calibration.txt', 'w')
def say(*a): print(*a); print(*a, file=out)

# ---- data
f = uproot.open('/sphenix/user/tmengel/eic/bic_testbeam/data/events_3062.root')['events']
a = f.arrays(['has_kcu0', 'has_kcu1', 'has_caen', 'kcu0_frag', 'pad', 'sfil_s', 'sfil_n', 'sfil5_s', 'caen_q'], library='np')
S = np.clip(a['sfil_s'], 0, None).astype(float); N = np.clip(a['sfil_n'], 0, None).astype(float)
e4 = np.sqrt(S * N)[:, [0, 2, 3], :].sum(axis=(1, 2))
s5 = np.clip(a['sfil5_s'], 0, None).astype(float); e5 = s5[:, 18:35].sum(axis=1) + s5[:, 36:53].sum(axis=1) + s5[:, 54:71].sum(axis=1)
bb = np.clip(a['caen_q'][:, :12], 0, None).astype(float).sum(axis=1)
base = (a['has_kcu0'] == 1) & (a['kcu0_frag'] == 0) & (a['has_kcu1'] == 1) & (a['has_caen'] == 1) & (e4 > 1200)
ELEM = ['SFIL1 c1+2', 'SFIL1 c3', 'SFIL1 c4', 'SFIL2 c1+2', 'SFIL2 c3', 'SFIL2 c4', 'SFIL3 c1+2', 'SFIL3 c3', 'SFIL3 c4', 'SFIL4 c1+2', 'SFIL4 c3', 'SFIL4 c4', 'SFIL5 S', 'BBCAL S']
def data_vec(m):
    v = []; err = []
    for L in range(4):
        for c in (0, 2, 3): v.append(S[m, c, L].mean()); err.append(S[m, c, L].std() / np.sqrt(m.sum()))
    v += [e5[m].mean(), bb[m].mean()]; err += [e5[m].std() / np.sqrt(m.sum()), bb[m].std() / np.sqrt(m.sum())]
    return np.array(v), np.array(err)
PBIT = {3: 1, 4: 2, 5: 4, 6: 8}   # 'pad' is a bitmask: P3 = 1, P4 = 2, P5 = 4, P6 = 8 (single-paddle events = exactly one bit)
data = {p: data_vec(base & (a['pad'] == PBIT[p])) for p in (5, 6)}; nd = {p: int((base & (a['pad'] == PBIT[p])).sum()) for p in (5, 6)}
say('data: P5 %d showers, P6 %d showers (triple-matched, SFIL1-4 sum > 1200)' % (nd[5], nd[6]))

# ---- simulation: South-end light per element, per beam position and energy (MeV)
def sim_vec(path):
    g = uproot.open(path)['sim'].arrays(['sfil_ls', 'bb_ls'], library='np'); ls = g['sfil_ls']
    v = []
    for L in range(4): v += [(ls[:, L, 0] + ls[:, L, 1]).mean(), ls[:, L, 2].mean(), ls[:, L, 3].mean()]   # cells 1+2 (one row on v2 boards), 3, 4
    v += [ls[:, 4, :].sum(axis=1).mean(), g['bb_ls'].sum(axis=(1, 2)).mean()]
    return np.array(v)
sim = {}   # (paddle, E) -> vector
for fpath in sorted(glob.glob(SIM + 'bic_sim_e*_s*.root')):
    tag = fpath.split('/')[-1]; E = int(tag.split('_e')[1].split('_')[0]); pos = 5 if '_x80' in tag else 6
    if E < 2500: continue
    try: sim[(pos, E)] = sim_vec(fpath)
    except Exception as ex: say('skipping', tag, '(not readable yet)'); continue
say('simulation points available:', sorted(sim))

# ---- gains g = data / sim (ADC or charge per MeV of South light) under energy hypotheses
hyp = [('scaled May map: P6 5 GeV, P5 6 GeV', 5000, 6000), ('equal energies 5 GeV', 5000, 5000), ('equal energies 6 GeV', 6000, 6000), ('reversed: P6 6 GeV, P5 5 GeV', 6000, 5000)]
say('\nPer-element gains data/sim, and the P5/P6 gain ratio (1 = consistent) per energy hypothesis')
res = {}
for name, E6, E5 in hyp:
    if (6, E6) not in sim or (5, E5) not in sim: say('  %s: simulation missing' % name); continue
    g6 = data[6][0] / sim[(6, E6)]; g5 = data[5][0] / sim[(5, E5)]; r = g5 / g6; res[name] = (g6, g5, r)
    say('\n  %s' % name); say('  %-12s %12s %12s %10s' % ('element', 'gain P6', 'gain P5', 'P5/P6'))
    for i, el in enumerate(ELEM): say('  %-12s %12.3f %12.3f %10.3f' % (el, g6[i], g5[i], r[i]))
    core = [1, 4, 7, 10, 12, 13]   # beam-row cells, SFIL5, BBCAL
    say('  beam-row + SFIL5 + BBCAL: mean P5/P6 ratio %.3f, rms %.3f' % (np.mean(r[core]), np.std(r[core])))
# ---- gain-free double ratios: (element / SFIL3 c3) in P5 over the same in P6, data vs simulation
say('\nGain-free double ratios [element/SFIL3 c3]_P5 / [element/SFIL3 c3]_P6: data, and simulation per hypothesis')
ref = 7
dd = (data[5][0] / data[5][0][ref]) / (data[6][0] / data[6][0][ref])
say('  %-12s %8s' % ('element', 'data') + ''.join(' %22s' % n[:22] for n, _, _ in hyp))
for i, el in enumerate(ELEM):
    line = '  %-12s %8.3f' % (el, dd[i])
    for name, E6, E5 in hyp:
        if (6, E6) in sim and (5, E5) in sim: ss = (sim[(5, E5)] / sim[(5, E5)][ref]) / (sim[(6, E6)] / sim[(6, E6)][ref]); line += ' %22.3f' % ss[i]
        else: line += ' %22s' % '-'
    say(line)
# total South response ratio P5/P6 in data (SFIL1-4 beam row + off-row) vs the simulated energy ratio
d_tot = {p: data[p][0][:12].sum() for p in (5, 6)}
say('\nSFIL1-4 South total (all cells): data P5/P6 = %.3f;  simulation South light P5(x=80)/P6(x=30) at equal energy: %s' % (d_tot[5] / d_tot[6], ', '.join('%d GeV %.3f' % (E, sim[(5, E)][:12].sum() / sim[(6, E)][:12].sum()) for E in (3000, 4000, 5000, 6000) if (5, E) in sim and (6, E) in sim)))
out.close()

# ---- figure: gains per element for the scaled map hypothesis, P5 and P6
name = hyp[0][0]
if name in res:
    g6, g5, r = res[name]; setup(); ROOT.gStyle.SetPadBottomMargin(0.22); ROOT.gStyle.SetPadLeftMargin(0.12); ROOT.gROOT.ForceStyle()
    c = ROOT.TCanvas('c', '', 1100, 700); c.Divide(1, 2, 0, 0); keep = []
    c.cd(1); ROOT.gPad.SetLogy(); fr = ROOT.gPad.DrawFrame(-0.5, 0.5, 13.5, 2000, ';;data / sim  [ADC per MeV]'); keep.append(fr)
    x = np.arange(14, dtype=float); ga = ROOT.TGraph(14, x - 0.1, g6); ga.SetMarkerStyle(21); ga.SetMarkerColor(4); ga.Draw('P'); gb = ROOT.TGraph(14, x + 0.1, g5); gb.SetMarkerStyle(20); gb.SetMarkerColor(2); gb.Draw('P'); keep += [ga, gb]
    lg = ROOT.TLegend(0.15, 0.7, 0.55, 0.92); lg.SetBorderSize(0); lg.SetFillStyle(0); lg.AddEntry(ga, 'P6 showers vs sim 5 GeV under P6', 'p'); lg.AddEntry(gb, 'P5 showers vs sim 6 GeV under P5', 'p'); lg.Draw(); keep.append(lg)
    c.cd(2); fr2 = ROOT.gPad.DrawFrame(-0.5, 0.4, 13.5, 1.6, ';;gain P5 / gain P6'); keep.append(fr2)
    for i, el in enumerate(ELEM): fr2.GetXaxis().SetBinLabel(fr2.GetXaxis().FindBin(i), el)
    fr2.GetXaxis().SetLabelSize(0.07); fr2.GetXaxis().LabelsOption('v'); gr = ROOT.TGraph(14, x, r); gr.SetMarkerStyle(20); gr.Draw('P'); keep.append(gr); l1 = ROOT.TLine(-0.5, 1, 13.5, 1); l1.SetLineStyle(2); l1.Draw(); keep.append(l1)
    c.SaveAs('south_gain_calibration.png'); print('wrote south_gain_calibration.png')
