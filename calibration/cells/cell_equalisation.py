"""Cell-to-cell equalisation within a layer.  The beam is fixed in y (cell 3 of 4), so cells 1, 2 and 4 see only the lateral
tails of the shower; their mean amplitude relative to the beam cell, divided by the model's lateral light fraction, is the
relative gain g_c / g_3 per layer and end.  The model fraction depends on the beam height, so every available beam-y run
(tag _yb<mm>, nominal -5 mm) is shown.  Data: events table raw window amplitudes, both cells unsaturated, showers > 600 ADC.
    python3 cell_equalisation.py [simdir] [kcu0 file] [tag] -> cell_equalisation<tag>.txt/.png   (default: run 3062, no tag)
"""
import numpy as np, uproot, glob, re, sys, ROOT
sys.path.insert(0, '/sphenix/user/tmengel/eic/bic_testbeam/common'); from plotstyle import setup, legend, label
TAG = sys.argv[3] if len(sys.argv) > 3 else ''
setup(); ROOT.gROOT.SetBatch(True); out = open('cell_equalisation%s.txt' % TAG, 'w')
def say(*a): print(*a, flush=True); print(*a, file=out)
SIMDIR = sys.argv[1] if len(sys.argv) > 1 else '/sphenix/tg/tg01/jets/tmengel/eic/sim'
NMAX = 800000
f = uproot.open(sys.argv[2] if len(sys.argv) > 2 else '/sphenix/user/tmengel/eic/bic_testbeam/data/kcu0_3062.root')
def peaks(tree, names):
    a = f[tree].arrays([n + '_adc' for n in names], entry_stop=NMAX, library='np'); P = np.full((NMAX, len(names)), np.nan)
    for j, n in enumerate(names):
        W = a[n + '_adc'].astype(float); good = W[:, 0] >= 0; ped = W[:, 1:4].mean(axis=1); pk = (W[:, 4:13] - ped[:, None]).max(axis=1); P[good, j] = pk[good]
    return P
Sraw = np.stack([peaks('sfil%d' % L, ['Sc1', 'Sc2', 'Sc3', 'Sc4']) for L in (1, 2, 3, 4)], axis=1)   # [n, layer, cell], cells by branch name
Nraw = np.stack([peaks('sfil%d' % L, ['Nc1', 'Nc2', 'Nc3', 'Nc4']) for L in (1, 2, 3, 4)], axis=1)
tot = np.nan_to_num(Sraw).sum(axis=(1, 2)); base = np.isfinite(Sraw).all(axis=(1, 2)); shower = base & (tot > 600)
say('first %d KCU0 triggers from the raw trees (cells by branch name Sc1-Sc4 / Nc1-Nc4): %d showers' % (NMAX, shower.sum()))
boards = {('S', 1): 'v2', ('S', 2): 'v2', ('S', 3): 'v2', ('S', 4): 'v2', ('N', 1): 'v1', ('N', 2): 'v2', ('N', 3): 'v2', ('N', 4): 'v2'}
data = {}
empty = base & (tot < 100)
say('noise bias of the raw peak in empty triggers (ADC), South / North per layer and cell: %s / %s' % (np.round(np.nanmean(Sraw[empty], axis=0), 1).tolist(), np.round(np.nanmean(Nraw[empty], axis=0), 1).tolist()))
for end, R in (('S', Sraw), ('N', Nraw)):
    bias = np.nanmean(R[empty], axis=0)
    for L in range(4):
        beam = shower & (R[:, L, 2] < 900) & (R[:, L, 2] > 100)
        row = {}
        for c in (0, 1, 3):
            m = beam & (R[:, L, c] < 900); num = R[m, L, c].mean() - bias[L, c]; den = R[m, L, 2].mean() - bias[L, 2]
            row[c] = (num / den, np.median((R[m, L, c] - bias[L, c]) / np.maximum(R[m, L, 2] - bias[L, 2], 1)), m.sum())
        # the v2 pair: cells 1 + 2 together against cell 3
        m = beam & (R[:, L, 0] < 900) & (R[:, L, 1] < 900); row['12'] = ((R[m, L, 0] + R[m, L, 1]).mean() - bias[L, 0] - bias[L, 1]) / (R[m, L, 2].mean() - bias[L, 2])
        data[(end, L)] = row
        say('  SFIL%d %s (%s board): (<A_c> - bias)/(<A_3> - bias) = cell1 %.4f  cell2 %.4f  cell4 %.4f  cells1+2 %.4f   (medians %.4f %.4f %.4f; n %d)' % (L + 1, end, boards[(end, L + 1)], row[0][0], row[1][0], row[3][0], row['12'], row[0][1], row[1][1], row[3][1], row[0][2]))
# model: per beam-y run, light fraction per cell / cell 3 per layer and end (4 GeV)
sims = {}
for fn in sorted(glob.glob(SIMDIR + '/bic_sim_e4000*_s*.root')):
    if '_x80' in fn or '_bby' in fn: continue
    m = re.search(r'_yb(-?[\d.]+)', fn); yb = float(m.group(1)) if m else -5.
    try: g = uproot.open(fn)['sim'].arrays(['sfil_ls', 'sfil_ln'], library='np')
    except Exception as e: say('  skipping %s (%s)' % (fn, e)); continue
    sims[yb] = (g['sfil_ls'][:, :4, :4].mean(axis=0), g['sfil_ln'][:, :4, :4].mean(axis=0))
say('\nmodel (4 GeV) light fraction per cell over cell 3, by beam y (mm):')
for yb in sorted(sims):
    LS, LN = sims[yb]
    say('  y = %+5.1f: South ' % yb + ' | '.join('SFIL%d %.3f %.3f %.3f' % (L + 1, LS[L, 0] / LS[L, 2], LS[L, 1] / LS[L, 2], LS[L, 3] / LS[L, 2]) for L in range(4)))
say('\nimplied relative gains g_c / g_3 = data ratio / model ratio (model at the nominal beam y = -5 mm; the spread over the beam-y runs is the lateral-position systematic):')
ybs = sorted(sims); nom = -5. if -5. in sims else ybs[0]
res = {}
for end in 'SN':
    for L in range(4):
        row = data[(end, L)]; vals = []
        for c in (0, 1, 3):
            g_nom = row[c][0] / (sims[nom][0 if end == 'S' else 1][L, c] / sims[nom][0 if end == 'S' else 1][L, 2])
            span = [row[c][0] / (sims[y][0 if end == 'S' else 1][L, c] / sims[y][0 if end == 'S' else 1][L, 2]) for y in ybs]
            vals.append((g_nom, min(span), max(span))); res[(end, L, c)] = vals[-1]
        LSN = sims[nom][0 if end == 'S' else 1]; m12 = (LSN[L, 0] + LSN[L, 1]) / LSN[L, 2]
        say('  SFIL%d %s: cell1 %.2f [%.2f-%.2f]  cell2 %.2f [%.2f-%.2f]  cell4 %.3f [%.3f-%.3f]   cells 1+2 pair: data %.3f / model %.3f = %.2f' % (L + 1, end, *vals[0], *vals[1], *vals[2], row['12'], m12, row['12'] / m12))
say('\nreading: cell 4 is empty in the data (ratio ~0.01 against a model 0.05-0.11), the v2-board South cells 1 and 2 of one layer are a summed pair on the summing boards (Bobae deck p18), so their individual gains are not separable; use the cell 1 + cell 2 sum where needed.')
# figure
c = ROOT.TCanvas('c', '', 2000, 800); c.Divide(2, 1); keep = []
for ip, end in enumerate('SN'):
    c.cd(ip + 1); fr = ROOT.gPad.DrawFrame(0.5, 0.0, 4.5, 0.4, '%s end;layer;<A_c>/<A_3>' % ('South' if end == 'S' else 'North')); keep.append(fr); lg = legend(0.15, 0.6, 0.6, 0.88)
    for c_, mk, nm in ((0, 20, 'cell 1'), (1, 21, 'cell 2'), (3, 22, 'cell 4')):
        gd = ROOT.TGraph(4, np.arange(1, 5, dtype=float), np.array([data[(end, L)][c_][0] for L in range(4)])); gd.SetMarkerStyle(mk); gd.SetMarkerSize(1.4); gd.Draw('p same'); keep.append(gd); lg.AddEntry(gd, 'data ' + nm, 'p')
        LSN = sims[nom][0 if end == 'S' else 1]; gm = ROOT.TGraph(4, np.arange(1, 5, dtype=float) + 0.15, np.array([LSN[L, c_] / LSN[L, 2] for L in range(4)])); gm.SetMarkerStyle(mk + 4); gm.SetMarkerSize(1.4); gm.SetMarkerColor(ROOT.kRed + 1); gm.Draw('p same'); keep.append(gm); lg.AddEntry(gm, 'model ' + nm + ' (y = %+.0f mm)' % nom, 'p')
    lg.Draw(); keep.append(lg)
c.SaveAs('cell_equalisation%s.png' % TAG); c.SaveAs('cell_equalisation%s.pdf' % TAG); say('wrote cell_equalisation%s.png' % TAG)
