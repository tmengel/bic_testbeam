"""Interim per-channel Baby BCAL South gains anchored to the v3 model.
Model: bb_ls[segment 0-9][cell 0-3] (South light, MeV) regrouped into the FTBF readout rows (segment 0 | 1-2 | 3-5 | 6-9 =
1 / 2 / 3 / 4 arrays deep) and columns (cells, two orientations tried).  Data: DRS4-cell-corrected CAEN charge qn per channel
(/sphenix/tg/.../caen/calib) for showers, by paddle.  The vertical placement of the Baby BCAL relative to the beam is unknown,
so files with a '_bby<mm>' tag (scan of the block offset) are compared and the offset/orientation that reproduces the
within-row patterns of rows 0 and 1 with equal gain per array is chosen; the gains g_k = <q_k> / <L_k> follow.
    python3 bbcal_gains.py [simdir] [tag]  -> bbcal_gains<tag>.txt/.png
"""
import numpy as np, uproot, glob, re, sys, os, ROOT
sys.path.insert(0, '/sphenix/user/tmengel/eic/bic_testbeam/common'); from plotstyle import setup, legend, label
setup()
SIMDIR = sys.argv[1] if len(sys.argv) > 1 else '/sphenix/tg/tg01/jets/tmengel/eic/sim'; TAG = sys.argv[2] if len(sys.argv) > 2 else ''
out = open('bbcal_gains%s.txt' % TAG, 'w')
def say(*a): print(*a, flush=True); print(*a, file=out)
labels = ['A0', 'A1', 'A3', 'A4', 'A2', 'A5', 'A6', 'A7', 'A8', 'A9', 'A14', 'A11']   # run-log sheet 'Changed on Aug 10': S A0 0, S A1 1, S A3 2, S A4 3, (ch 4 unlabeled: A2 by its pattern), S A5 5, 6-9 -> A6-A9, 14 -> 10, 11 -> 11
chan = {l: int(l[1:]) for l in labels}; row = {l: chan[l] // 4 for l in labels}; col = {l: chan[l] % 4 for l in labels}; narr = {l: row[l] + 1 for l in labels}
dead = {'A3', 'A4'}; live = [l for l in labels if l not in dead]
ROWSEG = {0: [0], 1: [1, 2], 2: [3, 4, 5], 3: [6, 7, 8, 9]}
# ---- data
NSEG = [687701, 694330, 125979]; OFF = np.cumsum([0] + NSEG); NC = OFF[-1]
qn = np.full((NC, 16), np.nan, np.float32)
for s in (0, 1, 2):
    d = np.load('/sphenix/tg/tg01/jets/tmengel/eic/caen/calib/caen_calib_seg%d.npz' % s); qn[OFF[s]:OFF[s + 1]] = d['qn']
TT = np.load('/sphenix/user/tmengel/eic/bic_testbeam/sync/out/triggers_3062.npz'); ic = TT['i_caen']
ev = uproot.open('/sphenix/user/tmengel/eic/bic_testbeam/data/events_3062.root')['events'].arrays(['has_kcu0', 'kcu0_frag', 'pad', 'sfil_s'], library='np')
hc = ic >= 0; S = np.clip(ev['sfil_s'], 0, None).astype(float).sum(axis=(1, 2))
base = (ev['has_kcu0'] == 1) & hc & (ev['kcu0_frag'] == 0) & np.isin(ev['pad'], [1, 2, 4, 8]) & (S > 600)
Q = np.zeros((len(ev['pad']), 12), np.float32); Q[hc] = qn[ic[hc]][:, :12]; Q = np.clip(Q, 0, None)
PB = {4: 2, 3: 1, 6: 8, 5: 4}
dmean = {p: Q[base & (ev['pad'] == PB[p])].mean(axis=0) for p in (4, 3, 6, 5)}; dmean['all'] = Q[base].mean(axis=0)
say('data: %d showers with CAEN; mean corrected charge per channel (all paddles): %s' % (base.sum(), ' '.join('%s %.0f' % (l, v) for l, v in zip(labels, dmean['all']))))
# ---- model files by energy and offset
files = {}
for fn in sorted(glob.glob(SIMDIR + '/bic_sim_e*_s*.root')):
    if '_x80' in fn or '_yb' in fn: continue   # x80: beam under P5; _yb: SFIL beam-height scan, not the Baby BCAL offset
    E = int(re.search(r'_e(\d+)', fn).group(1)); m = re.search(r'_bby(-?[\d.]+)', fn); y = float(m.group(1)) if m else 0.
    files.setdefault((E, y), []).append(fn)
def model_rows(fns):
    L = np.concatenate([uproot.open(f)['sim'].arrays(['bb_ls'], library='np')['bb_ls'][:, :10, :4] for f in fns])    # [n, seg, cell]
    R = np.stack([L[:, segs, :].sum(axis=1) for r, segs in sorted(ROWSEG.items())], axis=1)                            # [n, row, cell]
    return R
say('\nmodel files: ' + ', '.join('%d GeV y=%+.1f (%d)' % (E // 1000, y, len(v)) for (E, y), v in sorted(files.items())))
# ---- orientation / offset choice on rows 0 and 1 (equal gain per array within a row); data pattern from all paddles
def pattern_chi2(Rm, orient):
    chi = 0.; terms = []
    for r in (0, 1):
        ls = [l for l in live if row[l] == r]; dv = np.array([dmean['all'][labels.index(l)] for l in ls]); mv = np.array([Rm[r, col[l] if orient == 0 else 3 - col[l]] for l in ls])
        dv = dv / dv.sum(); mv = mv / mv.sum(); chi += ((dv - mv) ** 2 / np.maximum(mv, 0.02) ** 2).sum() * len(ls); terms.append((r, ls, dv.round(3), mv.round(3)))
    return chi, terms
best = None; results = {}
for (E, y), fns in sorted(files.items()):
    if E != 4000: continue
    R = model_rows(fns); Rm = R.mean(axis=0)
    for orient in (0, 1):
        chi, terms = pattern_chi2(Rm, orient); results[(y, orient)] = (chi, terms, Rm)
        say('  4 GeV, offset %+6.1f mm, orientation %s: pattern chi2 %.2f; ' % (y, 'col = cell' if orient == 0 else 'col = 3 - cell', chi) + ' | '.join('row %d %s data %s model %s' % (r, ls, dv, mv) for r, ls, dv, mv in terms))
        if best is None or chi < best[0]: best = (chi, y, orient)
chi, ybest, obest = best; say('\nchosen: offset %+.1f mm, orientation %s (chi2 %.2f)' % (ybest, 'col = cell' if obest == 0 else 'col = 3 - cell', chi))
# ---- gains at the chosen geometry, per energy available at that offset (else nominal-offset energy dependence as a relative correction)
def Lk(Rm, orient): return {l: Rm[row[l], col[l] if orient == 0 else 3 - col[l]] for l in labels}
Es = sorted({E for (E, y) in files if y == ybest}); Es0 = sorted({E for (E, y) in files if y == 0.})
say('\nper-channel gain g = <q> / <L> [CAEN charge units per MeV of South light], showers of all paddles; per array = g / arrays summed')
hdr = '%-5s %4s %6s ' % ('ch', 'row', 'arrays') + ' '.join('%10s' % ('%d GeV' % (E // 1000)) for E in (Es if len(Es) > 1 else Es0)) + '   per array @4 GeV'
say(hdr)
G = {}
for l in live:
    vals = []
    for E in (Es if len(Es) > 1 else Es0):
        if len(Es) > 1: Rm = model_rows(files[(E, ybest)]).mean(axis=0)
        else:   # scale the chosen-offset 4 GeV light by the nominal-offset energy dependence of that (row, col)
            Rm = results[(ybest, obest)][2] * model_rows(files[(E, 0.)]).mean(axis=0) / model_rows(files[(4000, 0.)]).mean(axis=0)
        vals.append(dmean['all'][labels.index(l)] / Lk(Rm, obest)[l])
    G[l] = vals[(Es if len(Es) > 1 else Es0).index(4000)]
    say('%-5s %4d %6d ' % (l, row[l], narr[l]) + ' '.join('%10.0f' % v for v in vals) + '   %8.0f' % (G[l] / narr[l]))
pa = np.array([G[l] / narr[l] for l in live]); say('per-array gain: mean %.0f, rms/mean %.2f (FTBF MIP peaks per array varied by about 15 %%)' % (pa.mean(), pa.std() / pa.mean()))
# per-paddle check: with the all-paddle gains, how do the per-paddle data patterns compare with the model (x = 30 only)
say('\nper-paddle ratio of the gain-corrected data to the model per channel (4 GeV, chosen offset; 1 = perfect):')
Rm = results[(ybest, obest)][2]; L4 = Lk(Rm, obest)
for p in (4, 3, 6, 5): say('  P%d: ' % p + ' '.join('%s %.2f' % (l, dmean[p][labels.index(l)] / G[l] / L4[l]) for l in live))
np.savez('bbcal_gains%s.npz' % TAG, labels=live, gain=np.array([G[l] for l in live]), yoffset=ybest, orient=obest)
# ---- figure: data vs model pattern per row at the chosen geometry; gains per channel
c = ROOT.TCanvas('c', '', 2000, 800); c.Divide(2, 1); keep = []
c.cd(1); fr = ROOT.gPad.DrawFrame(-0.5, 0, 11.5, 0.45, ';Baby BCAL South channel;share of the live-channel sum'); keep.append(fr)
ax = fr.GetXaxis(); [ax.SetBinLabel(ax.FindBin(i), labels[i]) for i in range(12)]; ax.SetLabelSize(0.05)
dv = dmean['all'].copy(); dv[[labels.index(l) for l in dead]] = 0; dv /= dv.sum(); mv = np.array([L4[l] * (0 if l in dead else 1) for l in labels]); mv /= mv.sum()
hd = ROOT.TH1F('hd', '', 12, -0.5, 11.5); hm = ROOT.TH1F('hm', '', 12, -0.5, 11.5)
for i in range(12): hd.SetBinContent(i + 1, dv[i]); hm.SetBinContent(i + 1, mv[i])
hd.SetMarkerStyle(20); hd.SetMarkerSize(1.4); hd.Draw('p same'); hm.SetLineColor(ROOT.kRed + 1); hm.SetLineWidth(2); hm.Draw('hist same'); keep += [hd, hm]
lg = legend(0.55, 0.7, 0.88, 0.88); lg.AddEntry(hd, 'data, all paddles', 'p'); lg.AddEntry(hm, 'model 4 GeV, equal gain per array, offset %+.0f mm' % ybest, 'l'); lg.Draw(); keep.append(lg)
c.cd(2); fr2 = ROOT.gPad.DrawFrame(-0.5, 0, 11.5, max(pa) * 1.6, ';Baby BCAL South channel;gain per array [CAEN charge units / MeV South light]'); keep.append(fr2)
ax2 = fr2.GetXaxis(); [ax2.SetBinLabel(ax2.FindBin(i), labels[i]) for i in range(12)]; ax2.SetLabelSize(0.05)
hg = ROOT.TH1F('hg', '', 12, -0.5, 11.5)
for l in live: hg.SetBinContent(labels.index(l) + 1, G[l] / narr[l])
hg.SetMarkerStyle(21); hg.SetMarkerSize(1.4); hg.Draw('p same'); keep.append(hg); keep.append(label('model-anchored interim gains; dead channels A3, A4 omitted', x=0.15, y=0.93, size=0.032))
c.SaveAs('bbcal_gains%s.png' % TAG); c.SaveAs('bbcal_gains%s.pdf' % TAG); say('wrote bbcal_gains%s.png' % TAG)
