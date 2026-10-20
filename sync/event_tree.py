"""Per-trigger event tree for run 3062 -> out/events_3062.root, tree "events".

One entry per row of the common-trigger table (trigger_table.py), in time order.
A stream that missed the trigger has has_<stream> = 0 and zeros in its branches.

  trig_ts            trigger time, KCU0 FPGA ticks (6.25 ns); t_s = seconds since the first trigger
  has_kcu0/kcu1/caen stream recorded the trigger (complete waveform for the KCUs)
  kcu0_frag          KCU0 saw the trigger only as decoder fragments (no amplitudes here)
  pad                paddle pattern, bits P3=1 P4=2 P5=4 P6=8 (KCU0, else CAEN; -1 none)
  pad_kcu0, pad_caen the pattern from each stream (-1 if that stream missed the trigger)
  sfil_s[4][4]       KCU0 ASIC0, SFIL1-4 South, [cell][layer], max - pedestal (ADC)
  sfil_n[4][4]       KCU0 ASIC2, SFIL1-4 North, same layout
                     cell 3 = SiPM3 (beam row), cell 4 = SiPM4; cells 1 and 2 both read SiPM1+2 on
                     v2 boards (all South, North SFIL2-4); North SFIL1 is v1 (one SiPM per cell)
  pad_dip[4]         KCU0 paddle P3-P6, pedestal - min(samples 8-16); a hit saturates at ~80
  sfil_s_w, sfil_n_w same layout, pedestal-subtracted window amplitude: mean(samples 6-9) - mean(samples 1-3)
                     - the ASIC common mode (median unconnected channel). Noise averages to 0 here,
                     while max - pedestal (sfil_s, sfil_n) adds +10-20 ADC per channel with no signal
  sfil5_s[72]        KCU1 ASIC2 ch 144-215, max - pedestal; arrays A..D = ch 144-160, 162-178, 180-196,
                     198-214 minus every 9th (152, 161, ...: not connected)
  sfil5_s_w[72]      same channels, mean(samples 5-8) - mean(samples 1-3) - common mode
  caen[16]           CAEN pulse height (pedestal - min, ADC): ch 0-11 Baby BCAL South, 12-15 = P3-P6
  caen_q[16]         CAEN pulse integral, sum over samples 20-159 of (pedestal - ADC): zero mean on noise
  n_astro            AstroPix clusters matched to this trigger (log10LR > 1)
  astro_*            best of them: x, y [mm, nominal pitch; NaN if no col/row hit], chip, n_hits,
                     log10LR, purity (per-match, from the control bands), latency_ms
"""
import numpy as np, uproot, glob, sys
sys.path.insert(0, '.')
from astro_geo_common import cluster_xy

TT = np.load('out/triggers_3062.npz'); N = len(TT['T'])
i0, i1, ic = TT['i_kcu0'], TT['i_kcu1'], TT['i_caen']
h0, h1, hc = i0 >= 0, i1 >= 0, ic >= 0
print("triggers:", N)

Kd = np.load('out/kcu0_amp.npz'); amp = Kd['amp']
a0 = np.zeros((N, 36), np.float32); a0[h0] = amp[i0[h0]]; del amp
# pedestal-subtracted window amplitudes, same row order as kcu0_amp.npz (sorted unique timestamps)
W = [uproot.open(f)['kcuadc'].arrays(['ts', 'win', 'cm'], library='np') for f in sorted(glob.glob('out/adc_beam_ROC_0-00003062-*.root'))]
wt = np.concatenate([w['ts'] for w in W]).astype(np.int64); o = np.argsort(wt, kind='stable'); u = np.r_[True, np.diff(wt[o]) > 0]
assert np.array_equal(wt[o][u], Kd['ts']), "kcuadc files differ from kcu0_amp.npz: rerun astro_physics.py"
ww = np.concatenate([w['win'] for w in W])[o][u]; wc = np.concatenate([w['cm'] for w in W])[o][u]; del W
side = np.r_[np.zeros(16, int), np.full(16, 2)]
w0 = np.zeros((N, 32), np.float32); w0[h0] = (ww[:, :32] - wc[:, side])[i0[h0]]; del ww, wc

# KCU1: rows of i_kcu1 are the unique complete-waveform timestamps (as trigger_table.py)
A1 = [uproot.open(f)['kcu'].arrays(['ts', 'nsamp'], library='np') for f in sorted(glob.glob('out/beam_ROC_1-00003062-*.root'))]
tt1 = np.concatenate([a['ts'][a['nsamp'] == 20] for a in A1]).astype(np.int64); tt1 = np.unique(tt1[tt1 > 1e8])
K1 = [uproot.open(f)['kcu1adc'].arrays(['ts', 'a1', 'w1', 'cm1'], library='np') for f in sorted(glob.glob('out/adc_beam_ROC_1-00003062-*.root'))]
t1 = np.concatenate([k['ts'] for k in K1]).astype(np.int64); a1 = np.concatenate([k['a1'] for k in K1])
b1 = np.concatenate([k['w1'] - k['cm1'][:, None] for k in K1]).astype(np.float32)
o = np.argsort(t1, kind='stable'); t1, a1, b1 = t1[o], a1[o], b1[o]; u = np.r_[True, np.diff(t1) > 0]; t1, a1, b1 = t1[u], a1[u], b1[u]
x = tt1[i1[h1]]; j = np.clip(np.searchsorted(t1, x), 0, len(t1) - 1); ok = t1[j] == x
s5 = np.zeros((N, 72), np.int16); r1 = np.nonzero(h1)[0]; s5[r1[ok]] = a1[j[ok]]
s5w = np.zeros((N, 72), np.float32); s5w[r1[ok]] = b1[j[ok]]
print("KCU1 rows with amplitudes: %d of %d" % (ok.sum(), h1.sum()))

C = [uproot.open(f)['caenadc'].arrays(['neg', 'qint'], library='np') for f in sorted(glob.glob('out/adc_beam_CAEN-00003062-*.root'))]
neg = np.concatenate([c['neg'] for c in C]); qint = np.concatenate([c['qint'] for c in C]); del C
cn = np.zeros((N, 16), np.int16); cn[hc] = neg[ic[hc]]; cq = np.zeros((N, 16), np.float32); cq[hc] = qint[ic[hc]]

# AstroPix: best match per trigger
M = np.load('out/astro_kcu0_matches.npz'); chip, row, col, X, Y = cluster_xy()
o = np.lexsort((-M['log10LR'], M['trig_row'])); tr = M['trig_row'][o]; first = np.r_[True, np.diff(tr) > 0]
best = o[first]; nA = np.bincount(M['trig_row'], minlength=N).astype(np.int16)
ax = np.full(N, np.nan, np.float32); ay = ax.copy(); alr = ax.copy(); apu = ax.copy(); alat = ax.copy()
achip = np.full(N, -1, np.int16); anh = np.zeros(N, np.int16)
cb = M['cluster'][best]; rb = M['trig_row'][best]
ax[rb], ay[rb], achip[rb], anh[rb] = X[cb], Y[cb], chip[cb], M['n_hits'][best]
alr[rb], apu[rb], alat[rb] = M['log10LR'][best], M['purity'][best], M['latency_ms'][best]
print("triggers with an AstroPix match: %d (expected true %.0f)" % (len(rb), apu[rb].sum()))

T = TT['T'].astype(np.int64)
cols = {'trig_ts': T, 't_s': ((T - T[0]) * 6.25e-9).astype(np.float64),
        'has_kcu0': h0.astype(np.int8), 'has_kcu1': h1.astype(np.int8), 'has_caen': hc.astype(np.int8),
        'kcu0_frag': TT['kcu0_frag'].astype(np.int8), 'pad': TT['pad'].astype(np.int8),
        'pad_kcu0': TT['pad_kcu0'].astype(np.int8), 'pad_caen': TT['pad_caen'].astype(np.int8),
        'sfil_s': a0[:, :16].reshape(N, 4, 4), 'sfil_n': a0[:, 16:32].reshape(N, 4, 4), 'pad_dip': a0[:, 32:36],
        'sfil_s_w': w0[:, :16].reshape(N, 4, 4), 'sfil_n_w': w0[:, 16:32].reshape(N, 4, 4),
        'sfil5_s': s5, 'sfil5_s_w': s5w, 'caen': cn, 'caen_q': cq, 'n_astro': nA, 'astro_x': ax, 'astro_y': ay, 'astro_chip': achip,
        'astro_nhits': anh, 'astro_log10LR': alr, 'astro_purity': apu, 'astro_latency_ms': alat}
with uproot.recreate('out/events_3062.root', compression=uproot.ZLIB(4)) as f:
    t = f.mktree('events', {k: (v.dtype if v.ndim == 1 else np.dtype((v.dtype, v.shape[1:]))) for k, v in cols.items()},
                 title='BIC run 3062 common-trigger events')
    for s in range(0, N, 500000):
        t.extend({k: v[s:s + 500000] for k, v in cols.items()})
print("wrote out/events_3062.root:", N, "entries")
