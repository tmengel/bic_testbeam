"""Common-trigger table, run 3062: one row per paddle-OR trigger recorded by any of
KCU0 (complete waveforms), KCU1 (complete waveforms) or CAEN, on the KCU0 clock.

The NIM paddle OR feeds the KCU trigger card and both CAENs, so every trigger carries
a paddle pattern (bits P3=1, P4=2, P5=4, P6=8):
  KCU0: amp[32:36] = pedestal - min(samples 8-16) > 65   (a hit saturates the ADC at 0)
  CAEN: neg[12:16] > 200
KCU0 triggers lost to the decoder as fragments are flagged when a fragment sits at
trigger + k*164 ticks (k = 0..19, +-4 ticks).
-> out/triggers_3062.npz"""
import numpy as np, uproot, glob
OFF1 = 1143665946                      # ts(KCU0) = ts(KCU1) + OFF1 (+0 or -4)
TOL = 16                               # ticks (100 ns) to group the three streams
bits = lambda h: (h * (1 << np.arange(4))).sum(1).astype(np.int8)

K = np.load('out/kcu0_amp.npz'); t0 = K['ts'].astype(np.int64); p0 = bits(K["amp"][:, 32:36] > 65)
A1 = [uproot.open(f)['kcu'].arrays(['ts', 'nsamp'], library='np') for f in sorted(glob.glob('out/beam_ROC_1-00003062-*.root'))]
t1 = np.concatenate([a['ts'][a['nsamp'] == 20] for a in A1]).astype(np.int64)
t1 = np.unique(t1[t1 > 1e8]) + OFF1
M = np.load('out/caen_kcu0_match.npz'); ic = np.nonzero(M['inspan'])[0]; tc = np.round(M['tpred'][ic] / 6.25).astype(np.int64)
neg = np.concatenate([uproot.open(f)['caenadc'].arrays(['neg'], library='np')['neg'] for f in sorted(glob.glob('out/adc_beam_CAEN-00003062-*.root'))])
pC = bits(neg[ic, 12:16] > 200)
print(f"recorded triggers: KCU0 complete {len(t0)}, KCU1 complete {len(t1)}, CAEN in KCU0 span {len(ic)}")

# group the three streams into triggers
t = np.concatenate([t0, t1, tc]); src = np.concatenate([np.zeros(len(t0), int), np.ones(len(t1), int), np.full(len(tc), 2)])
row = np.concatenate([np.arange(len(t0)), np.arange(len(t1)), np.arange(len(tc))])
o = np.argsort(t, kind='stable'); t, src, row = t[o], src[o], row[o]
start = np.r_[True, np.diff(t) > TOL]; g = np.cumsum(start) - 1; ng = g[-1] + 1
i0 = np.full(ng, -1); i1 = np.full(ng, -1); iC = np.full(ng, -1); dup = 0
for s, arr in ((0, i0), (1, i1), (2, iC)):
    k = src == s; dup += k.sum() - len(np.unique(g[k])); arr[g[k][::-1]] = row[k][::-1]   # keep the earliest
T = np.where(i0 >= 0, t0[np.maximum(i0, 0)], np.where(i1 >= 0, t1[np.maximum(i1, 0)], tc[np.maximum(iC, 0)]))
assert (np.diff(T) > 0).all(), 'trigger times not increasing'
pad0 = np.where(i0 >= 0, p0[np.maximum(i0, 0)], -1); padC = np.where(iC >= 0, pC[np.maximum(iC, 0)], -1)
print(f"distinct triggers: {ng}  (streams falling twice into one 100 ns group: {dup})")

# KCU0 fragments at T + k*164
A0 = [uproot.open(f)['kcu'].arrays(['ts', 'nsamp'], library='np') for f in sorted(glob.glob('out/beam_ROC_0-00003062-*.root'))]
fr = np.concatenate([a['ts'][a['nsamp'] < 20] for a in A0]).astype(np.int64); fr = np.sort(fr[fr > 1e8])
frag = np.zeros(ng, bool)
for k in range(20):
    a, b = np.searchsorted(fr, T + k * 164 - 4), np.searchsorted(fr, T + k * 164 + 5); frag |= b > a
frag &= i0 < 0

h0, h1, hC = i0 >= 0, i1 >= 0, iC >= 0
print("\nstream combination            triggers   share")
for name, m in (("KCU0 + KCU1 + CAEN", h0 & h1 & hC), ("KCU0 + CAEN", h0 & ~h1 & hC), ("KCU0 + KCU1", h0 & h1 & ~hC),
                ("KCU0 only", h0 & ~h1 & ~hC), ("KCU1 + CAEN", ~h0 & h1 & hC), ("KCU1 only", ~h0 & h1 & ~hC), ("CAEN only", ~h0 & ~h1 & hC)):
    print(f"  {name:24s} {m.sum():10d}   {m.mean():5.3f}   with KCU0 fragment: {(m & frag).sum()}")

# timing residuals inside a trigger
d10 = t0[i0[h0 & h1]] - t1[i1[h0 & h1]]; d0C = (M['tpred'][ic[iC[h0 & hC]]] / 6.25 - t0[i0[h0 & hC]]) * 6.25
d1C = (M['tpred'][ic[iC[h1 & hC]]] / 6.25 - t1[i1[h1 & hC]]) * 6.25
print("\nKCU0 - KCU1 [ticks]: 0: %.3f  -4: %.3f  other: %.4f" % ((d10 == 0).mean(), (d10 == -4).mean(), (~np.isin(d10, [0, -4])).mean()))
print("CAEN - KCU0 [ns] MAD %.1f;  CAEN - KCU1 [ns] median %.1f MAD %.1f  (KCU1 alone vs CAEN, %d triggers without KCU0: %d)"
      % (np.median(np.abs(d0C - np.median(d0C))), np.median(d1C), np.median(np.abs(d1C - np.median(d1C))), (h1 & hC).sum(), (~h0 & h1 & hC).sum()))

# paddle signature
both = h0 & hC
same = pad0[both] == padC[both]
dec = pad0[both] == np.roll(padC[both], 1)
print(f"\npaddle pattern, KCU0 vs CAEN on {both.sum()} shared triggers: identical {same.mean():.4f}   (decoy, neighbour trigger: {dec.mean():.4f})")
for name, p in (("KCU0", pad0[h0]), ("CAEN", padC[hC])):
    n = np.array([bin(int(x)).count('1') for x in range(16)])[p]
    print(f"  {name}: any paddle {np.mean(n > 0):.4f}   multiplicity 0..4 {np.round(np.bincount(n, minlength=5) / len(n), 4)}"
          f"   P3..P6 {np.round([(p >> b & 1).mean() for b in range(4)], 3)}")
pad = np.where(pad0 >= 0, pad0, padC)
print(f"triggers with a paddle pattern: {(pad >= 0).mean():.4f}   (from KCU0 {h0.mean():.3f}, CAEN only {(~h0 & hC).mean():.3f})")
span = (T[-1] - T[0]) * 6.25e-9
print(f"\nrates over {span:.0f} s: union {ng/span:.0f} Hz;  KCU0 complete {h0.sum()/span:.0f}, +fragments {(h0 | frag).sum()/span:.0f};"
      f"  KCU1 {h1.sum()/span:.0f};  CAEN {hC.sum()/span:.0f} Hz")
np.savez('out/triggers_3062.npz', T=T, i_kcu0=i0, i_kcu1=i1, i_caen=np.where(hC, ic[np.maximum(iC, 0)], -1), pad=pad,
         pad_kcu0=pad0, pad_caen=padC, kcu0_frag=frag)
