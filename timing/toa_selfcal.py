"""Self-calibration of the H2GCROC3 TOA on the KCU1 single-SiPM channels (no ADC reference).

For every channel k a table O_k[coarse, mid] maps the TOA fields to a time offset so that
  t_k = 25 x TOA sample + O_k[coarse, mid] (+ fine / 8 of the local mid step, ignored here)
agrees between all channels of the same waveform (they see the same shower; light-arrival differences
within one SFIL5 layer are small). The tables are fitted by iterated per-bin medians on the even
waveforms and the resolution is measured on the odd waveforms as the MAD of t_a - t_b between channel
pairs, timewalk-corrected with 1/sqrt(A). The ADC interpolation is not used as a reference: on these
fast pulses (rise < 1 sample) it is pinned to the sampling phase.

Input:  output/kcu1time_0000.root (ts_dump kcu1time, KCU1 ch 144-179)
Output: out/toa_selfcal.log
    python3 toa_selfcal.py
"""
import numpy as np, uproot

LOG = open('out/toa_selfcal.log', 'w')
def say(*a):
    s = ' '.join(str(x) for x in a); print(s); LOG.write(s + '\n')
def mad(x):
    x = x[np.isfinite(x)]
    return 1.4826 * np.median(np.abs(x - np.median(x))) if len(x) else np.nan

ch = next(uproot.iterate('../../output/kcu1time_0000.root:kcu1time', ['wf', 'toa'], library='np', step_size=100000))
wf = ch['wf'].astype(np.float32); toa = ch['toa'].astype(np.int32)
ped = wf[:, :, 1:4].mean(2); v = wf - ped[:, :, None]; A = v[:, :, 4:13].max(2)
N = len(A)
GOOD = [156, 157, 160, 162, 163, 166, 167, 168, 171, 172, 173, 174, 175, 176, 177, 178]
K = [c - 144 for c in GOOD]
AMIN, AMAX = 250, 900

# the TOA word of the pulse: first non-zero TOA in samples 5..9 (the pulse peaks at sample 6-7; thresholds are sane here)
ST = np.full((N, len(K)), -1); CO = np.zeros((N, len(K)), int); MI = np.zeros((N, len(K)), int); FI = np.zeros((N, len(K)), int); OK = np.zeros((N, len(K)), bool); AA = np.zeros((N, len(K)))
for j, k in enumerate(K):
    w = toa[:, k, 5:10]; nz = w != 0; have = nz.any(1); d = nz.argmax(1); code = w[np.arange(N), d]
    OK[:, j] = have & (A[:, k] > AMIN) & (A[:, k] < AMAX)
    ST[:, j] = 5 + d; CO[:, j] = code >> 8; MI[:, j] = (code >> 3) & 31; FI[:, j] = code & 7; AA[:, j] = A[:, k]
say(f'{N} waveforms, {len(K)} channels; waveforms with >= 2 usable channels: {(OK.sum(1) >= 2).sum()}')

train = (np.arange(N) % 2 == 0); test = ~train
NMB = 8; MB = MI // (32 // NMB)        # mid bins (4 codes each): ~45 training hits per (channel, coarse, mid bin)
O = np.zeros((len(K), 4, NMB)); CAL = np.zeros((len(K), 4, NMB), bool)   # offset tables and whether a bin got calibrated
def times(O, use_cal):
    T = np.full((N, len(K)), np.nan)
    for j in range(len(K)):
        m = OK[:, j] & (CAL[j, CO[:, j], MB[:, j]] if use_cal else True)
        T[m, j] = 25.0 * ST[m, j] + O[j, CO[m, j], MB[m, j]]
    return T
for it in range(40):
    T = times(O, False)
    cnt = np.isfinite(T).sum(1); S = np.nansum(T, 1)
    moved = 0.0
    for j in range(len(K)):
        m = OK[:, j] & train & (cnt >= 3)
        ref = (S[m] - T[m, j]) / (cnt[m] - 1)        # mean of the other channels
        r = T[m, j] - ref
        for q in range(4):
            for x in range(NMB):
                b = (CO[m, j] == q) & (MB[m, j] == x)
                if b.sum() >= 8:
                    d = np.median(r[b]); O[j, q, x] -= 0.7 * d; moved = max(moved, abs(d)); CAL[j, q, x] = True
    O -= np.nanmean(O[0][CAL[0]])                   # pin the global offset so the tables do not drift
    if it % 10 == 9 or moved < 0.01: say(f'   iteration {it + 1}: largest bin correction {moved:.3f} ns')
    if moved < 0.01: break
say(f'   calibrated bins: {CAL.sum()} of {CAL.size}; hits in uncalibrated bins are excluded from the test')
T = times(O, True)
say('\nfitted tables (one channel), offset [ns] per mid bin (4 codes each) for each coarse value:')
j0 = K.index(177 - 144)
for q in range(4): say(f'   ch177 coarse {q}: ' + ' '.join(f'{O[j0, q, x]:6.1f}' if CAL[j0, q, x] else '   -  ' for x in range(NMB)))
say('\nresolution on the odd (test) waveforms, pairs of channels, A 250-900 on both, timewalk-corrected:')
res = []
for a in range(len(K)):
    for b in range(a + 1, len(K)):
        m = OK[:, a] & OK[:, b] & test & np.isfinite(T[:, a]) & np.isfinite(T[:, b])
        if m.sum() < 150: continue
        dt = T[m, a] - T[m, b]; g = 1 / np.sqrt(AA[m, a]) - 1 / np.sqrt(AA[m, b])
        core = np.abs(dt - np.median(dt)) < 3 * max(mad(dt), 0.1); p = np.polyfit(g[core], dt[core], 1); dtc = dt - p[0] * g
        res.append((GOOD[a], GOOD[b], m.sum(), mad(dt), mad(dtc)))
res.sort(key=lambda r: r[2], reverse=True)
for a, b, n, m1, m2 in res[:10]:
    say(f'   ch{a} - ch{b}: n={n:5d}  MAD {m1:.2f} ns, timewalk-corrected {m2:.2f} ns -> {m2/np.sqrt(2)*1e3:.0f} ps per channel')
allm = np.array([r[4] for r in res]); alln = np.array([r[2] for r in res])
say(f'   {len(res)} pairs, {alln.sum()} test pairs in total; weighted mean corrected MAD {np.average(allm, weights=alln):.2f} ns -> {np.average(allm, weights=alln)/np.sqrt(2)*1e3:.0f} ps per channel')
# controls: (a) coarse bits only (25 x sample + 6.25 x coarse, no mid), (b) the ADC-calibrated model is in toa_timing.py
say('\ncontrol, coarse bits only (25 x sample + table per coarse value, no mid information):')
Oc = np.zeros((len(K), 4, 1)); MBc = np.zeros_like(MB); CALc = np.ones((len(K), 4, 1), bool)
for it in range(30):
    Tc = np.full((N, len(K)), np.nan)
    for j in range(len(K)):
        m = OK[:, j]; Tc[m, j] = 25.0 * ST[m, j] + Oc[j, CO[m, j], 0]
    cnt = np.isfinite(Tc).sum(1); S = np.nansum(Tc, 1)
    for j in range(len(K)):
        m = OK[:, j] & train & (cnt >= 3); ref = (S[m] - Tc[m, j]) / (cnt[m] - 1); r = Tc[m, j] - ref
        for q in range(4):
            b = CO[m, j] == q
            if b.sum() >= 8: Oc[j, q, 0] -= 0.7 * np.median(r[b])
mm = []
for a, b, n, m1, m2 in res[:10]:
    ia, ib = GOOD.index(a), GOOD.index(b)
    m = OK[:, ia] & OK[:, ib] & test & np.isfinite(Tc[:, ia]) & np.isfinite(Tc[:, ib])
    dt = Tc[m, ia] - Tc[m, ib]; g = 1 / np.sqrt(AA[m, ia]) - 1 / np.sqrt(AA[m, ib])
    core = np.abs(dt - np.median(dt)) < 3 * max(mad(dt), 0.1); p = np.polyfit(g[core], dt[core], 1); mm.append(mad(dt - p[0] * g))
say(f'   mean corrected MAD over the 10 largest pairs: coarse only {np.mean(mm):.2f} ns, coarse + mid bins {np.mean([r[4] for r in res[:10]]):.2f} ns')
LOG.close()
