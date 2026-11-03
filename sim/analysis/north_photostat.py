"""SFIL5 North: the four live single SiPMs (KCU1 ch 54, 60, 65, 66, D-array block of ASIC0) sit behind a light guide, so they share
the mixed light of their cell (pairwise correlation 0.74-0.83 against 0.10 for two bare South SiPMs). Photostatistics between them
gives the KCU1 single-SiPM channel gain in ADC per photoelectron (x excess noise), the KCU1 counterpart of the SFIL2/3 two-end numbers.
Their correlation with the South D-array sum is only 0.2-0.3: the bare South SiPMs cover about a third of the fibre face and their sum
is a noisy estimate of the cell light, so the South sum cannot be calibrated this way.
    python3 north_photostat.py [ntriggers=500000]  -> north_photostat.txt
"""
import sys, uproot, numpy as np
N = int(sys.argv[1]) if len(sys.argv) > 1 else 500000
out = open('north_photostat.txt', 'w')
def say(*a): print(*a); print(*a, file=out)
f = uproot.open('/sphenix/user/tmengel/eic/bic_testbeam/data/kcu1_3062.root'); tn = f['sfil5north']; ts = f['sfil5']
n = tn.arrays(['ch', 'adc'], entry_stop=N, library='np'); s = ts.arrays(['B_adc', 'D_adc'], entry_stop=N, library='np')
nch = list(n['ch'][0]); NCH = [54, 60, 65, 66]; idx = [nch.index(c) for c in NCH]
def peaks(w):
    w = w.astype(float); good = (w >= 0).all(axis=-1); ped = w[..., 1:4].mean(axis=-1); d = w - ped[..., None]; return d[..., 4:13].max(axis=-1), (w >= 1023).any(axis=-1), good
nA, nSat, nGood = peaks(n['adc'][:, idx, :]); bA, bSat, bGood = peaks(s['B_adc']); dA, dSat, dGood = peaks(s['D_adc'])
sumB, sumD = np.clip(bA, 0, None).sum(axis=1), np.clip(dA, 0, None).sum(axis=1)
shower = nGood.all(axis=1) & bGood.all(axis=1) & dGood.all(axis=1) & (sumB + sumD > 1000)
say('triggers %d, showers %d' % (N, shower.sum()))
for i, c in enumerate(NCH):
    m = shower & ~nSat[:, i] & (nA[:, i] > 60) & (nA[:, i] < 900); mD = m & ~dSat.any(axis=1); mB = m & ~bSat.any(axis=1)
    say('N ch %d: corr with South D sum %.2f (n=%d), with B sum %.2f' % (c, np.corrcoef(nA[mD, i], sumD[mD])[0, 1], mD.sum(), np.corrcoef(nA[mB, i], sumB[mB])[0, 1]))
edges = np.array([60, 90, 120, 160, 200, 250, 300, 360, 430, 500, 580, 670, 770])
def photostat(x, y):
    r = np.median(y / x); yy = y / r; sm = 0.5 * (x + yy); res = {}
    for end, (u, v, sc) in {'x': (x, yy, 1.), 'y': (yy, x, r)}.items():
        A, V, W = [], [], []
        for lo, hi in zip(edges[:-1], edges[1:]):
            m = (sm >= lo) & (sm < hi)
            if m.sum() < 200: continue
            d = u[m] - v[m]; med = np.median(d); mad = 1.4826 * np.median(np.abs(d - med)); k = np.abs(d - med) < 5 * max(mad, 1.)
            uu, vv = u[m][k], v[m][k]; var = uu.var(ddof=1) - np.cov(uu, vv)[0, 1]; A.append(uu.mean() * sc); V.append(var * sc * sc); W.append((len(uu) - 1) / (2 * max(var * sc * sc, 1.) ** 2))
        A, V, W = map(np.array, (A, V, W)); M = np.stack([np.ones_like(A), A], axis=1) * np.sqrt(W)[:, None]; coef, *_ = np.linalg.lstsq(M, V * np.sqrt(W), rcond=None); res[end] = coef[1]
    return r, res
a_est = {c: [] for c in NCH}
for i in range(4):
    for j in range(i + 1, 4):
        m = shower & ~nSat[:, i] & ~nSat[:, j] & (nA[:, i] > 60) & (nA[:, j] > 60) & (nA[:, i] < 900) & (nA[:, j] < 900)
        r, res = photostat(nA[m, i], nA[m, j]); a_est[NCH[i]].append(res['x']); a_est[NCH[j]].append(res['y'])
        say('N%d / N%d: n=%d, corr %.2f, amplitude ratio %.2f, a_%d = %.1f, a_%d = %.1f ADC/pe (x F)' % (NCH[i], NCH[j], m.sum(), np.corrcoef(nA[m, i], nA[m, j])[0, 1], r, NCH[i], res['x'], NCH[j], res['y']))
for c in NCH: say('N ch %d: a = %.1f +- %.1f ADC/pe (x F), mean and spread over its three partners' % (c, np.mean(a_est[c]), np.std(a_est[c])))
say('KCU0 summed channels (two-end photostatistics): SFIL2 S c3 14.1, SFIL3 S c3 9.5 ADC/pe (x F) -> KCU1 single-SiPM chain gain 2.5-3.5 x the KCU0 summed-channel chain')
out.close()
