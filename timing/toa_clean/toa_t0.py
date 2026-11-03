"""Independent event time for the SFIL2 cell 3 TOA North-South study.
The KCU latches the trigger on the 25 ns sample clock (ts0 is always 3 mod 4), so where the particle falls inside the sample
is not recorded by the KCU, and t_N and t_S share that phase: the N-vs-S scatter lies on the +45 deg diagonal. The CAEN
group trigger time (8.533 ns ticks), mapped onto the KCU0 clock by the run-3062 sync (caen_kcu0.py), does see it:
    res = 6.25 ts0 - t_CAEN   (uniform over ~25 ns, MAD 6.4 ns)
so -res is the trigger time relative to the KCU window, with ~2.5 ns rms from the CAEN tick. The particle is not at a fixed
time from the trigger: the CAEN paddle pulse sits at a spread of +-6.6 ns in the DRS4 window, which the trigger stops
(caen_padtime.py). The event time is therefore
    t0 = -res + t_pad(CAEN)     (CAEN paddle 50 % crossing, ns in the window)
The KCU0 paddles do not help: a fired paddle drives the input to 0 ADC within one sample and never writes a TOA word.
    python3 toa_t0.py -> toa_t0.txt, toa_t0.npz (same entries as toa_clean.npz)"""
import numpy as np, uproot
out = open('toa_t0.txt', 'w')
def say(*a): print(*a, flush=True); print(*a, file=out)
z = np.load('toa_clean.npz'); NT = len(z['pad']); even = z['even']; pad = z['pad']
PB = {4: 2, 3: 1, 6: 8, 5: 4}; PX = {4: -80., 3: -30., 6: 30., 5: 80.}
ts0 = uproot.open('/sphenix/user/tmengel/eic/bic_testbeam/data/kcu0_3062.root')['header']['ts0'].array(library='np', entry_stop=NT).astype(np.int64)
S = '/sphenix/user/tmengel/eic/bic_testbeam/sync/out/'
m = np.load(S + 'caen_kcu0_match.npz'); kts = np.load(S + 'kcu0_ts.npy')[m['ikcu0']]
o = np.argsort(kts); i = np.clip(np.searchsorted(kts[o], ts0), 0, len(kts) - 1); hit = kts[o][i] == ts0
ph = np.where(hit, -m['res'][o][i], np.nan)                      # trigger phase in the KCU window
cp = np.load('caen_padtime.npz')['tp']; ic = np.where(hit, m['icaen'][o][i], 0)
tpc = np.where(hit, cp[ic], np.nan); tpc -= np.nanmedian(tpc)      # particle time relative to the trigger
t0 = ph + tpc
say('KCU0 entries %d, with a CAEN partner within 1 us: %.3f; ts0 mod 4 values %s' % (NT, hit.mean(), np.unique(ts0 % 4).tolist()))
b = z['both'] & ~even & np.isfinite(t0); tN = z['tN'][b]; tS = z['tS'][b]; tt = t0[b]; pb = pad[b]
for lab, x in (('trigger phase only', ph[b]), ('CAEN paddle in window only', tpc[b]), ('both', tt)):
    say('  %-27s: slope of (t_N + t_S)/2 %.3f, correlation %.3f, spread of (t_N + t_S)/2 - x %.2f ns'
        % (lab, np.polyfit(x, (tN + tS) / 2, 1)[0], np.corrcoef(x, (tN + tS) / 2)[0, 1], (lambda v: (np.percentile(v, 84) - np.percentile(v, 16)) / 2)((tN + tS) / 2 - x)))
r = lambda x: (np.percentile(x, 84) - np.percentile(x, 16)) / 2
s = (tN + tS) / 2; d = (tN - tS) / 2
sl = np.polyfit(tt, s, 1)[0]
say('odd-trigger showers with TOA at both ends and a CAEN partner: %d of %d' % (b.sum(), (z['both'] & ~even).sum()))
say('slope of (t_N + t_S)/2 against t0: %.3f; correlation %.3f; (t_N - t_S)/2 against t0: correlation %.3f' % (sl, np.corrcoef(tt, s)[0, 1], np.corrcoef(tt, d)[0, 1]))
say('spread (half 16-84 %%): (t_N + t_S)/2 %.2f ns -> minus t0 %.2f ns; (t_N - t_S)/2 %.2f ns' % (r(s), r(s - tt), r(d)))
say('correlation t_N vs t_S %.3f -> t_N - t0 vs t_S - t0 %.3f' % (np.corrcoef(tN, tS)[0, 1], np.corrcoef(tN - tt, tS - tt)[0, 1]))
for p in (4, 3, 6, 5):
    mm = pb == PB[p]
    say('P%d x %+4.0f mm: median t_N - t0 %7.2f, t_S - t0 %7.2f, (N + S)/2 - t0 %7.2f, N - S %+.2f ns (%d showers)'
        % (p, PX[p], np.median((tN - tt)[mm]), np.median((tS - tt)[mm]), np.median((s - tt)[mm]), np.median((tN - tS)[mm]), mm.sum()))
np.savez('toa_t0.npz', t0=t0, ph=ph, tpc=tpc)
