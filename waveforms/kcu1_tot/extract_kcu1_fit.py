"""Chunked pass over the full-run KCU1 SFIL5 fit tree: compact per-trigger arrays for the SFIL5 rebuild and
per-channel (afit, TOT) histograms for the TOT calibration.  Output: kcu1_fit_arrays.npz (on /sphenix/tg) and
kcu1_tot_hists.npz (here)."""
import uproot, numpy as np, sys, os
fn = '/sphenix/tg/tg01/jets/tmengel/eic/kcu1/fit_full/kcu1_sfil5_fit_full.root'
outdir = '/sphenix/tg/tg01/jets/tmengel/eic/kcu1/fit_full'
t = uproot.open(fn)['fits']; n = t.num_entries; print('entries', n, flush=True)
ts0 = np.zeros(n, np.int64); trig = np.zeros(n, np.int64)
afit = np.zeros((n, 64), np.float32); amp = np.zeros((n, 64), np.float32); nsat = np.zeros((n, 64), np.int8); tot = np.zeros((n, 64), np.int16); ok = np.zeros((n, 64), bool)
ae = np.arange(0, 4001, 20.); te = np.arange(-4, 4100, 8.)        # afit bins 20 ADC, TOT bins 8 (one coarse TOT unit)
H0 = np.zeros((64, len(ae) - 1, len(te) - 1), np.int32)            # unsaturated pulses with a TOT word
H1 = np.zeros((4, 64, len(ae) - 1, len(te) - 1), np.int32)         # nsat = 1,2,3,>=4 pulses with a TOT word (afit = template extrapolation)
i0 = 0
for arr in t.iterate(['trig', 'ts0', 'afit', 'amp', 'nsat', 'tot', 'ok'], step_size='200 MB', library='np'):
    m = len(arr['trig']); sl = slice(i0, i0 + m)
    trig[sl] = arr['trig']; ts0[sl] = arr['ts0']; afit[sl] = arr['afit']; amp[sl] = arr['amp']; nsat[sl] = arr['nsat']; tot[sl] = arr['tot']; ok[sl] = arr['ok']
    a, s, w, o = arr['afit'], arr['nsat'], arr['tot'], arr['ok']
    for k in range(64):
        has = o[:, k] & (w[:, k] > 0)
        for cls, sel in ((None, has & (s[:, k] == 0)), (0, has & (s[:, k] == 1)), (1, has & (s[:, k] == 2)), (2, has & (s[:, k] == 3)), (3, has & (s[:, k] >= 4))):
            if not sel.any(): continue
            h, _, _ = np.histogram2d(a[sel, k], w[sel, k], bins=(ae, te))
            if cls is None: H0[k] += h.astype(np.int32)
            else: H1[cls, k] += h.astype(np.int32)
    i0 += m; print('  %d / %d' % (i0, n), flush=True)
np.savez(os.path.join(outdir, 'kcu1_fit_arrays.npz'), trig=trig, ts0=ts0, afit=afit, amp=amp, nsat=nsat, tot=tot, ok=ok)
np.savez('kcu1_tot_hists.npz', ae=ae, te=te, H0=H0, H1=H1)
# quick per-channel summary
print('%4s %7s %7s %7s %7s %7s' % ('k', 'ok', 'tot&ok', 'sat', 'sat&tot', 'frac'))
for k in range(64):
    o = ok[:, k]; s = o & (nsat[:, k] > 0); w = o & (tot[:, k] > 0)
    if o.sum() > 1000: print('%4d %7d %7d %7d %7d %7.3f' % (k, o.sum(), (w & (nsat[:, k] == 0)).sum(), s.sum(), (s & w).sum(), (s & w).sum() / max(1, s.sum())))
print('done')
