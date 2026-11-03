"""CAEN paddle pulse time in the DRS4 window (1 GS/s, window stopped by the trigger): the particle time relative to the
trigger, per CAEN entry, for the largest paddle dip. 50 % crossing of the leading edge, linear between samples.
    python3 caen_padtime.py -> caen_padtime.npz (tp [ns], amp, which paddle 3-6; NaN if no dip > 200 ADC)"""
import numpy as np, uproot
f = uproot.open('/sphenix/user/tmengel/eic/bic_testbeam/data/caen_3062.root')
n = f['paddles'].num_entries; tp = np.full(n, np.nan); amp = np.zeros(n); which = np.zeros(n, np.int8); i0 = 0
for a in f['paddles'].iterate(['P%d_adc' % p for p in (3, 4, 5, 6)], step_size=50000, library='np'):
    W = np.stack([a['P%d_adc' % p].astype(float) for p in (3, 4, 5, 6)], 1)      # (n, 4, samples)
    y = np.median(W[:, :, :20], axis=2)[:, :, None] - W
    A = y[:, :, 20:230].max(axis=2); j = A.argmax(axis=1); m = len(j); r = np.arange(m)
    yy = y[r, j]; AA = A[r, j]; ipk = yy[:, 20:230].argmax(axis=1) + 20
    above = yy >= 0.5 * AA[:, None]; idx = np.arange(yy.shape[1])[None, :]
    k = np.where(above & (idx <= ipk[:, None]) & (idx >= 20), idx, 9999).min(axis=1); k = np.clip(k, 21, yy.shape[1] - 1)
    y0 = yy[r, k - 1]; y1 = yy[r, k]; t = (k - 1) + np.clip((0.5 * AA - y0) / np.maximum(y1 - y0, 1e-3), 0, 1)
    ok = AA > 200; tp[i0:i0 + m] = np.where(ok, t, np.nan); amp[i0:i0 + m] = AA; which[i0:i0 + m] = np.array((3, 4, 5, 6))[j]
    i0 += m
np.savez('/sphenix/user/tmengel/eic/bic_testbeam/timing/toa_clean/caen_padtime.npz', tp=tp, amp=amp, which=which)
print('entries %d, with a paddle dip > 200: %.3f; crossing sample percentiles 1/16/50/84/99: %s' % (n, np.isfinite(tp).mean(), np.round(np.nanpercentile(tp, [1, 16, 50, 84, 99]), 2)))
