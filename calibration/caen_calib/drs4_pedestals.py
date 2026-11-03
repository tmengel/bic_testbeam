"""Per-channel, per-DRS4-cell pedestal table for the CAEN DT5742 (run 3062) from the data themselves.
Sample s of an event whose group start cell is c0 was stored in physical cell (c0 + s) mod 1024.  For each of the 16
channels (12 Baby BCAL + 4 paddles) the pedestal of every cell is the clipped mean of that cell's values over the
pre-pulse samples 0-19 of many events (the start cell is uniformly distributed, so every cell is covered).
    python3 drs4_pedestals.py [nevents]  -> drs4_pedestals.npz (ped[16,1024], n[16,1024]) and a residual-noise summary
"""
import uproot, numpy as np, sys
N = int(sys.argv[1]) if len(sys.argv) > 1 else 150000
f = uproot.open('/sphenix/tg/tg01/jets/tmengel/eic/caen/caen_3062_seg1.root')
hd = f['header'].arrays(['g_index_cell'], entry_stop=N, library='np')['g_index_cell']        # [n,2] start cell per group
bb = f['bbcal'].arrays(['adc'], entry_stop=N, library='np')['adc']                             # [n,12,256]
pd = f['paddles'].arrays(['P3_adc', 'P4_adc', 'P5_adc', 'P6_adc'], entry_stop=N, library='np')
W = np.concatenate([bb, np.stack([pd[k] for k in ('P3_adc', 'P4_adc', 'P5_adc', 'P6_adc')], axis=1)], axis=1).astype(np.int32)   # [n,16,256]
grp = np.array([0] * 8 + [1] * 8)                                                               # DT5742: ch 0-7 group 0, 8-15 group 1
n, nch, ns = W.shape
PRE = 20
S = np.zeros((nch, 1024)); S2 = np.zeros((nch, 1024)); C = np.zeros((nch, 1024))
cells = (hd[:, :, None] + np.arange(PRE)[None, None, :]) % 1024                                # [n,2,PRE]
for ch in range(nch):
    cc = cells[:, grp[ch], :].ravel(); v = W[:, ch, :PRE].ravel().astype(float)
    np.add.at(S[ch], cc, v); np.add.at(S2[ch], cc, v * v); np.add.at(C[ch], cc, 1)
ped = S / np.maximum(C, 1); rms = np.sqrt(np.maximum(S2 / np.maximum(C, 1) - ped ** 2, 0))
# second pass with 3-sigma clipping against the first-pass pedestal (removes pulses that start inside the pre-window)
S[:] = 0; S2[:] = 0; C[:] = 0
for ch in range(nch):
    cc = cells[:, grp[ch], :].ravel(); v = W[:, ch, :PRE].ravel().astype(float); keep = np.abs(v - ped[ch][cc]) < 3 * np.maximum(rms[ch][cc], 5)
    np.add.at(S[ch], cc[keep], v[keep]); np.add.at(S2[ch], cc[keep], v[keep] ** 2); np.add.at(C[ch], cc[keep], 1)
ped = S / np.maximum(C, 1); rms = np.sqrt(np.maximum(S2 / np.maximum(C, 1) - ped ** 2, 0))
print('events %d, entries per cell: min %d mean %.0f' % (n, C.min(), C.mean()))
print('per-channel pedestal spread across cells (std of ped over cells) and mean within-cell rms:')
for ch in range(nch): print('  ch%2d  ped mean %7.1f  cell-to-cell std %5.1f  min %6.0f max %6.0f   within-cell rms %5.1f' % (ch, ped[ch].mean(), ped[ch].std(), ped[ch].min(), ped[ch].max(), rms[ch].mean()))
# validation on the tail region (samples 200-255) of the SAME events: raw vs cell-corrected noise, and a sample-index dependence check
tailcells = (hd[:, :, None] + np.arange(200, 256)[None, None, :]) % 1024
raw_rms = []; cor_rms = []; cor_mean = []
for ch in range(nch):
    cc = tailcells[:, grp[ch], :]; v = W[:, ch, 200:256].astype(float); corr = v - ped[ch][cc]
    raw_rms.append((v - np.median(v, axis=1, keepdims=True)).std()); cor_rms.append(corr.std()); cor_mean.append(corr.mean())
print('tail samples 200-255: rms around the per-event median, raw: %s' % np.round(raw_rms, 1)); print('                      rms after cell pedestal subtraction: %s' % np.round(cor_rms, 1)); print('                      mean residual after subtraction:    %s' % np.round(cor_mean, 1))
# sample-index dependence of the residual (full waveform, channel 0 and a paddle), to see whether a per-sample term is also needed
allcells = (hd[:, :, None] + np.arange(256)[None, None, :]) % 1024
for ch in (0, 7, 12):
    cc = allcells[:, grp[ch], :]; corr = W[:, ch, :].astype(float) - ped[ch][cc]; prof = np.median(corr, axis=0)
    print('ch%2d median residual vs sample index (every 16th sample): %s' % (ch, np.round(prof[::16], 1)))
np.savez('drs4_pedestals.npz', ped=ped, rms=rms, n=C, nevents=n)
print('wrote drs4_pedestals.npz')
