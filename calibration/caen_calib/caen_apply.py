"""Apply the per-cell DRS4 pedestals (drs4_pedestals.npz) to one CAEN segment of run 3062 and write compact per-event
pulse quantities for the 16 channels (12 Baby BCAL South + 4 paddles):
  base[16]  median of the corrected samples 0-19 (residual per-event offset, subtracted from everything below)
  rms0[16]  rms of the corrected samples 0-19 (noise monitor)
  peak[16]  max over samples 20-230 of (base - corrected)      tpk[16] its sample
  q140[16]  sum over samples 20-159 of (base - corrected)      (same window as events_3062 caen_q)
  qn[16]    sum over samples 35-135 (tight window around the Baby BCAL pulses, which rise at sample ~45)
keyed by header (seg, idx, evt, gtt, g_index_cell).
    python3 caen_apply.py SEG                 -> caen_calib_segSEG.npz           (run 3062 segment)
    python3 caen_apply.py FILE TAG            -> caen_calib_TAG.npz              (any CAEN file, e.g. caen_3088.root 3088, caen_3064_part00.root 3064_part00;
                                                                                  the per-cell pedestals of run 3062 are used, the same DRS4 cells)
"""
import uproot, numpy as np, sys, os
OUT = '/sphenix/tg/tg01/jets/tmengel/eic/caen/calib'; os.makedirs(OUT, exist_ok=True)
if len(sys.argv) > 2: FILE, TAG = sys.argv[1], sys.argv[2]; seg = -1
else: seg = int(sys.argv[1]); FILE, TAG = '/sphenix/tg/tg01/jets/tmengel/eic/caen/caen_3062_seg%d.root' % seg, 'seg%d' % seg
P = np.load('/sphenix/user/tmengel/eic/bic_testbeam/calibration/caen_calib/drs4_pedestals.npz'); ped = P['ped']            # [16,1024]
grp = np.array([0] * 8 + [1] * 8)
f = uproot.open(FILE)
hb, bb, pdl = f['header'], f['bbcal'], f['paddles']; n = hb.num_entries; print('%s: %d events' % (FILE, n), flush=True)
keys = ['seg', 'idx', 'evt', 'gtt']; H = hb.arrays(keys + ['g_index_cell'], library='np')
base = np.zeros((n, 16), np.float32); rms0 = np.zeros((n, 16), np.float32); peak = np.zeros((n, 16), np.float32); tpk = np.zeros((n, 16), np.int16); q140 = np.zeros((n, 16), np.float32); qn = np.zeros((n, 16), np.float32)
CH = 20000; samp = np.arange(256)
for i0 in range(0, n, CH):
    i1 = min(n, i0 + CH)
    W = np.concatenate([bb.arrays(['adc'], entry_start=i0, entry_stop=i1, library='np')['adc'],
                        np.stack([pdl.arrays([k], entry_start=i0, entry_stop=i1, library='np')[k] for k in ('P3_adc', 'P4_adc', 'P5_adc', 'P6_adc')], axis=1)], axis=1).astype(np.float32)   # [m,16,256]
    ic = H['g_index_cell'][i0:i1]                                                     # [m,2]
    cells = (ic[:, grp, None] + samp[None, None, :]) % 1024                            # [m,16,256]
    C = W - ped[np.arange(16)[None, :, None], cells]
    b = np.median(C[:, :, :20], axis=2); D = b[:, :, None] - C                           # positive pulses, per-event offset removed
    base[i0:i1] = b; rms0[i0:i1] = C[:, :, :20].std(axis=2)
    peak[i0:i1] = D[:, :, 20:230].max(axis=2); tpk[i0:i1] = D[:, :, 20:230].argmax(axis=2) + 20
    q140[i0:i1] = D[:, :, 20:160].sum(axis=2); qn[i0:i1] = D[:, :, 35:135].sum(axis=2)
    print('  %d / %d' % (i1, n), flush=True)
np.savez(os.path.join(OUT, 'caen_calib_%s.npz' % TAG), seg=H['seg'], idx=H['idx'], evt=H['evt'], gtt=H['gtt'], icell=H['g_index_cell'], base=base, rms0=rms0, peak=peak, tpk=tpk, q140=q140, qn=qn)
print('wrote', os.path.join(OUT, 'caen_calib_%s.npz' % TAG), ' mean rms0 per channel:', rms0.mean(axis=0).round(1))
