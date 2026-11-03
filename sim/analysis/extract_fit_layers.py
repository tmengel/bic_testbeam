"""Per-event SFIL layer sums from the full-run template fits (kcu0_sfil_fit_full.root) -> data_fit_layers.npz.
    python3 extract_fit_layers.py [fit file] [output npz]     (default: run 3062 full fits -> data_fit_layers.npz; e.g. the 3064 parts -> data_fit_layers_3064_partPP.npz)
Per channel the amplitude is the template fit (saturation recovered) when the fit was accepted, else the raw
pedestal-subtracted peak minus that channel's noise-peak bias (the mean raw peak in triggers without any pulse,
~5 ADC: the maximum of nine noise samples). North SFIL1 is masked (0).
Arrays: pad; S[n,4], N[n,4] layer sums over the four cells; Sraw/Nraw raw peaks only (no fit, no bias correction);
Ssat/Nsat the part of S/N from saturated channels (peak sample at the 10-bit ceiling, fitted); Snfit the number of South channels taken from the raw peak;
satS/satN[n,4] number of saturated channels per layer and end; S3/N3 the beam-row (cell 3) amplitude per layer and end, Sraw3 its raw peak;
fullfit (every channel with a raw peak >= 60 ADC has an accepted fit); nsatfail (raw peak >= 1000 without a fit); bias[32].
"""
import uproot, numpy as np, sys
F = sys.argv[1] if len(sys.argv) > 1 else '/sphenix/tg/tg01/jets/tmengel/eic/kcu0/fit_full/kcu0_sfil_fit_full.root'
OUT = sys.argv[2] if len(sys.argv) > 2 else 'data_fit_layers.npz'
t = uproot.open(F)['fits']
# noise-peak bias per channel from the first 400k triggers: triggers where no unmasked channel has a raw peak >= 25 ADC
a = t.arrays(['amp', 'masked'], entry_stop=400000, library='np'); raw = np.nan_to_num(a['amp']); m = a['masked'] == 1; raw[m] = 0
empty = ~np.any(raw >= 25, axis=1); bias = raw[empty].mean(axis=0); bias[m[0]] = 0
print('empty triggers for the bias: %d of 400k; bias per channel (ADC):' % empty.sum()); print(bias.round(1).reshape(4, 2, 4))
out = {k: [] for k in ('trig', 'pad', 'S', 'N', 'Sraw', 'Nraw', 'Ssat', 'Nsat', 'Snfit', 'fullfit', 'nsatfail', 'satS', 'satN', 'S3', 'N3', 'Sraw3')}
for a in t.iterate(['trig', 'pad', 'ok', 'amp', 'afit', 'masked', 'ped'], step_size=200000, library='np'):
    ok = a['ok']; raw = np.nan_to_num(a['amp']); fit = np.nan_to_num(a['afit']); m = a['masked'] == 1; raw[m] = 0
    satur = (raw + np.nan_to_num(a['ped'])) >= 1020.   # the peak sample sits at the 10-bit ceiling
    use = np.where(ok, fit, raw - bias); use[m] = 0
    big = (raw >= 60) & ~m
    out['trig'].append(a['trig']); out['pad'].append(a['pad']); out['fullfit'].append(~np.any(big & ~ok, axis=1)); out['nsatfail'].append(np.sum(satur & ~ok, axis=1))
    k = np.arange(32); L = k // 8; e = (k // 4) % 2
    sat = np.where(ok & satur, fit, 0.)                      # fitted amplitude of saturated channels (for the recovery systematic)
    nfit = np.where(~ok & ~m, 1., 0.)                                # channels taken from the raw peak (for the noise-bias systematic)
    satf = satur.astype(float)                                       # 1 where the channel's peak sample is at the ceiling (fitted or not)
    for name, arr in (('S', use), ('N', use), ('Sraw', raw), ('Nraw', raw), ('Ssat', sat), ('Nsat', sat), ('Snfit', nfit), ('satS', satf), ('satN', satf)):
        sel = e == (0 if name in ('S', 'Sraw', 'Ssat', 'Snfit', 'satS') else 1)   # end by name (a startswith('S') test put satS on the North end once)
        out[name].append(np.stack([arr[:, sel & (L == l)].sum(axis=1) for l in range(4)], axis=1))
    out['S3'].append(np.stack([use[:, 8 * l + 2] for l in range(4)], axis=1)); out['N3'].append(np.stack([use[:, 8 * l + 6] for l in range(4)], axis=1)); out['Sraw3'].append(np.stack([raw[:, 8 * l + 2] for l in range(4)], axis=1))
    print('.', end='', flush=True)
out = {k: np.concatenate(v) for k, v in out.items()}; out['bias'] = bias; np.savez(OUT, **out); print('\nsaved', OUT, len(out['pad']), 'events')
tot = out['S'].sum(axis=1)
for cut in (300, 600, 1200):
    mm = tot > cut; print('cut S total > %4d: %7d events (%.2f), South shares' % (cut, mm.sum(), mm.mean()), (out['S'][mm].mean(axis=0) / out['S'][mm].mean(axis=0).sum()).round(3), 'mean per layer', out['S'][mm].mean(axis=0).round(1))
for p in (3, 4, 5, 6):
    mm = (tot > 600) & (out['pad'] == p); print('P%d: %6d events, South mean per layer' % (p, mm.sum()), out['S'][mm].mean(axis=0).round(1), 'shares', (out['S'][mm].mean(axis=0) / out['S'][mm].mean(axis=0).sum()).round(3), ' sigma/mu of the total %.3f' % (tot[mm].std() / tot[mm].mean()))
