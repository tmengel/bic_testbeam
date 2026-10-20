"""Trigger lists for wf_acc (waveform sums), one file per stream: out/wfkeys_{kcu0,kcu1,caen}.txt

Classes 0-7, background subtraction at waveform level: every AstroPix candidate trigger used by
the geometry plots (astro_geo_common: right latency window), split by the cluster's chip region
  group 0 = chip 0 (P4), 1 = chips 2-3 (P3), 2 = chips 5-6 (P6), 3 = chip 8 (P5)
and by band: class 2*group = phase-matched (signal) band, 2*group+1 = control bands. Weight =
number of candidate rows. True-partner sum = s(signal) - 0.5 s(control).
Classes 100+, event displays: single triggers (weight 1), listed in out/evdisp_events.npz: two per
paddle region: the two most confident single AstroPix matches with all three streams (purity
is low on chips 0 and 8, see astro_final.py; the display prints it)."""
import numpy as np, uproot, glob, sys
sys.path.insert(0, '.')
from astro_geo_common import load

G = load(); TT = np.load('out/triggers_3062.npz')
grp = np.full(9, -1); grp[0] = 0; grp[[2, 3]] = 1; grp[[5, 6]] = 2; grp[8] = 3
g = grp[G['chip'][G['cl']]]; ok = g >= 0
cls = 2 * g[ok] + (G['w'][ok] < 0)
key = G['ki'][ok].astype(np.int64) * 16 + cls
u, cnt = np.unique(key, return_counts=True); rows, cls = u // 16, u % 16

# event displays: clean three-stream triggers with one confident AstroPix match, two per paddle region
ev = uproot.open('out/events_3062.root')['events'].arrays(['has_kcu0', 'has_kcu1', 'has_caen', 'n_astro', 'astro_purity',
                                                         'astro_chip', 'astro_x', 'astro_y', 'pad', 'sfil_s', 'sfil_n'], library='np')
good = (ev['has_kcu0'] == 1) & (ev['has_kcu1'] == 1) & (ev['has_caen'] == 1) & (ev['n_astro'] == 1) & np.isfinite(ev['astro_x'])
nsat = ((ev['sfil_s'] > 900) | (ev['sfil_n'] > 900)).sum((1, 2))
pick = []
for p, chips in ((2, [0]), (1, [2, 3]), (8, [5, 6]), (4, [8])):
    c = np.nonzero(good & (ev['pad'] == p) & np.isin(ev['astro_chip'], chips) & (nsat <= 2))[0]
    c = c[np.argsort(-ev["astro_purity"][c], kind="stable")][:2]                  # the two most confident matches
    pick += list(np.sort(c))
pick = np.array(pick); print("event-display triggers:", pick, " chips", ev['astro_chip'][pick], " pads", ev['pad'][pick])
np.savez('out/evdisp_events.npz', rows=pick)

A1 = [uproot.open(f)['kcu'].arrays(['ts', 'nsamp'], library='np') for f in sorted(glob.glob('out/beam_ROC_1-00003062-*.root'))]
tt1 = np.concatenate([a['ts'][a['nsamp'] == 20] for a in A1]).astype(np.int64); tt1 = np.unique(tt1[tt1 > 1e8])
allrows = np.r_[rows, pick]; allcls = np.r_[cls, 100 + np.arange(len(pick))]; allw = np.r_[cnt, np.ones(len(pick), int)]
for name, idx, conv in (("kcu0", TT['i_kcu0'], lambda r: TT['T'][r]), ("kcu1", TT['i_kcu1'], lambda r: tt1[TT['i_kcu1'][r]]),
                        ("caen", TT['i_caen'], lambda r: TT['i_caen'][r])):
    m = idx[allrows] >= 0
    np.savetxt(f'out/wfkeys_{name}.txt', np.c_[conv(allrows[m]), allcls[m], allw[m]], fmt='%d')
    print(f"{name}: {m.sum()} (trigger, class) lines, {len(np.unique(allrows[m]))} triggers")
