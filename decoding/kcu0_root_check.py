"""Checks on the merged KCU0 file: entry counts per tree, trigger order, duplicates, completeness, segment coverage.
    python3 kcu0_root_check.py output/kcu0_root/kcu0_3062.root
"""
import sys, numpy as np, uproot
f = uproot.open(sys.argv[1])
trees = ['header', 'sfil1', 'sfil2', 'sfil3', 'sfil4', 'paddles', 'bbcal', 'unused', 'packets', 'map']
print('entries:', {t: f[t].num_entries for t in trees if t in f})
h = f['header']
a = h.arrays(['ts0', 'seg', 'nsamp', 'nlines', 'complete', 'backward', 'evt_first', 'evt_last', 'seg_last', 'udp_first'], library='np')
ts0 = a['ts0']; n = len(ts0)
if n == 0: sys.exit('no triggers in the header tree')
d = np.diff(ts0)
print('triggers', n, '; complete', int(a['complete'].sum()), '(%.1f%%)' % (100 * a['complete'].mean()), '; nsamp==20:', int((a['nsamp'] == 20).sum()), '; single-sample fragments', int((a['nsamp'] == 1).sum()))
print('ts0 order: backward steps', int((d < 0).sum()), '; duplicate ts0', int((d == 0).sum()), '; flagged backward', int(a['backward'].sum()))
print('ts0 range %.3f .. %.3f s' % (ts0[0] * 6.25e-9, ts0[-1] * 6.25e-9))
for s in np.unique(a['seg']):
    m = a['seg'] == s
    if m.sum(): print('  segment %d: %7d triggers, ts0 %.3f..%.3f s, complete %.1f%%, triggers spanning two RCDAQ events %d, spanning two segments %d' % (s, m.sum(), ts0[m][0] * 6.25e-9, ts0[m][-1] * 6.25e-9, 100 * a['complete'][m].mean(), int((a['evt_first'][m] != a['evt_last'][m]).sum()), int((a['seg_last'][m] != s).sum())))
# boundaries: the trigger before and after each segment change
chg = np.flatnonzero(np.diff(a['seg']) != 0)
for i in chg:
    print('  boundary %d->%d: entries %d/%d ts0 %d -> %d (dt %.1f us), nsamp %d/%d' % (a['seg'][i], a['seg'][i + 1], i, i + 1, ts0[i], ts0[i + 1], (ts0[i + 1] - ts0[i]) * 6.25e-3, a['nsamp'][i], a['nsamp'][i + 1]))
# sample timestamps inside a trigger must step 164 ticks
ts = h['ts'].array(library='np', entry_stop=200000); ns = a['nsamp'][:200000]
bad = 0
for i in range(len(ts)):
    k = min(ns[i], 20)
    if k > 1 and np.any(np.diff(ts[i][:k]) != 164): bad += 1
print('first 200k triggers with a sample step != 164 ticks:', bad)
p = f['packets'].arrays(['seg', 'evt', 'evt_type', 'has_kcu', 'udp_jumps', 'udp_lost', 'nomarker', 'crc_bad', 'ham_flagged', 'dup_lines', 'badlink', 'nwf_emitted'], library='np')
print('packets:', len(p['evt']), '; data events with the KCU packet', int(p['has_kcu'].sum()), '; event types', dict(zip(*np.unique(p['evt_type'], return_counts=True))))
ev = p['evt'][p['has_kcu']]; print('RCDAQ event sequence: min', ev.min(), 'max', ev.max(), 'gaps', int((np.diff(ev) != 1).sum()), 'repeats', int((np.diff(ev) <= 0).sum()))
print('UDP: sequence jumps', int(p['udp_jumps'].sum()), 'datagrams lost', int(p['udp_lost'].sum()), '; marker-less lines', int(p['nomarker'].sum()), '; bad link ids', int(p['badlink'].sum()), '; duplicate lines', int(p['dup_lines'].sum()), '; CRC failures', int(p['crc_bad'].sum()), '; hamming-flagged lines', int(p['ham_flagged'].sum()))
for k, v in f.iteritems(recursive=False):
    if k.startswith('info') or k.startswith('merge_note'): print(k, ':', str(v))
