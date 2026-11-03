"""Decoding statistics of one run from the merged ROOT files: the converter summaries (info_seg*), the boundary
de-duplication notes (merge_note_*), the per-RCDAQ-event counters (packets) and the per-trigger completeness.

    python3 run_decode_stats.py RUN    -> /sphenix/user/tmengel/eic/output/sync_runs/run<RUN>/decode.txt, decode.npz
"""
import sys, os, re, glob, numpy as np, uproot
RUN = int(sys.argv[1])
D = '/sphenix/user/tmengel/eic/bic_testbeam/data'
OUT = '/sphenix/user/tmengel/eic/output/sync_runs/run%d' % RUN; os.makedirs(OUT, exist_ok=True)
log = open(OUT + '/decode.txt', 'w'); R = {}
def say(*a):
    s = ' '.join(str(x) for x in a); print(s, flush=True); log.write(s + '\n')
def files(stream):
    f = sorted(glob.glob('%s/%s_%d_part*.root' % (D, stream, RUN)))
    return f or [p for p in ['%s/%s_%d.root' % (D, stream, RUN)] if os.path.exists(p)]
def notes(fs, prefix):
    out = {}
    for f in fs:
        F = uproot.open(f)
        for k, obj in F.items(filter_classname='TNamed', recursive=False):     # items(): names with '/' are not paths
            name = k.split(';')[0]
            if name.startswith(prefix) and not name.endswith('_head'): out[f + name] = obj.member('fTitle')
    return out

KRE = re.compile(r'Triggers kept (\d+) \(\+(\d+) in the head file\), left to the next segment (\d+), waveforms seen (\d+); pool: CRC-flagged lines (\d+), '
                 r'hamming-flagged (\d+), bad lines (\d+), late lines (\d+), resets (\d+); ADC cross-check mismatches (\d+); backward ts0 steps (\d+); '
                 r'waveforms without any line found (\d+)')
KF = ['kept', 'head', 'left', 'seen', 'crc_lines', 'ham_lines', 'bad_lines', 'late_lines', 'resets', 'adc_mismatch', 'backward', 'no_line']
say('==== run %d decoding' % RUN)
for st in ('kcu0', 'kcu1'):
    fs = files(st)
    if not fs: say(st, ': no file'); continue
    info = notes(fs, 'info_seg'); tot = np.zeros(len(KF), np.int64)
    for k, v in info.items():
        m = KRE.search(v)
        if m: tot += np.array([int(x) for x in m.groups()])
    mn = notes(fs, 'merge_note'); dropped = [tuple(map(int, re.search(r'(\d+) of (\d+) triggers dropped', v).groups())) for v in mn.values()]
    say('%s: %d segment summaries; %s' % (st, len(info), ', '.join('%s %d' % (a, b) for a, b in zip(KF, tot))))
    say('  segment boundaries: %d head files de-duplicated, %d of %d head triggers dropped as already complete in the previous segment'
        % (len(dropped), sum(d[0] for d in dropped), sum(d[1] for d in dropped)))
    pk = {k: [] for k in ('evt', 'evt_type', 'has_kcu', 'udp_jumps', 'udp_lost', 'nomarker', 'crc_bad', 'ham_flagged', 'dup_lines', 'badlink', 'nlines')}
    hd = {k: [] for k in ('nsamp', 'nlines', 'complete')}
    for f in fs:
        F = uproot.open(f)
        for b in F['packets'].iterate(list(pk), step_size=2000000, library='np'):
            for k in pk: pk[k].append(b[k])
        for b in F['header'].iterate(list(hd), step_size=5000000, library='np'):
            for k in hd: hd[k].append(b[k])
    if not any(len(x) for x in hd['nsamp']): say('%s: files present but no triggers' % st); continue
    pk = {k: np.concatenate(v) for k, v in pk.items()}; hd = {k: np.concatenate(v) for k, v in hd.items()}
    ev = np.unique(pk['evt'][pk['has_kcu']])
    lines = int(pk['nlines'].sum())
    say('  RCDAQ events %d (with the KCU packet %d, types %s); sequence %d..%d, gaps %d'
        % (len(pk['evt']), pk['has_kcu'].sum(), dict(zip(*[x.tolist() for x in np.unique(pk['evt_type'], return_counts=True)])), ev.min(), ev.max(), (np.diff(ev) != 1).sum()))
    say('  UDP: sequence jumps %d, datagrams lost %d; lines %d: marker-less %d, bad link %d, duplicate %d, CRC failures %d (%.3f %%), Hamming-flagged %d (%.3f %%)'
        % (pk['udp_jumps'].sum(), pk['udp_lost'].sum(), lines, pk['nomarker'].sum(), pk['badlink'].sum(), pk['dup_lines'].sum(),
           pk['crc_bad'].sum(), 100 * pk['crc_bad'].sum() / max(lines, 1), pk['ham_flagged'].sum(), 100 * pk['ham_flagged'].sum() / max(lines, 1)))
    n = len(hd['nsamp'])
    say('  triggers %d: complete (20 samples, 160 lines) %.2f %%, 20 samples %.2f %%, single-sample entries %d (%.2f %%), mean lines %.1f of 160'
        % (n, 100 * hd['complete'].mean(), 100 * (hd['nsamp'] == 20).mean(), (hd['nsamp'] == 1).sum(), 100 * (hd['nsamp'] == 1).mean(), hd['nlines'].mean()))
    R[st] = dict(zip(KF, tot.tolist())) | dict(n=n, complete=float(hd['complete'].mean()), nsamp20=float((hd['nsamp'] == 20).mean()),
                 single=int((hd['nsamp'] == 1).sum()), udp_lost=int(pk['udp_lost'].sum()), udp_jumps=int(pk['udp_jumps'].sum()), lines=lines,
                 crc_bad=int(pk['crc_bad'].sum()), ham=int(pk['ham_flagged'].sum()), nomarker=int(pk['nomarker'].sum()), rcdaq=int(pk['has_kcu'].sum()),
                 evt_gaps=int((np.diff(ev) != 1).sum()), dedup_dropped=sum(d[0] for d in dropped))

fs = files('caen')
if fs:
    info = notes(fs, 'info_seg'); tr = cmp_ = mis = bs = bm = 0
    for v in info.values():
        m = re.search(r'(\d+) triggers; bad size (\d+), bad marker (\d+); stock-decoder cross-check: (\d+) sample comparisons, (\d+) mismatches', v)
        tr += int(m.group(1)); bs += int(m.group(2)); bm += int(m.group(3)); cmp_ += int(m.group(4)); mis += int(m.group(5))
    h = [uproot.open(f)['header'].arrays(['evt', 'event_counter', 'gtt', 'ok_size', 'ok_marker'], library='np') for f in fs]
    h = {k: np.concatenate([x[k] for x in h]) for k in h[0]}
    ec = np.diff(h['event_counter'].astype(np.int64))
    say('caen: %d segment summaries; triggers %d; bad size %d, bad marker %d; stock-decoder cross-check %d samples, %d mismatches'
        % (len(info), tr, bs, bm, cmp_, mis))
    say('  event counter steps %s; gtt backward steps %d; RCDAQ sequence gaps %d'
        % (dict(zip(*[x.tolist() for x in np.unique(ec, return_counts=True)])) if len(ec) else {}, (np.diff(h['gtt']) < 0).sum(), (np.diff(h['evt']) != 1).sum()))
    R['caen'] = dict(n=tr, bad_size=bs, bad_marker=bm, compared=cmp_, mismatches=mis, ec_ok=bool(len(ec) == 0 or (ec == 1).all()), gtt_back=int((np.diff(h['gtt']) < 0).sum()))
else:
    say('caen: no file')

fs = files('astropix')
if fs:
    F = uproot.open(fs[0]); v = F['info'].member('fTitle') if 'info' in F else ''
    m = re.search(r'(\d+) RCDAQ events, (\d+) hit records, (\d+) other records; stock-decoder hit-count cross-check on (\d+) packets, (\d+) mismatches', v)
    hv = F['hits'].arrays(['valid', 'id', 'layer'], library='np'); p = F['packets'].arrays(['truncated', 'nhit'], library='np')
    say('astropix: RCDAQ events %s, hit records %s, other (filler) records %s; cross-check on %s packets, %s mismatches' % m.groups() if m else 'astropix: no info')
    say('  valid hits (layer 1, chip 0-8) %d of %d (%.2f %%); chips %s; truncated packets %d'
        % (hv['valid'].sum(), len(hv['valid']), 100 * hv['valid'].mean() if len(hv['valid']) else 0, np.unique(hv['id'][hv['valid']]).tolist(), int(p['truncated'].sum())))
    R['astropix'] = dict(hits=len(hv['valid']), valid=int(hv['valid'].sum()), other=int(m.group(3)) if m else 0, mismatches=int(m.group(5)) if m else -1)
else:
    say('astropix: no file')
np.savez(OUT + '/decode.npz', R=np.array(R, dtype=object))
say('done')
