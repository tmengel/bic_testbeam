"""Quick per-channel QA of a short KCU conversion (kcu_to_root ... max_packets=300) of any run: per channel, the fraction of
triggers with a pulse (> 40 ADC above the samples 1-3 pedestal), its median size, peak sample and clipped fraction, and the
number of SFIL1-4 layers with a pulse per trigger. Used for the 2026-10-28 survey of the new runs (docs/new_runs.md).
    python3 newrun_qa.py sample.root [...]"""
import uproot, numpy as np, sys
def channels(f):
    u = uproot.open(f); m = u['map'].arrays(library='np'); n = u['header'].num_entries
    W = np.full((n, 288, 20), -1, np.int16); cache = {}
    for i in range(len(m['ch'])):
        tr, br, ch = m['tree'][i], m['branch'][i], m['ch'][i]
        if '[' in br:                                   # unused / bbcal: adc[j] rows with a ch list
            if tr not in cache: cache[tr] = u[tr].arrays(['ch', 'adc'], library='np')
            j = int(br.split('[')[1][:-1]); W[:, ch] = cache[tr]['adc'][:, j]
        else:
            W[:, ch] = u[tr][br].array(library='np')
    return m, W
for f in sys.argv[1:]:
    m, W = channels(f); good = (W[:, :, 0] >= 0); y = W.astype(float) - W[:, :, 1:4].mean(axis=2, keepdims=True)
    A = y[:, :, 4:16].max(axis=2); pk = y[:, :, 4:16].argmax(axis=2) + 4; hit = good & (A > 40)
    comp = good.all(axis=1).mean()
    print('\n== %s: %d triggers, all 288 channels present in %.2f' % (f, len(W), comp))
    lab = {int(c): (s, int(l), int(ce), e) for c, s, l, ce, e in zip(m['ch'], m['subsystem'], m['layer'], m['cell'], m['end'])}
    rows = []
    for c in range(288):
        s, l, ce, e = lab[c]
        if s == 'unused' and hit[:, c].mean() < 0.01: continue
        h = hit[:, c]; rows.append((s, l, ce, e, c, good[:, c].mean(), h.mean(), np.median(A[h, c]) if h.sum() > 5 else np.nan, np.median(pk[h, c]) if h.sum() > 5 else np.nan, (W[:, c, :] >= 1023).any(axis=1)[good[:, c]].mean()))
    print('  %-8s %3s %4s %3s %4s  present  hit>40  median A  peak smp  clipped' % ('subsys', 'lay', 'cell', 'end', 'ch'))
    for r in sorted(rows, key=lambda r: (r[0], r[1], r[2], r[3])):
        print('  %-8s %3d %4d %3s %4d  %6.2f  %6.3f  %8.0f  %8.1f  %6.3f' % r)
    # layers hit per trigger (SFIL1-4, any cell/end), a track-like signature for cosmics
    sf = np.array([lab[c][0].startswith('SFIL') for c in range(288)]); ly = np.array([lab[c][1] for c in range(288)])
    nl = np.array([[hit[i, sf & (ly == L)].any() for L in (1, 2, 3, 4)] for i in range(len(W))]).sum(axis=1)
    print('  SFIL1-4 layers with a pulse per trigger: ' + ', '.join('%d: %.3f' % (k, (nl == k).mean()) for k in range(5)))
