"""Run-generic stream synchronization from the merged ROOT trees (kcu_to_root / caen_to_root / astro_to_root output).

    python3 run_sync.py RUN            -> /sphenix/user/tmengel/eic/output/sync_runs/run<RUN>/

Same method as the run-3062 chain (reproduce.sh: kcu1_align, caen_kcu0, trigger_table, astro_clusters, astro_onoff,
astro_ts10, astro_latency), but every constant is measured on the run itself:
  0 KCU triggers: header entries after the stale-clock prefix (the entries before the one flagged backward hold the
    counter from before the reset at run start); an entry with nsamp < 20 sitting k x 164 ticks (k = 1..19, +-4) after
    another entry is the continuation of that trigger (decoder fragment), not a trigger
  1 KCU1 -> KCU0: coarse offset from the host-clock epochs L = rt - 6.25 ns x ts0, exact offset = the most common
    ts(KCU0) - ts(KCU1) within +-20 ms of it; checked in 10 time slices
  2 CAEN -> KCU0: t_CAEN = gtt x 8.533 ns; coarse from the epochs, then the coincidence peak in slices of 300 CAEN
    triggers along the run (anchor: the most significant of 15 slices over the middle 90 % of the run, searched over +-10 s; outward from it, a slice without a peak within
    +-2 ms of the local line is skipped), a line through the slice offsets (clock drift), a running
    median of matched-pair residuals; control = the same with CAEN moved by +10 ms
  3 common trigger table: KCU0, KCU1 (+offset) and CAEN (on the KCU0 clock) grouped within 16 ticks (100 ns)
  4 paddle check: KCU0 pattern (pedestal - min(samples 6-19) > 65) vs CAEN pattern (pedestal - min(samples 20-229) > 200)
    on shared triggers, against the neighbouring trigger as a decoy
  5 AstroPix: hits of layer 1, chips 0-8; clusters; coarse offset from the readout envelopes; beam on/off lag against
    the CAEN rate; 10 MHz phase ((K + s) // 16 - ts) mod 256 against a ts-shuffled control; readout latency of
    phase-matched minus control pairs
Outputs: sync.txt (log), sync.npz (constants + histograms for run_sync_figs.py), triggers.npz (the table, with entry
numbers into the merged files, which are listed in the npz).
"""
import sys, os, glob, numpy as np, uproot
RUN = int(sys.argv[1])
D = '/sphenix/user/tmengel/eic/bic_testbeam/data'
OUT = '/sphenix/user/tmengel/eic/output/sync_runs/run%d' % RUN; os.makedirs(OUT, exist_ok=True)
TICK = 6.25; CT = 1e3 / 117.1875                       # ns per KCU tick, per CAEN GROUPTRIGGERTIME tick
log = open(OUT + '/sync.txt', 'w')
def say(*a):
    s = ' '.join(str(x) for x in a); print(s, flush=True); log.write(s + '\n'); log.flush()
R = {}                                                 # everything that goes into sync.npz

def files(stream):
    f = sorted(glob.glob('%s/%s_%d_part*.root' % (D, stream, RUN)))
    return f or [p for p in ['%s/%s_%d.root' % (D, stream, RUN)] if os.path.exists(p)]

def bits(h): return (h * (1 << np.arange(4))).sum(1).astype(np.int8)

# ---------------------------------------------------------------- 0 KCU triggers
def load_kcu(stream, paddles=False):
    fs = files(stream); cols = {k: [] for k in ('ts0', 'rt', 'ev', 'nsamp', 'complete', 'backward', 'pad')}; n0 = 0
    for f in fs:
        F = uproot.open(f); h = F['header']
        a = h.arrays(['ts0', 'rt_first', 'evt_first', 'nsamp', 'complete', 'backward'], library='np')
        cols['ts0'].append(a['ts0'].astype(np.int64)); cols['rt'].append(a['rt_first'].astype(np.int64)); cols['ev'].append(a['evt_first'])
        cols['nsamp'].append(a['nsamp']); cols['complete'].append(a['complete']); cols['backward'].append(a['backward'])
        if paddles:
            for b in F['paddles'].iterate(['P%d_adc' % i for i in (3, 4, 5, 6)], step_size=500000, library='np'):
                dip = np.stack([np.median(b['P%d_adc' % i][:, :4], 1) - b['P%d_adc' % i][:, 6:20].min(1) for i in (3, 4, 5, 6)], 1)
                cols['pad'].append(bits(dip > 65))
        n0 += h.num_entries
    if n0 == 0: return None
    K = {k: np.concatenate(v) for k, v in cols.items() if v}; K['entry'] = np.arange(n0); K['files'] = fs
    bw = np.flatnonzero(K['backward'])
    stale = bw[0] if len(bw) and bw[0] < 1000 else 0
    say('%s: %d entries in %d file(s); backward-flagged entries %s; stale-clock prefix dropped: %d' % (stream, n0, len(fs), bw[:5].tolist(), stale))
    for k in ('ts0', 'rt', 'ev', 'nsamp', 'complete', 'entry') + (('pad',) if paddles else ()): K[k] = K[k][stale:]
    if (np.diff(K['ts0']) < 0).any():
        say('  WARNING: %d further backward steps, sorting' % (np.diff(K['ts0']) < 0).sum())
        o = np.argsort(K['ts0'], kind='stable')
        for k in ('ts0', 'rt', 'ev', 'nsamp', 'complete', 'entry') + (('pad',) if paddles else ()): K[k] = K[k][o]
    t = K['ts0']; cont = np.zeros(len(t), bool)
    for k in range(1, 20):
        a, b = np.searchsorted(t, t - k * 164 - 4), np.searchsorted(t, t - k * 164 + 5); cont |= b > a
    cont &= K['nsamp'] < 20
    K['start'] = ~cont
    st = K['start']
    say('  triggers %d (complete %.1f %%, nsamp 20 %.1f %%), continuation fragments %d; span %.1f s' % (st.sum(), 100 * K['complete'][st].mean(),
        100 * (K['nsamp'][st] == 20).mean(), cont.sum(), (t[-1] - t[0]) * TICK * 1e-9))
    R[stream + '_ntrig'] = st.sum(); R[stream + '_ncomplete'] = K['complete'][st].sum(); R[stream + '_nfrag'] = cont.sum()
    R[stream + '_span_s'] = (t[-1] - t[0]) * TICK * 1e-9; R[stream + '_stale'] = stale
    return K

say('==== run %d' % RUN)
K0 = load_kcu('kcu0', paddles=True)
K1 = load_kcu('kcu1')
s0 = K0['start']; t0 = K0['ts0'][s0]; e0 = K0['entry'][s0]; p0 = K0['pad'][s0]; c0 = K0['complete'][s0]
L0 = np.median(K0['rt'] - TICK * K0['ts0']); R['L_kcu0_s'] = L0 * 1e-9
say('KCU0 host-clock epoch L0 = %.6f s (p1-p99 spread %.2f ms)' % (L0 * 1e-9, np.ptp(np.percentile(K0['rt'] - TICK * K0['ts0'], [1, 99])) * 1e-6))

# ---------------------------------------------------------------- 1 KCU1 -> KCU0
OFF1 = None
if K1 is not None and K1['start'].sum() > 100:
    s1 = K1['start']; t1 = K1['ts0'][s1]; e1 = K1['entry'][s1]; c1 = K1['complete'][s1]
    L1 = np.median(K1['rt'] - TICK * K1['ts0']); coarse = (L1 - L0) / TICK; R['L_kcu1_s'] = L1 * 1e-9
    W = int(20e6 / TICK)
    def diffs(sub):
        lo = np.searchsorted(t0, sub + int(coarse) - W); hi = np.searchsorted(t0, sub + int(coarse) + W)
        return np.concatenate([t0[a:b] - s for a, b, s in zip(lo, hi, sub)])
    rng = np.random.default_rng(5); sub = np.sort(rng.choice(t1, min(100000, len(t1)), replace=False))
    d = diffs(sub); vals, cnt = np.unique(d, return_counts=True); OFF1 = int(vals[cnt.argmax()])
    bg = np.median(np.bincount(((d - int(coarse) + W) // 1600).astype(np.int64)))    # random pairs per 10 us (same int window as diffs)
    say('KCU1: coarse offset %.0f ticks (%.4f s); exact %d ticks (%.6f s); coarse - exact %.2f ms; peak %d pairs over %.1f random per 10 us'
        % (coarse, coarse * TICK * 1e-9, OFF1, OFF1 * TICK * 1e-9, (coarse - OFF1) * TICK * 1e-6, cnt.max(), bg))
    res = d - OFF1; m = np.abs(res) <= 8
    sh = [np.mean(res[m] == 0), np.mean(res[m] == -4), np.mean(~np.isin(res[m], [0, -4]))]
    say('  residual shares within +-8 ticks: 0: %.4f  -4: %.4f  other: %.4f' % tuple(sh))
    R.update(off1=OFF1, off1_coarse=coarse, off1_shares=sh, off1_res=np.bincount(res[m] + 8, minlength=17),
             off1_spike=np.histogram((d - coarse) * TICK * 1e-6, np.linspace(-20, 20, 4001))[0])
    sl = []
    for part in np.array_split(t1, 10):
        if len(part) < 200: continue
        dd = diffs(np.sort(rng.choice(part, min(20000, len(part)), replace=False))); v, c = np.unique(dd, return_counts=True)
        sl.append(int(v[c.argmax()]))
    say('  exact offset in 10 time slices:', sorted(set(sl)), '(all equal)' if len(set(sl)) == 1 else '(NOT constant)')
    R['off1_slices'] = sl
else:
    say('KCU1: no triggers, skipped'); t1 = np.zeros(0, np.int64); e1 = np.zeros(0, np.int64); c1 = np.zeros(0, bool)

# ---------------------------------------------------------------- 2 CAEN -> KCU0
tk = t0 * TICK                                          # KCU0 trigger starts in ns
fc = files('caen'); have_caen = bool(fc) and sum(uproot.open(f)['header'].num_entries for f in fc) > 100
if have_caen:
    H = [uproot.open(f)['header'].arrays(['gtt', 'rt'], library='np') for f in fc]
    gtt = np.concatenate([h['gtt'] for h in H]).astype(np.int64); rtc = np.concatenate([h['rt'] for h in H]).astype(np.int64)
    padc = []
    for f in fc:
        for b in uproot.open(f)['paddles'].iterate(['P%d_adc' % i for i in (3, 4, 5, 6)], step_size=200000, library='np'):
            neg = np.stack([np.median(b['P%d_adc' % i][:, :20], 1) - b['P%d_adc' % i][:, 20:230].min(1) for i in (3, 4, 5, 6)], 1)
            padc.append(bits(neg > 200))
    padc = np.concatenate(padc); del H
    say('CAEN: %d triggers in %d file(s); gtt backward steps %d' % (len(gtt), len(fc), (np.diff(gtt) < 0).sum()))
    tc = gtt * CT; Lc = np.percentile(rtc - tc, 1); L0p = np.percentile(K0['rt'] - TICK * K0['ts0'], 1)
    x0 = tc + (Lc - L0p)
    def pairs(x, guess, win):
        lo = np.searchsorted(tk, x + guess - win); hi = np.searchsorted(tk, x + guess + win)
        return np.concatenate([tk[a:b] - xi - guess for a, b, xi in zip(lo, hi, x)] or [np.zeros(0)])
    def peak(x, guess, win, binw):
        d = pairs(x, guess, win); h = np.bincount(((d + win) / binw).astype(np.int64), minlength=int(2 * win / binw))[:int(2 * win / binw)]
        k = h.argmax(); return guess + (k + 0.5) * binw - win, h[k], max(np.median(h), 1)
    inK = np.flatnonzero((x0 > tk[0] + 2e9) & (x0 < tk[-1] - 2e9))
    nsl = int(np.clip((tk[-1] - tk[0]) * 1e-9 / 20, 40, 2500)); starts = inK[np.linspace(0, len(inK) - 301, nsl).astype(int)]
    def slc(s): x = x0[s:s + 300]; return x[x < x[0] + 5e9]                   # <= 300 triggers and <= 5 s (clock drift)
    # anchor: the most significant +-10 s peak among 15 slices spread over the middle 90 % of the run, scored as
    # (n - b) / sqrt(b + 1) with at least 50 pairs: the former score n / b picked the empty run start of 3064 (b = 1) and the walk
    # outward then lost the line (2026-10-30)
    best = None
    for ia in np.linspace(0.05 * (len(starts) - 1), 0.95 * (len(starts) - 1), 15).astype(int):
        g_, n_, b_ = peak(slc(starts[ia]), 0.0, 10e9, 1e5)                      # +-10 s, 100 us bins
        g_, n_, b_ = peak(slc(starts[ia]), g_, 2e6, 1e4)                        # +-2 ms, 10 us bins
        sc = (n_ - b_) / np.sqrt(b_ + 1) if n_ >= 50 else -1
        if best is None or sc > best[4]: best = (ia, g_, n_, b_, sc)
    ia, g, n, b, _ = best
    say('  anchor slice %d of %d: offset %.6f s after the epochs, %d pairs in the 10 us peak over %.1f' % (ia, len(starts), g * 1e-9, n, b))
    T, O = [x0[starts[ia]]], [g]
    for k in list(range(ia + 1, len(starts))) + list(range(ia - 1, -1, -1)):  # outward from the anchor
        x = slc(starts[k]); Ta, Oa = np.array(T), np.array(O); near = np.argsort(np.abs(Ta - x[0]))[:20]
        guess = np.polyval(np.polyfit(Ta[near], Oa[near], 1), x[0]) if len(near) >= 5 else Oa[near[0]]
        o, n, b = peak(x, guess, 2e6, 1e4)                                      # a slice without partners is skipped
        if n >= 5 * b and n >= 8: T.append(x[0]); O.append(o)
    o_ = np.argsort(T); T, O = np.array(T)[o_], np.array(O)[o_]; keep = np.ones(len(T), bool)
    for it in range(5):                                                        # line through the slice offsets, outliers out
        p = np.polyfit(T[keep] - tk[0], O[keep], 1); dev = np.abs(O - np.polyval(p, T - tk[0]))
        keep = dev < max(50e3, 5 * np.median(dev[keep]))
    say('  slice offsets kept in the line: %d of %d' % (keep.sum(), len(T))); T, O = T[keep], O[keep]
    p = np.polyfit(T - tk[0], O, 1)
    say('  slices with a peak: %d of %d; drift %+.3f ppm; KCU ticks per CAEN tick %.9f' % (len(T), nsl, p[0] * 1e6, (1 + p[0]) * CT / TICK))
    sres = (O - np.polyval(p, T - tk[0])) * 1e-3
    say('  slice offset - line [us] p05/p50/p95:', np.round(np.percentile(sres, [5, 50, 95]), 2), ' max |.| %.1f' % np.abs(sres).max())
    def nearest(x):
        i = np.clip(np.searchsorted(tk, x), 1, len(tk) - 1); l, r = tk[i - 1], tk[i]; j = np.where(x - l < r - x, i - 1, i)
        return tk[j] - x, j
    pred = x0 + np.polyval(p, x0 - tk[0]); inspan = (pred > tk[0]) & (pred < tk[-1])
    r, j = nearest(pred); idx = np.flatnonzero(inspan & (np.abs(r) < 50e3))
    knots = idx[::50]; rm = np.array([np.median(r[idx[max(0, k - 100):k + 101]]) for k in range(0, len(idx), 50)])
    corr = np.interp(pred, pred[knots], rm); tpred = pred + corr
    r2, j2 = nearest(tpred)
    for name, shf in (('matched', 0.0), ('shifted +10 ms', 10e6)):
        rr, _ = nearest(tpred + shf)
        fr = [100 * (np.abs(rr[inspan]) < v).mean() for v in (1e3, 100, 50e3)]
        say('  %-15s: of %d CAEN in the KCU0 span, |dt| < 1 us %.1f %%, < 100 ns %.1f %%, < 50 us %.1f %%' % ((name, inspan.sum()) + tuple(fr)))
        R['caen_frac_' + name.split()[0]] = fr
    good = inspan & (np.abs(r2) < 1e3)
    mad = np.median(np.abs(r2[good] - np.median(r2[good])))
    say('  matched residual [ns] p05/p50/p95:', np.round(np.percentile(r2[good], [5, 50, 95]), 1), ' MAD %.1f ns; running-median correction p01/p99 %.0f/%.0f ns'
        % (mad, *np.percentile(corr[inspan], [1, 99])))
    rr, _ = nearest(tpred + 10e6)
    R.update(caen_drift_ppm=p[0] * 1e6, caen_line=p, caen_slice_t=(T - tk[0]) * 1e-9, caen_slice_res_us=sres, caen_mad_ns=mad,
             caen_res_hist=np.histogram(r2[inspan] * 1e-3, np.linspace(-5, 5, 201))[0], caen_ctl_hist=np.histogram(rr[inspan] * 1e-3, np.linspace(-5, 5, 201))[0],
             caen_zoom_hist=np.histogram(r2[good], np.linspace(-50, 50, 101))[0], caen_n=len(gtt), caen_inspan=inspan.sum(), caen_good=good.sum())
    ic = np.flatnonzero(inspan); tcK = np.round(tpred[ic] / TICK).astype(np.int64); pC = padc[ic]
else:
    say('CAEN: no triggers, skipped'); ic = np.zeros(0, np.int64); tcK = np.zeros(0, np.int64); pC = np.zeros(0, np.int8)

# ---------------------------------------------------------------- 3 common trigger table
TOL = 16
t1K = t1 + (OFF1 or 0)
t = np.concatenate([t0, t1K, tcK]); src = np.concatenate([np.zeros(len(t0), np.int8), np.ones(len(t1K), np.int8), np.full(len(tcK), 2, np.int8)])
row = np.concatenate([np.arange(len(t0)), np.arange(len(t1K)), np.arange(len(tcK))])
o = np.argsort(t, kind='stable'); t, src, row = t[o], src[o], row[o]; del o
gid = np.cumsum(np.r_[True, np.diff(t) > TOL]) - 1; ng = gid[-1] + 1
i0 = np.full(ng, -1, np.int64); i1 = np.full(ng, -1, np.int64); iC = np.full(ng, -1, np.int64); dup = 0
for sidx, arr in ((0, i0), (1, i1), (2, iC)):
    k = src == sidx; dup += k.sum() - len(np.unique(gid[k])); arr[gid[k][::-1]] = row[k][::-1]
del t, src, row, gid
TT = np.where(i0 >= 0, t0[np.maximum(i0, 0)], np.where(i1 >= 0, t1K[np.maximum(i1, 0)] if len(t1K) else 0, tcK[np.maximum(iC, 0)] if len(tcK) else 0))
h0, h1, hC = i0 >= 0, i1 >= 0, iC >= 0
say('trigger table: %d distinct triggers (same stream twice in one 100 ns group: %d)' % (ng, dup))
say('  stream combination        triggers   share')
comb = {}
for name, m in (('KCU0 + KCU1 + CAEN', h0 & h1 & hC), ('KCU0 + CAEN', h0 & ~h1 & hC), ('KCU0 + KCU1', h0 & h1 & ~hC), ('KCU0 only', h0 & ~h1 & ~hC),
                ('KCU1 + CAEN', ~h0 & h1 & hC), ('KCU1 only', ~h0 & h1 & ~hC), ('CAEN only', ~h0 & ~h1 & hC)):
    say('  %-22s %10d   %.4f' % (name, m.sum(), m.mean())); comb[name] = (int(m.sum()), float(m.mean()))
R['comb'] = np.array(list(comb.items()), dtype=object); R['ntrig'] = ng
span = (TT[-1] - TT[0]) * TICK * 1e-9
say('  rates over %.0f s: union %.0f Hz; KCU0 %.0f, KCU1 %.0f, CAEN %.0f Hz' % (span, ng / span, h0.sum() / span, h1.sum() / span, hC.sum() / span))
if h0.any() and h1.any():
    d10 = t0[i0[h0 & h1]] - t1K[i1[h0 & h1]]
    say('  KCU0 - KCU1 in shared triggers [ticks]: 0: %.4f  -4: %.4f  other: %.4f' % ((d10 == 0).mean(), (d10 == -4).mean(), (~np.isin(d10, [0, -4])).mean()))

# ---------------------------------------------------------------- 4 paddle patterns
pad0 = np.where(h0, p0[np.maximum(i0, 0)], -1); padC = np.where(hC, pC[np.maximum(iC, 0)] if len(pC) else -1, -1)
both = h0 & hC
if both.sum() > 100:
    same = (pad0[both] == padC[both]).mean(); dec = (pad0[both] == np.roll(padC[both], 1)).mean()
    say('paddle pattern KCU0 vs CAEN on %d shared triggers: identical %.4f (decoy, neighbouring trigger: %.4f)' % (both.sum(), same, dec))
    R.update(pad_same=same, pad_decoy=dec)
for name, pp in (('KCU0', pad0[h0]), ('CAEN', padC[hC])):
    if len(pp) == 0: continue
    nb = np.array([bin(x).count('1') for x in range(16)])[pp]
    fr = [float((pp >> b & 1).mean()) for b in range(4)]
    say('  %s: any paddle %.4f; P3..P6 %s' % (name, (nb > 0).mean(), np.round(fr, 3))); R['pad_any_' + name] = (nb > 0).mean(); R['pad_frac_' + name] = fr
pad = np.where(pad0 >= 0, pad0, padC)
np.savez(OUT + '/triggers.npz', T=TT, e_kcu0=np.where(h0, e0[np.maximum(i0, 0)], -1), e_kcu1=np.where(h1, e1[np.maximum(i1, 0)] if len(e1) else -1, -1),
         e_caen=np.where(hC, ic[np.maximum(iC, 0)] if len(ic) else -1, -1), complete_kcu0=np.where(h0, c0[np.maximum(i0, 0)], False),
         complete_kcu1=np.where(h1, c1[np.maximum(i1, 0)] if len(c1) else False, False), pad=pad, pad_kcu0=pad0, pad_caen=padC,
         files_kcu0=np.array(K0['files']), files_kcu1=np.array(K1['files'] if K1 else []), files_caen=np.array(fc), off1=OFF1 or 0)
say('wrote triggers.npz (T in KCU0 ticks; e_* = entry in the chained merged files, -1 = stream missed the trigger)')

# ---------------------------------------------------------------- 5 AstroPix
fa = files('astropix')
nh = sum(uproot.open(f)['hits'].num_entries for f in fa) if fa else 0
if nh > 1000:
    a = uproot.open(fa[0])['hits'].arrays(['ts', 'fpga_ts', 'rt', 'evt', 'layer', 'id'], library='np')
    ok = (a['layer'] == 1) & (a['id'] <= 8)
    ts = a['fpga_ts'][ok].astype(np.int64); dd = np.diff(ts)
    F = (ts + np.cumsum(np.r_[0, (dd < -2**31).astype(np.int64) - (dd > 2**31)]) * 2**32) * 4      # 40 MHz -> KCU ticks
    TS = a['ts'][ok].astype(np.int64); rta = a['rt'][ok].astype(np.int64); ev = a['evt'][ok]; ID = a['id'][ok]
    say('AstroPix: %d hits, %d on layer 1 chips 0-8 (chips %s)' % (nh, ok.sum(), np.unique(ID).tolist()))
    start = np.zeros(len(F), bool); start[0] = True; ref = TS[0]; f0 = F[0]
    for i in range(1, len(F)):
        if abs((TS[i] - ref + 128) % 256 - 128) > 2 or F[i] - f0 > 160000: start[i] = True; ref = TS[i]; f0 = F[i]
    first = np.flatnonzero(start); ncl = np.diff(np.r_[first, len(F)])
    Fmin = np.minimum.reduceat(F, first); TSc = TS[first]
    say('  clusters %d; size 1/2/3+: %d/%d/%d' % (len(first), (ncl == 1).sum(), (ncl == 2).sum(), (ncl >= 3).sum()))
    def env(tt, rt, e):
        last = np.r_[np.flatnonzero(np.diff(e)), len(e) - 1]; tt = tt[last]; L = rt[last] - tt
        b = ((tt - tt[0]) // 10e9).astype(int); v = [np.percentile(L[b == k], 1) for k in np.unique(b) if (b == k).sum() >= 10]
        return np.median(v) if v else np.percentile(L, 1)
    Lk = env(K0['ts0'] * TICK, K0['rt'], K0['ev']); La = env(F * TICK, rta, ev); coarse = (La - Lk) / TICK
    say('  envelope epochs: KCU0 %.6f s, AstroPix %.6f s -> coarse %+.6f s' % (Lk * 1e-9, La * 1e-9, (La - Lk) * 1e-9))
    # wide scan: time offset (10 ms bins, +-60 s around the envelope offset) x 10 MHz phase, on a subset of clusters;
    # true pairs pile up in one (offset, phase) cell, random pairs are flat in phase
    sel = ncl >= 2 if (ncl >= 2).sum() > 500 else np.ones(len(ncl), bool)
    rng = np.random.default_rng(7); P = 16
    rate = len(t0) / ((t0[-1] - t0[0]) * TICK * 1e-9); W = int(60e9 / TICK); NB = 12000
    nsub = int(min(max(3e8 / (rate * 120), 300), sel.sum()))
    sub = np.sort(rng.choice(np.flatnonzero(sel), nsub, replace=False)); Fw = Fmin[sub] + int(coarse); Tw = TSc[sub]
    b0, b1 = np.searchsorted(t0, Fw - W), np.searchsorted(t0, Fw + W)
    Hw = np.zeros((4, NB, 256), np.int64); step = max(1, int(4e6 / max(rate * 120, 1)))
    for s_ in range(0, nsub, step):
        sl = slice(s_, s_ + step); cnt = b1[sl] - b0[sl]
        if cnt.sum() == 0: continue
        Kx = t0[np.concatenate([np.arange(p_, q_) for p_, q_ in zip(b0[sl], b1[sl])])]
        rep = np.repeat(Tw[sl], cnt); dtb = np.clip((Kx - np.repeat(Fw[sl], cnt) + W) * NB // (2 * W), 0, NB - 1)
        for q in range(4): Hw[q] += np.bincount(dtb * 256 + ((Kx + 4 * q) // P - rep) % 256, minlength=NB * 256).reshape(NB, 256)
    mu = Hw.mean(2, keepdims=True); Zw = np.where(mu > 0, (Hw - mu) / np.sqrt(mu + 1e-9), 0)
    qw, bw_, pw = np.unravel_index(Zw.argmax(), Zw.shape); dt_pk = ((bw_ + 0.5) * 2 * W / NB - W)
    say('  wide scan (%d clusters, +-60 s x 256 phases, 10 ms bins): max %.1f sigma at %+.3f s from the envelope offset, phase %d, sub-offset %d; '
        'z quantiles 0.5/0.9999 of all cells %.1f/%.1f' % (nsub, Zw.max(), dt_pk * TICK * 1e-9, pw, 4 * qw, *np.percentile(Zw[qw], [50, 99.99])))
    R.update(astro_wide_z=Zw.max(), astro_wide_dt_s=dt_pk * TICK * 1e-9, astro_wide_nsub=nsub)
    cands = [('wide scan', int(coarse) + int(dt_pk))]; del Hw, Zw
    # beam on/off (the run-3062 method, astro_onoff.py): CAEN on the KCU0 clock marks beam-on (it has no readout loss);
    # its 10 ms on/off pattern cross-correlated with the AstroPix cluster rate gives the lag. Works where the wide scan
    # drowns: at kHz trigger rates only a few % of clusters have a KCU0 partner.
    if have_caen and len(tcK) > 10000:
        bwt = int(10e6 / TICK); e = np.arange(t0[0] - int(80e9 / TICK), t0[-1] + int(80e9 / TICK), bwt)
        rc = np.histogram(tcK, e)[0].astype(float); ra = np.histogram(Fmin + int(coarse), e)[0].astype(float)
        sm = np.convolve(rc, np.ones(10) / 10, 'same'); on = (sm > 0.1 * np.median(sm[sm > 0])).astype(float)
        a_, b_ = on - on.mean(), ra - ra.mean(); n_ = len(a_)
        xc = np.fft.irfft(np.conj(np.fft.rfft(a_, 2 * n_)) * np.fft.rfft(b_, 2 * n_), 2 * n_)
        xc = np.r_[xc[-n_ + 1:], xc[:n_]] / np.sqrt((a_ * a_).sum() * (b_ * b_).sum() + 1e-9); lag = np.arange(-n_ + 1, n_)
        w_ = np.abs(lag) <= 6000; xs, ls = xc[w_], lag[w_]; k_ = xs.argmax()
        say('  beam on/off: AstroPix cluster rate vs CAEN beam-on, max r %.3f at %+.2f s from the envelope offset (median |r| %.3f, beam-on fraction %.2f)'
            % (xs[k_], -ls[k_] * 0.01, np.median(np.abs(xs)), on.mean()))
        R.update(astro_onoff_r=xs[k_], astro_onoff_dt_s=-ls[k_] * 0.01)
        cands.append(('beam on/off', int(coarse) - int(ls[k_]) * bwt))
    TSs = TSc[sel]; lo, hi = int(-250e6 / TICK), int(250e6 / TICK)
    rng = np.random.default_rng(7); TSh = TSs[rng.permutation(len(TSs))]
    def phase_scan(O):
        Fc = Fmin[sel] + O; a0, a1 = np.searchsorted(t0, Fc + lo), np.searchsorted(t0, Fc + hi)
        Hh = np.zeros((2, 4, 256), np.int64)
        for s_ in range(0, len(Fc), 2000):
            sl = slice(s_, s_ + 2000); cnt = a1[sl] - a0[sl]
            if cnt.sum() == 0: continue
            Kx = t0[np.concatenate([np.arange(p_, q_) for p_, q_ in zip(a0[sl], a1[sl])])]
            for ti, Tz in enumerate((TSs, TSh)):
                rep = np.repeat(Tz[sl], cnt)
                for q in range(4): Hh[ti, q] += np.bincount(((Kx + 4 * q) // P - rep) % 256, minlength=256)
        return Fc, a0, a1, (Hh - Hh.mean(2, keepdims=True)) / np.sqrt(Hh.mean(2, keepdims=True) + 1e-9)
    best = None
    for name, O_ in cands:
        Fc_, a0_, a1_, Z_ = phase_scan(O_)
        q_, ph_ = np.unravel_index(Z_[0].argmax(), Z_[0].shape)
        say('  10 MHz phase scan around the %s offset (%d clusters, +-250 ms): max %.1f sigma at phase %d, sub-offset %d; ts-shuffled control max %.1f sigma'
            % (name, len(Fc_), Z_[0].max(), ph_, 4 * q_, Z_[1].max()))
        if best is None or Z_[0].max() > best[4][0].max(): best = (name, O_, Fc_, (a0_, a1_), Z_)
    src, O, Fc, (a0, a1), Z = best
    q, ph = np.unravel_index(Z[0].argmax(), Z[0].shape); zs = Z[1].max()
    say('  used: %s offset; profile: ' % src + ' '.join('%d:%.1f' % ((ph + k) % 256, Z[0, q, (ph + k) % 256]) for k in range(-4, 5)))
    R.update(astro_z=Z[0, q], astro_zmax=Z[0].max(), astro_phase=ph, astro_sub=4 * q, astro_zctl=zs, astro_coarse=coarse, astro_O=O, astro_ncl=len(first),
             astro_nsel=len(Fc), astro_src=src)
    if Z[0].max() > 6:
        dm, dx = [], []; cph = (ph + 128) % 256
        for s_ in range(0, len(Fc), 5000):
            sl = slice(s_, s_ + 5000); cnt = a1[sl] - a0[sl]
            if cnt.sum() == 0: continue
            ii = np.concatenate([np.arange(p_, q_) for p_, q_ in zip(a0[sl], a1[sl])]); Kx = t0[ii]
            rep = np.repeat(TSs[sl], cnt); Fr = np.repeat(Fc[sl], cnt); phx = ((Kx + 4 * q) // P - rep) % 256
            dm.append((Kx - Fr)[np.abs((phx - ph + 128) % 256 - 128) <= 1]); dx.append((Kx - Fr)[np.abs((phx - cph + 128) % 256 - 128) <= 1])
        dm = np.concatenate(dm) * TICK * 1e-6; dx = np.concatenate(dx) * TICK * 1e-6
        e = np.arange(-250, 251, 5); hm = np.histogram(dm, e)[0]; hx = np.histogram(dx, e)[0]
        say('  latency K - (AstroPix stamp + O), O = envelope offset + the used coarse offset [ms], phase-matched vs control (phase + 128), bins with a > 3 sigma excess:')
        for k in range(len(hm)):
            if hm[k] - hx[k] > 3 * np.sqrt(hx[k] + 1): say('    %+5d..%+5d   %7d   %7d   %+6d' % (e[k], e[k + 1], hm[k], hx[k], hm[k] - hx[k]))
        R.update(astro_lat_m=hm, astro_lat_x=hx, astro_lat_edges=e)
    else:
        say('  no phase lock above 6 sigma: latency not measured')
else:
    say('AstroPix: %d hits, skipped' % nh)
np.savez(OUT + '/sync.npz', **R)
say('done')
