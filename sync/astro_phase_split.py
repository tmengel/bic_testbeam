"""AstroPix on-chip ts phase lock to KCU0, split by chip, record type, latency population,
ToT and run time. Same model as sync/astro_ts10.py: phase = ((K+4)//16 - ts) mod 256 over all
complete KCU0 triggers K in the two latency cores before the readout stamp. Also the
background-subtracted latency distribution per chip for phase-matched triggers.
Run from sync/ after the dump, prep and astro stages. Results of 2026-10-21 are in docs/method.md 4.3."""
import numpy as np, sys
tk = np.load('out/kcu0_ts.npy').astype(np.int64)
c = np.load('out/astro_clusters.npz')
a = np.load('out/astro_unwrapped.npz')
O = int(c['coarse']) - int(0.43e9 / 6.25)
F = c['Fmin'] + O; TS = c['TS'].astype(np.int64); n = c['n']
chip = np.load('out/astro_cluster_chip.npy'); fh = np.load('out/astro_cluster_firsthit.npy')
iscol = a['iscol'][fh]; tot = a['tot'][fh]
tsec = (F * 6.25e-9); tbin = np.clip((tsec // 200).astype(int), 0, 7)
totbin = np.digitize(tot, [1, 300, 800, 1500])          # 0: tot==0, 1: 1-299, 2: 300-799, 3: 800-1499, 4: >=1500
N = len(F)
# keys: chip(9) x iscol(2) x pop(2) x tbin(8) x totbin(5) x nbin(2)
dims = (9, 2, 2, 8, 5, 2)
H = np.zeros(dims + (256,), np.int64)
# latency histogram per chip, 1 ms bins from -120..0 ms, phase-hit (|d|<=1) vs all
LB = 120
L = np.zeros((9, 2, LB), np.int64)
ph0 = 189
for pop, (lo, hi) in enumerate(((-30, -5), (-100, -75))):
    i0 = np.searchsorted(tk, F + int(lo * 1e6 / 6.25)); i1 = np.searchsorted(tk, F + int(hi * 1e6 / 6.25))
    for s0 in range(0, N, 50000):
        sl = slice(s0, min(s0 + 50000, N)); cnt = i1[sl] - i0[sl]
        idx = np.concatenate([np.arange(x, y) for x, y in zip(i0[sl], i1[sl])])
        cl = np.repeat(np.arange(sl.start, sl.stop), cnt)
        K = tk[idx]
        ph = ((K + 4) // 16 - TS[cl]) % 256
        key = np.ravel_multi_index((chip[cl], iscol[cl], np.full(len(cl), pop), tbin[cl], totbin[cl], (n[cl] >= 2).astype(int), ph), dims + (256,))
        H += np.bincount(key, minlength=H.size).reshape(H.shape)
        d = (ph - ph0 + 128) % 256 - 128
        lat = ((F[cl] - K) * 6.25e-6).astype(int)   # ms
        ok = (lat >= 0) & (lat < LB)
        hit = (np.abs(d) <= 1).astype(int)
        L += np.bincount(np.ravel_multi_index((chip[cl][ok], hit[ok], lat[ok]), (9, 2, LB)), minlength=L.size).reshape(L.shape)
np.savez(sys.argv[1] if len(sys.argv) > 1 else 'astro_phase_split.npz', H=H, L=L)

def report(name, h):
    h = h.astype(float)
    if h.sum() < 1000: print(f"  {name:28s} too few ({int(h.sum())})"); return
    bg = np.median(h); sd = np.sqrt(max(bg, 1)); z = (h - bg) / sd; pk = int(z.argmax())
    tri = (h[(pk - 1) % 256] + h[pk] + h[(pk + 1) % 256] - 3 * bg) / h.sum()
    print(f"  {name:28s} N={int(h.sum()):>10d}  peak {pk:3d}  z={z[pk]:6.1f}  neighbours z={z[(pk-1)%256]:5.1f},{z[(pk+1)%256]:5.1f}  excess/cand={tri:.4f}")

print("ALL"); report("all", H.sum(axis=(0, 1, 2, 3, 4, 5)))
print("BY CHIP"); [report(f"chip {k}", H[k].sum(axis=(0, 1, 2, 3, 4))) for k in range(9)]
print("BY RECORD TYPE"); [report(f"{['row', 'column'][k]} first record", H[:, k].sum(axis=(0, 1, 2, 3, 4))) for k in range(2)]
print("BY LATENCY POPULATION"); [report(f"pop {k} ({['5-30', '75-100'][k]} ms)", H[:, :, k].sum(axis=(0, 1, 2, 3, 4))) for k in range(2)]
print("BY RUN TIME (200 s bins)"); [report(f"t {k*200}-{k*200+200} s", H[:, :, :, k].sum(axis=(0, 1, 2, 3, 4))) for k in range(8)]
print("BY ToT"); [report(f"tot bin {k} {['0', '1-299', '300-799', '800-1499', '>=1500'][k]}", H[:, :, :, :, k].sum(axis=(0, 1, 2, 3, 4))) for k in range(5)]
print("BY CLUSTER SIZE"); [report(f"n {['==1', '>=2'][k]}", H[:, :, :, :, :, k].sum(axis=(0, 1, 2, 3, 4))) for k in range(2)]
print("BY CHIP x RECORD TYPE");
for k in range(9):
    for r in range(2): report(f"chip {k} {['row', 'col'][r]}", H[k, r].sum(axis=(0, 1, 2, 3)))
print("LATENCY per chip: background-subtracted (phase-hit - 3/256 * all) in ms, peaks")
for k in range(9):
    hit = L[k, 1].astype(float); allc = L[k].sum(0).astype(float)
    ex = hit - allc * 3 / 256
    if allc.sum() < 1000: continue
    top = np.argsort(ex)[::-1][:6]
    print(f"  chip {k}: excess total {ex.sum():8.0f}; top ms bins {sorted(top.tolist())}; "
          f"sum 12-26 ms {ex[12:27].sum():7.0f}, sum 80-96 ms {ex[80:97].sum():7.0f}, elsewhere {ex.sum()-ex[12:27].sum()-ex[80:97].sum():7.0f}")
