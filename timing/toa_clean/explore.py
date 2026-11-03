import numpy as np, uproot
f = uproot.open('/sphenix/user/tmengel/eic/bic_testbeam/data/kcu0_3062.root'); NMAX = 1200000
for br in ('Sc3', 'Nc3'):
    a = f['sfil2'].arrays([br + '_adc', br + '_toa'], entry_stop=NMAX, library='np'); W = a[br + '_adc'].astype(float); T = a[br + '_toa'].astype(int)
    n = len(W); good = W[:, 0] >= 0; y = W - W[:, 1:4].mean(axis=1)[:, None]
    ipk = y[:, 4:16].argmax(axis=1) + 4; A = y[np.arange(n), ipk]; k = ipk.copy()
    for d in range(1, 6):
        kk = np.clip(ipk - d, 0, 19); k = np.where(y[np.arange(n), kk] >= 0.5 * A, kk, k)
    sel = good & (A > 150) & (A < 950) & (ipk >= 5) & (ipk <= 14) & (W.max(axis=1) < 1000)
    noped = good & (A < 15)
    print(br, 'pulses', sel.sum(), ' no-pulse triggers', noped.sum())
    print('  word rate per sample, no-pulse triggers:', np.round((T[noped] > 0).mean(axis=0), 3))
    for d in range(-3, 4):
        s = np.clip(k + d, 0, 19); w = T[np.arange(n), s]; print('  pulses: word on k%+d: %.3f' % (d, (w[sel] > 0).mean()))
    B = T[noped][:, 1:19].ravel(); hb = np.bincount(B[B > 0], minlength=1024)[1:1024]; wk = T[np.arange(n), k]; hk = np.bincount(wk[sel & (wk > 0)], minlength=1024)[1:1024]
    wk1 = T[np.arange(n), np.clip(k + 1, 0, 19)]; hk1 = np.bincount(wk1[sel & (wk1 > 0)], minlength=1024)[1:1024]
    print('  both k and k+1: %.3f, k only %.3f, k+1 only %.3f, neither %.3f' % (((wk>0)&(wk1>0))[sel].mean(), ((wk>0)&(wk1<=0))[sel].mean(), ((wk<=0)&(wk1>0))[sel].mean(), ((wk<=0)&(wk1<=0))[sel].mean()))
    print('  negative toa values present:', (T < 0).mean())
    for nm, h in (('no-pulse', hb), ('k', hk), ('k+1', hk1)):
        c = np.add.reduceat(h, np.arange(0, 1023, 64)); print('  %-8s codes per 64:' % nm, ' '.join('%d' % v for v in c))
    np.savez('explore_%s.npz' % br, hb=hb, hk=hk, hk1=hk1)
