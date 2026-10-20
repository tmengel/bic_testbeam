"""ts-phase at P=16 KCU ticks (10 MHz): true pairing vs ts-shuffled control,
full phase profile, and offset (K - stamp) distribution of phase-matched pairs."""
import numpy as np, sys
tk=np.load('out/kcu0_ts.npy').astype(np.int64)
c=np.load('out/astro_clusters.npz')
O=int(c['coarse'])-int(0.43e9/6.25)
sel=c['n']>=2; F=c['Fmin'][sel]+O; TS=c['TS'][sel].astype(np.int64)
P=16; lo,hi=int(-200e6/6.25),int(100e6/6.25)
rng=np.random.default_rng(7); TSs=TS[rng.permutation(len(TS))]
i0=np.searchsorted(tk,F+lo); i1=np.searchsorted(tk,F+hi)
H=np.zeros((2,P,256),np.int64)
for s0 in range(0,len(F),2000):
    sl=slice(s0,s0+2000); cnt=i1[sl]-i0[sl]
    idx=np.concatenate([np.arange(a,b) for a,b in zip(i0[sl],i1[sl])]); K=tk[idx]
    for t,T in enumerate((TS,TSs)):
        rep=np.repeat(T[sl],cnt)
        for s in range(0,P,4):                       # KCU triggers sit on a 4-tick grid
            H[t,s]+=np.bincount(((K+s)//P-rep)%256,minlength=256)
for t,name in enumerate(("true","shuffled ts")):
    h=H[t,::4]; z=(h-h.mean(1,keepdims=True))/np.sqrt(h.mean(1,keepdims=True)); q,ph=np.unravel_index(z.argmax(),z.shape)
    print(f"{name:12s}: max {z.max():5.1f} sigma at phase {ph}, sub-offset {4*q};  min {z.min():5.1f}")
h=H[0,::4]; z=(h-h.mean(1,keepdims=True))/np.sqrt(h.mean(1,keepdims=True)); q,ph=np.unravel_index(z.argmax(),z.shape)
print("profile at sub-offset",4*q,": phase:z ", " ".join(f"{(ph+k)%256}:{z[q,(ph+k)%256]:.1f}" for k in range(-6,7)))
np.save('out/astro_ts10_H.npy',H)
