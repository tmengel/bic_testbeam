"""AstroPix readout latency: K - readout stamp for phase-matched (187-190) and control (60-63) cluster-trigger pairs
-> out/astro_latency.npz (dm, dx in ms), used by astro_summary_plot.py"""
import numpy as np
tk=np.load('out/kcu0_ts.npy').astype(np.int64)
c=np.load('out/astro_clusters.npz')
O=int(c['coarse'])-int(0.43e9/6.25)
sel=c['n']>=2; F=c['Fmin'][sel]+O; TS=c['TS'][sel].astype(np.int64)
P=16; s=4; lo,hi=int(-400e6/6.25),int(300e6/6.25)
i0=np.searchsorted(tk,F+lo); i1=np.searchsorted(tk,F+hi)
dm,dx=[],[]
for s0 in range(0,len(F),5000):
    sl=slice(s0,s0+5000); cnt=i1[sl]-i0[sl]
    idx=np.concatenate([np.arange(a,b) for a,b in zip(i0[sl],i1[sl])]); K=tk[idx]
    rep=np.repeat(TS[sl],cnt); Fr=np.repeat(F[sl],cnt)
    ph=((K+s)//P-rep)%256
    m=(ph>=187)&(ph<=190); x=(ph>=60)&(ph<=63)          # matched phases vs a control band
    dm.append((K-Fr)[m]); dx.append((K-Fr)[x])
dm=np.concatenate(dm)*6.25e-6; dx=np.concatenate(dx)*6.25e-6     # ms
e=np.arange(-400,301,5); hm=np.histogram(dm,e)[0]; hx=np.histogram(dx,e)[0]
print(" K - stamp [ms]   matched   control   excess")
for k in range(len(hm)):
    ex=hm[k]-hx[k]
    if ex>3*np.sqrt(hx[k]+1): print("  %+5d..%+5d   %7d   %7d   %+6d"%(e[k],e[k+1],hm[k],hx[k],ex))
np.savez('out/astro_latency.npz',dm=dm,dx=dx)
