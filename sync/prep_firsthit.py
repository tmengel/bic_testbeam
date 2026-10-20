"""First-hit index of every AstroPix cluster (same clustering as astro_clusters.py)
-> out/astro_cluster_firsthit.npy. Needed by astro_geo_common.py."""
import numpy as np, uproot
a=np.load('out/astro_unwrapped.npz')
raw=uproot.open('out/astro_3062.root')['astro'].arrays(['ts','layer','id'],library='np')
ok=(raw['layer']==1)&(raw['id']<=8)
TS=raw['ts'][ok].astype(np.int64); F=a['tu'].astype(np.int64)*4
start=np.zeros(len(F),bool); start[0]=True; ref=TS[0]; f0=F[0]
for i in range(1,len(F)):
    d=(TS[i]-ref+128)%256-128
    if abs(d)>2 or F[i]-f0>160000: start[i]=True; ref=TS[i]; f0=F[i]
fi=np.nonzero(start)[0]
assert len(fi)==len(np.load('out/astro_clusters.npz')['n'])
np.save('out/astro_cluster_firsthit.npy',fi); print("clusters:",len(fi))
