"""Group AstroPix hits into clusters (one particle shower) by on-chip ts, and get
the coarse AstroPix->KCU0 offset from the minimum-latency envelopes."""
import numpy as np, uproot, glob
a=np.load('out/astro_unwrapped.npz')
raw=uproot.open('out/astro_3062.root')['astro'].arrays(['ts','layer','id'],library='np')
ok=(raw['layer']==1)&(raw['id']<=8)
TS=raw['ts'][ok].astype(np.int64); F=a['tu'].astype(np.int64)*4; rt=a['rt']; ev=a['ev']; ID=a['id']; tot=a['tot']
# clusters: consecutive hits whose ts stays within +-2 of the cluster's first hit, read within 1 ms
start=np.zeros(len(F),bool); start[0]=True; ref=TS[0]; f0=F[0]
for i in range(1,len(F)):
    d=(TS[i]-ref+128)%256-128
    if abs(d)>2 or F[i]-f0>160000: start[i]=True; ref=TS[i]; f0=F[i]
cid=np.cumsum(start)-1; first=np.nonzero(start)[0]
n=np.bincount(cid); nchip=np.array([len(np.unique(ID[s:e])) for s,e in zip(first,np.r_[first[1:],len(F)])])
print("clusters:",len(first)," size distribution:",{k:int((n==k).sum()) for k in range(1,9)},">8:",int((n>8).sum()))
# coarse offset from envelopes (last row of each RCDAQ event = least latency)
def env(t,rt,ev):
    last=np.r_[np.nonzero(np.diff(ev))[0],len(ev)-1]; tt=t[last]; L=rt[last]-tt
    b=((tt-tt[0])//10e9).astype(int); return np.median([np.percentile(L[b==k],1) for k in np.unique(b) if (b==k).sum()>=10])
fs=sorted(glob.glob('out/beam_ROC_0-00003062-*.root'))
A=[uproot.open(f)['kcu'].arrays(['evt','rt','ts','nsamp'],library='np') for f in fs]
k={x:np.concatenate([q[x] for q in A]) for x in A[0]}; m=(k['nsamp']==20)&(k['ts']>1e8)
Lk=env(k['ts'][m]*6.25,k['rt'][m],k['evt'][m]); La=env(F*6.25,rt,ev)
coarse=(La-Lk)/6.25
print("envelope epochs: KCU0 %.6f s, Astro %.6f s -> coarse %+.6f s"%(Lk/1e9,La/1e9,(La-Lk)/1e9))
np.savez('out/astro_clusters.npz',Ffirst=F[first],TS=TS[first],n=n,nchip=nchip,coarse=coarse,
         Fmin=np.array([F[s:e].min() for s,e in zip(first,np.r_[first[1:],len(F)])]))
