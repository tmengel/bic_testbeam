"""Per-cluster attributes + signal/control KCU0 candidates over a broad latency
range, for studying what predicts the readout latency."""
import numpy as np, uproot
a=np.load('out/astro_unwrapped.npz')
raw=uproot.open('out/astro_3062.root')['astro'].arrays(['ts','layer','id','loc','iscol','tot'],library='np')
ok=(raw['layer']==1)&(raw['id']<=8)
TSh=raw['ts'][ok].astype(np.int64); Fh=a['tu'].astype(np.int64)*4; ev=a['ev']
ID=raw['id'][ok]; LOC=raw['loc'][ok]; COL=raw['iscol'][ok]; TOT=raw['tot'][ok]
start=np.zeros(len(Fh),bool); start[0]=True; ref=TSh[0]; f0=Fh[0]
for i in range(1,len(Fh)):
    d=(TSh[i]-ref+128)%256-128
    if abs(d)>2 or Fh[i]-f0>160000: start[i]=True; ref=TSh[i]; f0=Fh[i]
fi=np.nonzero(start)[0]
c=np.load('out/astro_clusters.npz'); assert len(fi)==len(c['n'])
# readout windows = RCDAQ events
wf=np.r_[0,np.nonzero(np.diff(ev))[0]+1]; wl=np.r_[wf[1:]-1,len(Fh)-1]
ws,we=Fh[wf],Fh[wl]; gap_before=np.r_[np.nan,(ws[1:]-we[:-1])*6.25e-6]
wi=np.searchsorted(wf,fi,side='right')-1
A=dict(u=(Fh[fi]-ws[wi])*6.25e-6, to_end=(we[wi]-Fh[fi])*6.25e-6, nbefore=fi-wf[wi],
       gap=gap_before[wi], span=(we[wi]-ws[wi])*6.25e-6,
       recent=(fi-np.searchsorted(Fh,Fh[fi]-int(25e6/6.25)))*1.0,
       spacing=np.r_[np.nan,np.diff(Fh)][fi]/5700.0, tot=TOT[fi]*1.0, loc=LOC[fi]*1.0, iscol=COL[fi]*1.0,
       chip=ID[fi]*1.0, n=c['n']*1.0, ts=c['TS']*1.0)
np.savez('out/astro_attr.npz',**A)
tk=np.load('out/kcu0_ts.npy').astype(np.int64); O=int(c['coarse'])-int(0.43e9/6.25)
F=c['Fmin']+O; TS=c['TS'].astype(np.int64)
keep=np.nonzero((c['n']>=2)|np.isin(ID[fi],[3,5,6,8]))[0]
i0=np.searchsorted(tk,F[keep]+int(-110e6/6.25)); i1=np.searchsorted(tk,F[keep]+int(-3e6/6.25))
rows=[]
for s0 in range(0,len(keep),50000):
    sl=slice(s0,s0+50000); n=i1[sl]-i0[sl]
    idx=np.concatenate([np.arange(p,q) for p,q in zip(i0[sl],i1[sl])]); cl=np.repeat(keep[sl],n)
    ph=((tk[idx]+4)//16-TS[cl])%256
    band=np.full(len(ph),-1,np.int8); band[np.abs(ph-189)<=1]=0; band[np.abs(ph-60)<=1]=1; band[np.abs(ph-120)<=1]=2
    k=band>=0
    rows.append(np.c_[cl[k].astype(np.int32),((tk[idx[k]]-F[cl[k]])*6.25e-6).astype(np.float32),band[k]])
R=np.concatenate(rows)
np.savez('out/astro_cand_bands.npz',cl=R[:,0].astype(np.int32),lat=R[:,1].astype(np.float32),band=R[:,2].astype(np.int8))
print("clusters kept",len(keep),"rows",len(R))
