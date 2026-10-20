"""Which KCU0 observables distinguish true AstroPix partners (phase-matched excess)
from random triggers (control phases)?"""
import numpy as np, uproot, glob
fs=sorted(glob.glob('out/adc_beam_ROC_0-00003062-*.root'))
T=[uproot.open(f)['kcuadc'].arrays(['ts','amp'],library='np') for f in fs]
ts=np.concatenate([t['ts'] for t in T]).astype(np.int64); amp=np.concatenate([t['amp'] for t in T]).astype(np.float32)
o=np.argsort(ts,kind='stable'); ts,amp=ts[o],amp[o]
u=np.r_[True,np.diff(ts)>0]; ts,amp=ts[u],amp[u]
np.savez('out/kcu0_amp.npz',ts=ts,amp=amp)
print("KCU0 complete triggers with amplitudes:",len(ts))
c=np.load('out/astro_clusters.npz'); O=int(c['coarse'])-int(0.43e9/6.25)
A=np.load('out/astro_attr.npz'); C=np.load('out/astro_cand_bands.npz')
cl,lat,band=C['cl'],C['lat'].astype(np.float64),C['band']
toend=A['to_end'][cl]; chip=A['chip'][cl]
inwin=((lat>-24)&(lat<-12)&(toend<20))|((lat>-94)&(lat<-82)&(toend>13))
K=np.round(c['Fmin'][cl]+O+lat/6.25e-6).astype(np.int64)
j=np.clip(np.searchsorted(ts,K),0,len(ts)-1); found=ts[j]==K
print("candidate rows in windows:",inwin.sum()," with amplitude lookup:",(inwin&found).mean()/inwin.mean())
sel=inwin&found
np.savez('out/astro_cand_amp_idx.npz',rows=np.nonzero(sel)[0],j=j[sel])
names=[f"S{L+1}c{r+1}" for r in range(4) for L in range(4)]+[f"N{L+1}c{r+1}" for r in range(4) for L in range(4)]+["P3","P4","P5","P6"]
def mean_ex(x, s):                                          # excess mean = (sum_sig - sum_ctl)/(n_sig - n_ctl)
    sg=s&(band==0); ct=s&(band>0); w=0.5                    # two control bands
    n=sg.sum()-w*ct.sum(); return (x[sg].sum()-w*x[ct].sum())/n, x[ct].mean(), n
a=amp[j]
print("\nper observable: true-partner mean vs random mean (all chips)")
for k in [32,33,34,35]+[0,4,8,12,16,20,24,28]+[1,9,17,25]:
    tm,rm,n=mean_ex(a[:,k],sel); print(f"  {names[k]:6s} true {tm:7.1f}  random {rm:7.1f}  ratio {tm/max(rm,1e-3):5.2f}")
print("\npaddle fire fraction (amp>50) for true partners, by AstroPix chip:")
for ch in range(9):
    s=sel&(chip==ch); row=[]
    for k in (32,33,34,35):
        tm,rm,n=mean_ex((a[:,k]>50).astype(float),s); row.append(f"{names[k]} {tm:5.2f}/{rm:4.2f}")
    print(f"  chip {ch}: n_true~{n:6.0f}  "+"  ".join(row))
