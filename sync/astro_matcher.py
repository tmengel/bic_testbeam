"""Event-level AstroPix -> KCU0 matching, run 3062.

Candidates: common-trigger table rows K (trigger_table.py: KCU0, KCU1 or CAEN; on the
KCU0 clock; TRIG=kcu0 restricts to KCU0 complete waveforms as before) with
  * 10 MHz phase  ph = ((K+4)//16 - ts) mod 256  within +-1 of 189 (signal band)
    or of 60 / 120 (control bands: identical geometry, purely random)
  * K - stamp in the late (-26..-7 ms) or early (-97..-76 ms) latency window
Likelihood ratio per candidate (true partner vs random), factorised:
  paddle pattern (KCU0 or CAEN; 16 = none) | chip,  beam-row energy class | chip
  (3 = no KCU0 waveform),  phase residual,
  latency bin | read-position class (time to end of readout window)
p_true = signal - control/2 (background-subtracted), p_rand = control/2.
Trained on even clusters, evaluated on odd ones. For each cluster the best
candidate is kept if log10(LR) > cut; contamination from the control bands.
"""
import numpy as np, os
TT=np.load('out/triggers_3062.npz'); keep=TT['i_kcu0']>=0 if os.environ.get('TRIG')=='kcu0' else np.ones(len(TT['T']),bool)
tsK=TT['T'][keep]; k0map=TT['i_kcu0'][keep]; padT=TT['pad'][keep]
amp=np.load('out/kcu0_amp.npz')['amp']
c=np.load('out/astro_clusters.npz'); O=int(c['coarse'])-int(0.43e9/6.25)
A=np.load('out/astro_attr.npz')
F=c['Fmin']+O; TS=c['TS'].astype(np.int64); N=len(F); chip=A['chip'].astype(int); toend=A['to_end']
rows=[]
for lo,hi in ((-26,-7),(-97,-76)):
    i0=np.searchsorted(tsK,F+int(lo*1e6/6.25)); i1=np.searchsorted(tsK,F+int(hi*1e6/6.25))
    for s0 in range(0,N,50000):
        sl=slice(s0,min(s0+50000,N)); n=i1[sl]-i0[sl]
        idx=np.concatenate([np.arange(p,q) for p,q in zip(i0[sl],i1[sl])]); cl=np.repeat(np.arange(sl.start,sl.stop),n)
        ph=((tsK[idx]+4)//16-TS[cl])%256
        for b,cen in enumerate((189,60,120)):
            d=(ph-cen+128)%256-128; k=np.abs(d)<=1
            rows.append(np.c_[cl[k],idx[k],d[k],np.full(k.sum(),b)])
R=np.concatenate(rows); cl,ki,dph,band=R.T
lat=(tsK[ki]-F[cl])*6.25e-6
h0=k0map[ki]>=0; a=amp[np.maximum(k0map[ki],0)]
pad=np.where(padT[ki]>=0,padT[ki],16).astype(int)                                     # 16 paddle patterns + none
E=a[:,[9,13,25,29]].sum(1)                                                          # SFIL2-3 beam row (cell 3), S+N
eq=np.quantile(E[(band>0)&h0],[1/3,2/3]); ecl=np.where(h0,np.digitize(E,eq),3)        # 3 energy classes + no KCU0
tcl=np.digitize(toend[cl],[13,30])                                                    # read-position class
lbin=np.where(lat>-40,np.clip(((lat+26)//2).astype(int),0,9),10+np.clip(((lat+97)//2).astype(int),0,10))   # 21 latency bins
ch=chip[cl]
train=(cl%2==0); test=~train
def lr_table(keys, nkeys):
    s=np.bincount(keys[train&(band==0)],minlength=nkeys).astype(float)
    b=np.bincount(keys[train&(band>0)],minlength=nkeys)/2.0
    t=np.clip(s-b,0,None)+0.5; r=b+0.5
    return np.log10((t/t.sum())/(r/r.sum()))
L =lr_table(ch*17+pad,9*17)[ch*17+pad]
L+=lr_table(ch*4+ecl,9*4)[ch*4+ecl]
L+=lr_table(dph+1,3)[dph+1]
L+=lr_table(tcl*21+lbin,3*21)[tcl*21+lbin]
np.savez('out/astro_matcher_rows.npz',cl=cl,ki=np.nonzero(keep)[0][ki],band=band,L=L,lat=lat,dph=dph)
def best(sel):
    o=np.lexsort((-L[sel],cl[sel])); cs=cl[sel][o]; first=np.r_[True,np.diff(cs)>0]
    bl=np.full(N,-np.inf); bl[cs[first]]=L[sel][o][first]; bk=np.full(N,-1); bk[cs[first]]=ki[sel][o][first]; return bl,bk
bS,kS=best(band==0); bC1,_=best(band==1); bC2,_=best(band==2)
odd=np.arange(N)%2==1
print("test half: %d clusters"%odd.sum())
print(" cut log10LR | accepted (signal) | random (control) | est. true | purity | true / cluster")
for cut in (0,0.5,1.0,1.5,2.0,2.5,3.0):
    fs=(bS[odd]>cut).mean(); fc=0.5*((bC1[odd]>cut).mean()+(bC2[odd]>cut).mean())
    print("   %4.1f      |   %.4f          |   %.4f         |  %.4f   | %4.2f   | "%(cut,fs,fc,fs-fc,(fs-fc)/max(fs,1e-9)))
for sub,s in (("chips 3,5,6,8",np.isin(chip,[3,5,6,8])&odd),("n>=2",(c['n']>=2)&odd)):
    for cut in (1.0,2.0):
        fs=(bS[s]>cut).mean(); fc=0.5*((bC1[s]>cut).mean()+(bC2[s]>cut).mean())
        print(f"  {sub:14s} cut {cut}: accepted {fs:.4f}, random {fc:.4f}, purity {(fs-fc)/max(fs,1e-9):.2f}, matched clusters {int((fs-fc)*s.sum())}")
np.savez('out/astro_best.npz',bestL=bS,bestK=kS,bestL_ctrl1=bC1,bestL_ctrl2=bC2)
