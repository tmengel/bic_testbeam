"""Final AstroPix -> KCU0 event matching, run 3062 (stock decoders only).
Cross-trained likelihood ratio (train on one half of the clusters, apply to the
other), best candidate per cluster, purity vs log10(LR) calibrated on the control
phase bands. Candidates are common-trigger table rows (trigger_table.py); the KCU0
features fall into a 'no KCU0 waveform' class when only KCU1/CAEN saw the trigger.
Output: out/astro_kcu0_matches.npz"""
import numpy as np
TT=np.load('out/triggers_3062.npz'); tsK=TT['T']; amp=np.load('out/kcu0_amp.npz')['amp']
R=np.load('out/astro_matcher_rows.npz'); cl,ki,band,lat,dph=R['cl'],R['ki'],R['band'],R['lat'],R['dph']
c=np.load('out/astro_clusters.npz'); A=np.load('out/astro_attr.npz'); N=len(c['n']); O=int(c['coarse'])-int(0.43e9/6.25)
chip=A['chip'].astype(int)[cl]; toend=A['to_end'][cl]; loc=A['loc'].astype(int)[cl]; isc=A['iscol'].astype(int)[cl]
k0=TT['i_kcu0'][ki]; h0=k0>=0; a=amp[np.maximum(k0,0)]; S=lambda L,r: a[:,r*4+L]; Nn=lambda L,r: a[:,16+r*4+L]
pad=np.where(TT['pad'][ki]>=0,TT['pad'][ki],16).astype(int)
q=lambda x,n: np.where(h0,np.digitize(x,np.quantile(x[(band>0)&h0],np.linspace(0,1,n+1)[1:-1])),n)
tcl=np.digitize(toend,[13,30])
lbin=np.where(lat>-40,np.clip(((lat+26)//2).astype(int),0,9),10+np.clip(((lat+97)//2).astype(int),0,10))
feats=[((chip*8+isc*4+np.clip(loc//16,0,3))*17+pad,9*8*17),
       (chip*4+q(S(1,2)+Nn(1,2)+S(2,2)+Nn(2,2),3),9*4),
       (chip*5+np.where(h0,np.argmax(np.stack([S(1,r)+Nn(1,r) for r in range(4)],1),1),4),9*5),
       (dph+1,3),(tcl*21+lbin,63)]
L=np.zeros(len(cl))
for half in (0,1):
    tr=(cl%2==half); ap=~tr
    for keys,n in feats:
        s=np.bincount(keys[tr&(band==0)],minlength=n).astype(float); b=np.bincount(keys[tr&(band>0)],minlength=n)/2.0
        t=np.clip(s-b,0,None)+0.5; r=b+0.5; L[ap]+=np.log10((t/t.sum())/(r/r.sum()))[keys[ap]]
def best(sel):
    o=np.lexsort((-L[sel],cl[sel])); cs=cl[sel][o]; f=np.r_[True,np.diff(cs)>0]
    bl=np.full(N,-np.inf); bl[cs[f]]=L[sel][o][f]; bk=np.full(N,-1); bk[cs[f]]=ki[sel][o][f]; return bl,bk
bS,kS=best(band==0); bC1,_=best(band==1); bC2,_=best(band==2)
# purity as a function of the cut (integral) and per LR bin (differential)
cuts=np.arange(0,5.01,0.1); acc=np.array([(bS>x).sum() for x in cuts]); rnd=np.array([0.5*((bC1>x).sum()+(bC2>x).sum()) for x in cuts])
print(" log10LR >   matches   est. random   purity   est. true")
for x,s_,r_ in zip(cuts,acc,rnd):
    if abs(x*10%5)<1e-6 and s_>0: print("   %4.1f     %7d     %7.0f     %5.2f    %7.0f"%(x,s_,r_,(s_-r_)/s_,s_-r_))
edges=np.arange(0,5.01,0.25); hs=np.histogram(bS,edges)[0]; hc=0.5*(np.histogram(bC1,edges)[0]+np.histogram(bC2,edges)[0])
pur_bin=np.clip((hs-hc)/np.maximum(hs,1),0,1)
m=bS>1.0
pur=np.zeros(N); pur[m]=pur_bin[np.clip(np.digitize(bS[m],edges)-1,0,len(pur_bin)-1)]
idx=np.nonzero(m)[0]
np.savez('out/astro_kcu0_matches.npz', cluster=idx, astro_stamp=c['Fmin'][idx], astro_ts=c['TS'][idx], chip=A['chip'][idx].astype(int),
         n_hits=c['n'][idx], trig_ts=tsK[kS[idx]], trig_row=kS[idx], i_kcu0=TT['i_kcu0'][kS[idx]],
         i_kcu1=TT['i_kcu1'][kS[idx]], i_caen=TT['i_caen'][kS[idx]], pad=TT['pad'][kS[idx]], log10LR=bS[idx], purity=pur[idx],
         latency_ms=(c['Fmin'][idx]+O-tsK[kS[idx]])*6.25e-6, offset_O=O)
print("saved",len(idx),"matches with log10LR>1; expected true among them %.0f"%pur[idx].sum())
for p in (0.5,0.7,0.8,0.9):
    s=pur[idx]>=p; print("  per-match purity >= %.1f: %6d matches, expected true %6.0f"%(p,s.sum(),pur[idx][s].sum()))
for name,m in (("with KCU0 waveform",TT['i_kcu0'][kS[idx]]>=0),("KCU1/CAEN only",TT['i_kcu0'][kS[idx]]<0)):
    print(f"  {name:18s}: {m.sum():6d} matches (log10LR>1), expected true {pur[idx][m].sum():6.0f};  log10LR>2.5: {(m&(bS[idx]>2.5)).sum()}")
