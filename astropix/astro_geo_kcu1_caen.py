"""AstroPix position vs KCU1 (SFIL5 South, individual SiPMs) and CAEN channels.
True partners: background-subtracted as in astro_geo_common (no physics selection).
KCU1 and CAEN rows come from the common-trigger table (trigger_table.py), so triggers
without a complete KCU0 waveform contribute too."""
import numpy as np, uproot, glob, ROOT, sys
sys.path.insert(0,'.')
from astro_geo_common import load, wmean
from plotstyle import setup, finish, legend, headroom, colz
setup()
G=load(); w=G['w']
chip=G['chip'][G['cl']]; col=G['col'][G['cl']]; row=G['row'][G['cl']]
X=np.where(col>=0,chip*20.0+col*0.5+0.25,np.nan); Y=np.where(row>=0,row*0.5+0.25,np.nan)
# KCU1
T=[uproot.open(f)['kcu1adc'].arrays(['ts','a1'],library='np') for f in sorted(glob.glob('out/adc_beam_ROC_1-00003062-*.root'))]
t1=np.concatenate([t['ts'] for t in T]).astype(np.int64); a1=np.concatenate([t['a1'] for t in T]).astype(np.float32)
o=np.argsort(t1,kind='stable'); t1,a1=t1[o],a1[o]; u=np.r_[True,np.diff(t1)>0]; t1,a1=t1[u],a1[u]
A1=[uproot.open(f)['kcu'].arrays(['ts','nsamp'],library='np') for f in sorted(glob.glob('out/beam_ROC_1-00003062-*.root'))]
tt1=np.concatenate([a['ts'][a['nsamp']==20] for a in A1]).astype(np.int64); tt1=np.unique(tt1[tt1>1e8])   # rows of i_kcu1
x=tt1[np.maximum(G['i_kcu1'],0)]; j1=np.clip(np.searchsorted(t1,x),0,len(t1)-1); j1=np.where((G['i_kcu1']>=0)&(t1[j1]==x),j1,-1)
h1=j1>=0
arr={'A':[0,1,2,3,4,5,6,7,9,10,11,12,13,14,15,16],'B':list(range(18,26))+list(range(27,35)),
     'C':list(range(36,44))+list(range(45,53)),'D':list(range(54,62))+list(range(63,71))}
AE={k:a1[j1[h1]][:,v].clip(0,None).sum(1) for k,v in arr.items()}; tot=sum(AE.values())
# CAEN
hc=G['i_caen']>=0
neg=np.concatenate([uproot.open(f)['caenadc'].arrays(['neg'],library='np')['neg'] for f in sorted(glob.glob('out/adc_beam_CAEN-00003062-*.root'))]).astype(np.float32)
cneg=neg[G['i_caen'][hc]]
print("true partners with KCU1: %.0f, with CAEN: %.0f"%(w[h1].sum(),w[hc].sum()))
keep=[]; cols=(ROOT.kRed+1,ROOT.kBlue+1,ROOT.kGreen+2,ROOT.kMagenta+1)
cv=ROOT.TCanvas("cv","",1900,1150); cv.Divide(3,2,0.004,0.006)
ye=np.arange(0,17.6,1.25); xe=np.arange(0,181,10.0)
def curves(pad,title,xv,xedges,series,ymin,ymax,legpos,minw=200):
    cv.cd(pad); ROOT.gPad.SetGrid()
    fr=ROOT.TH1D(f"fr{pad}",title,1,xedges[0],xedges[-1]); fr.SetMinimum(ymin); fr.SetMaximum(ymax); fr.Draw(); keep.append(fr)
    lg=legend(*legpos, 0.04, 2 if len(series) > 2 else 1); keep.append(lg)
    for k,(name,xx,ww,oo_) in enumerate(series):
        m,e,sw=wmean(xx,ww,oo_,xedges); ok=(sw>minw)&np.isfinite(m)
        g=ROOT.TGraphErrors(int(ok.sum()),(0.5*(xedges[1:]+xedges[:-1]))[ok].astype(float),m[ok].astype(float),np.zeros(ok.sum()),e[ok].astype(float))
        g.SetMarkerStyle(20+k); g.SetMarkerColor(cols[k%4]); g.SetLineColor(cols[k%4]); g.Draw("P same"); keep.append(g); lg.AddEntry(g,name,"p")
    lg.Draw()
okE=tot>200
curves(1,"SFIL5 South array share vs AstroPix y;AstroPix y [mm];fraction of SFIL5 S energy",Y[h1][okE],ye,
       [(f"{k}: ch {v[0]+144}-{v[-1]+144}",Y[h1][okE],w[h1][okE],(AE[k]/np.maximum(tot,1))[okE]) for k,v in arr.items()],0,1.15,(0.2,0.7,0.93,0.86))
def relmap(pad,title,xv,xedges,obs,ww,nch,chlabel,minw=200,zr=0.6):
    cv.cd(pad); colz()
    h=ROOT.TH2D(f"h{pad}",title,len(xedges)-1,xedges[0],xedges[-1],nch,-0.5,nch-0.5); keep.append(h)
    for ch in range(nch):
        m0,_,_=wmean(np.zeros(len(ww)),ww,obs[:,ch],np.array([-1,1]))
        if not np.isfinite(m0[0]) or m0[0]<5: continue
        m,e,sw=wmean(xv,ww,obs[:,ch],xedges)
        for i in range(len(xedges)-1):
            if sw[i]>minw and np.isfinite(m[i]): h.SetBinContent(i+1,ch+1,max(-zr,min(zr,m[i]/m0[0]-1))+1e-6)
    h.SetMinimum(-zr); h.SetMaximum(zr); h.GetYaxis().SetTitle(chlabel); h.Draw("colz")
ROOT.gStyle.SetPalette(ROOT.kTemperatureMap)
fy=np.isfinite(Y[h1]); relmap(2,"SFIL5 S channel response / channel mean - 1;AstroPix y [mm]",Y[h1][fy],ye,a1[j1[h1]][fy].clip(0,None),w[h1][fy],72,"KCU1 channel - 144")
# (3) same vertical boundary seen by SFIL2 (KCU0) and SFIL5 (KCU1)
h0=G['i_kcu0']>=0; a0=np.load('out/kcu0_amp.npz')['amp'][np.maximum(G['i_kcu0'],0)]
E2=[a0[:,r*4+1]+a0[:,16+r*4+1] for r in range(4)]; up2=(E2[0]+E2[1]); lo2=E2[2]; ok2=h0&((up2+lo2)>100)
okd=(AE['B']+AE['D'])>100
curves(3,"Vertical sharing, SFIL2 and SFIL5 S;AstroPix y [mm];upper / (upper + lower)",Y,ye,
       [("SFIL2: cells 1+2 / (cells 1+2 + cell 3)",Y[ok2],w[ok2],(up2/np.maximum(up2+lo2,1))[ok2]),
        ("SFIL5 S: array D / (D + B)",Y[h1][okd],w[h1][okd],(AE['D']/np.maximum(AE['B']+AE['D'],1))[okd])],0,1.4,(0.2,0.72,0.93,0.86))
fx=np.isfinite(X[hc]); relmap(4,"CAEN channel response / channel mean - 1;AstroPix x [mm]",X[hc][fx],xe,cneg[fx],w[hc][fx],16,"CAEN channel",zr=0.8)
xe2=np.arange(0,181,20.0)
curves(5,"CAEN paddle channels vs AstroPix x;AstroPix x [mm];mean pulse height [ADC]",X[hc],xe2,
       [(f"ch {c} = {pn}",X[hc],w[hc],cneg[:,c]) for c,pn in ((12,'P3'),(13,'P4'),(14,'P5'),(15,'P6'))],0,4600,(0.2,0.72,0.93,0.86),minw=800)
# (6) common-trigger paddle signature: P(CAEN paddle hit | KCU0 paddle hit) on all shared triggers
TT=np.load('out/triggers_3062.npz'); sh=(TT['pad_kcu0']>=0)&(TT['pad_caen']>=0); p0=TT['pad_kcu0'][sh]; pc=TT['pad_caen'][sh]
cv.cd(6); colz()
hm=ROOT.TH2D("hm","P(CAEN | KCU0 paddle): %.0fk triggers, %.1f%% same;KCU0 paddle (dip > 65 ADC);CAEN paddle (neg > 200 ADC)"%(sh.sum()/1e3,100*np.mean(p0==pc)),4,0.5,4.5,4,0.5,4.5); keep.append(hm)
for k in range(4):
    hm.GetXaxis().SetBinLabel(k+1,("P3","P4","P5","P6")[k]); hm.GetYaxis().SetBinLabel(k+1,("P3","P4","P5","P6")[k])
    for b in range(4): hm.SetBinContent(k+1,b+1,((p0>>k)&(pc>>b)&1).sum()/max(1,((p0>>k)&1).sum()))
hm.SetMinimum(0); hm.SetMaximum(1); ROOT.gStyle.SetPaintTextFormat(".3f"); hm.SetMarkerSize(1.6); hm.Draw("colz text")
finish(cv); cv.SaveAs("out/astro_geo_kcu1_caen.png")
