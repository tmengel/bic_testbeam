"""One-figure summary: every BIC stream of run 3062 on the KCU0 clock (stock decoders).
Each panel: real pairing (filled) vs a null control (red)."""
import numpy as np, uproot, glob, ROOT
from plotstyle import setup, finish, legend, headroom
setup()
tk=np.load('out/kcu0_ts.npy').astype(np.int64)
def nearest(ref, x):
    i=np.clip(np.searchsorted(ref,x),1,len(ref)-1); l,r=ref[i-1],ref[i]
    return np.where(x-l<r-x,l-x,r-x)
# --- KCU1: fixed tick offset
k1=np.concatenate([uproot.open(f)['kcu'].arrays(['ts','nsamp'],library='np')['ts'][uproot.open(f)['kcu'].arrays(['nsamp'],library='np')['nsamp']==20] for f in sorted(glob.glob('out/beam_ROC_1-00003062-*.root'))]).astype(np.int64)
k1=np.unique(k1[k1>1e8]); off=-1143665946
x=k1-off; r1=nearest(tk,x)*6.25; r1c=nearest(tk,x+int(10e6/6.25))*6.25
# --- CAEN: clock model from the matched pairs, residual of every CAEN trigger
g=uproot.open('out/caen_3062.root')['caen'].arrays(['gtt'],library='np')['gtt'].astype(float)
M=np.load('out/caen_kcu0_match.npz'); gi,ki=M['icaen'],M['ikcu0']
p=np.polyfit(g[gi],tk[ki].astype(float),1); pred=np.polyval(p,g)
o=np.argsort(g[gi]); med=np.array([np.median((tk[ki]-np.polyval(p,g[gi]))[o][max(0,k-100):k+101]) for k in range(0,len(o),50)])
pred=pred+np.interp(g,g[gi][o][::50],med)
rc=nearest(tk,np.round(pred).astype(np.int64))*6.25; rcc=nearest(tk,np.round(pred).astype(np.int64)+int(10e6/6.25))*6.25
# --- AstroPix: 10 MHz on-chip ts phase
H=np.load('out/astro_ts10_H.npy')
cv=ROOT.TCanvas("cv","",2100,760); cv.Divide(3,1,0.002,0.002); keep=[]
def panel(i,title,xt,sig,ctl,nb,lo,hi,txt,yt="entries"):
    cv.cd(i); ROOT.gPad.SetLogy()
    hs=ROOT.TH1D(f"s{i}",f"{title};{xt};{yt}",nb,lo,hi); hc=ROOT.TH1D(f"c{i}","",nb,lo,hi)
    for h,v in ((hs,sig),(hc,ctl)):
        x,w=v if isinstance(v,tuple) else (v,np.ones(len(v)))   # (values, weights) or values
        h.FillN(len(x),np.asarray(x,float),np.asarray(w,float))
    hs.SetFillColorAlpha(ROOT.kAzure+1,0.6); hs.SetLineColor(ROOT.kAzure+2); hc.SetLineColor(ROOT.kRed+1); hc.SetLineWidth(3)
    headroom(hs,0.5,max(hs.GetMaximum(),hc.GetMaximum()),0.38,log=True)
    hs.GetXaxis().SetNdivisions(505)
    hs.Draw("hist"); hc.Draw("hist same"); keep.extend([hs,hc])
    t=ROOT.TLatex(); t.SetNDC(); t.SetTextFont(42); t.SetTextSize(0.042)
    for k,line in enumerate(txt): t.DrawLatex(0.21,0.81-0.052*k,line)
    keep.append(t)
    lg=legend(0.62,0.50,0.94,0.61,0.04); lg.AddEntry(hs,"real pairing","f"); lg.AddEntry(hc,"null control","l"); lg.Draw(); keep.append(lg)
f1=(np.abs(r1)<50).mean(); f1c=(np.abs(r1c)<50).mean()
panel(1,"KCU1 #rightarrow KCU0 (H2GCROC, SFIL5)","t_{KCU0} - t_{KCU1} - offset [ns]",r1[np.abs(r1)<1e3],r1c[np.abs(r1c)<1e3],80,-1000,1000,
      ["counter offset 7.148 s, fixed","#bf{%.0f%%} of KCU1 triggers within 50 ns"%(100*f1),"control (+10 ms): %.1f%%"%(100*f1c),"residual: 0 or +25 ns (1 bunch)"])
fc=(np.abs(rc)<100).mean(); fcc=(np.abs(rcc)<100).mean()
panel(2,"CAEN #rightarrow KCU0 (V1742, paddles + Baby BCAL)","t_{KCU0} - t_{CAEN,pred} [ns]",rc[np.abs(rc)<1e3],rcc[np.abs(rcc)<1e3],200,-1000,1000,
      ["clock ratio 1.36534068, +5.4 ppm","#bf{%.0f%%} of CAEN triggers within 100 ns"%(100*fc),"control (+10 ms): %.1f%%"%(100*fcc),"residual MAD 6.4 ns"])
sig=(np.arange(256),H[0,4]/1e3); ctl=(np.arange(256),H[1,4]/1e3)
cv.cd(3)
panel(3,"AstroPix #rightarrow KCU0 (9-chip module)","(t_{KCU0} - t_{AstroPix}) mod 25.6 #mus [100 ns]",sig,ctl,256,-0.5,255.5,
      ["chip ts = 10 MHz, locked to KCU","#bf{19.3#sigma} phase lock, shuffled ts: 3.0#sigma","latency 19 / 88 #pm 5 ms","51k events matched, 81% pure"],"10^{3} entries")
ROOT.gPad.SetLogy(0); ROOT.gPad.SetLeftMargin(0.21); keep[-4].GetYaxis().SetTitleOffset(1.9); lo3=min(H[0,4].min(),H[1,4].min())*0.998/1e3; headroom(keep[-4],lo3,H[0,4].max()/1e3,0.5); keep[-4].SetFillColorAlpha(ROOT.kAzure+1,0.25)
finish(cv)
cv.SaveAs("out/sync_banger.png"); cv.SaveAs("out/sync_banger.pdf")
print("KCU1 %.3f/%.4f  CAEN %.3f/%.4f"%(f1,f1c,fc,fcc))
