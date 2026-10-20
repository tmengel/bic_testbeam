"""AstroPix position vs KCU0 (paddles, SFIL1-4), background-subtracted true partners.
Paddle panels use the common-trigger paddle pattern (KCU0, else CAEN); SFIL panels need a KCU0 waveform."""
import numpy as np, ROOT, sys
sys.path.insert(0,'.')
from astro_geo_common import load, wmean
from plotstyle import setup, finish, legend, headroom, colz
setup()
G=load(); amp=np.load('out/kcu0_amp.npz')['amp']
cl,ki,w=G['cl'],G['ki'],G['w']; chip=G['chip'][cl]; col=G['col'][cl]; row=G['row'][cl]
X=np.where(col>=0, chip*20.0+col*0.5+0.25, np.nan)          # mm along the module, chip pitch 20 mm, 0.5 mm pixels
Y=np.where(row>=0, row*0.5+0.25, np.nan)                       # mm across the module
h0=G['i_kcu0']>=0; a=np.where(h0[:,None],amp[np.maximum(G['i_kcu0'],0)],np.nan); hp=G['pad']>=0
hit=lambda p: ((G['pad']>>p)&1).astype(float)
cv=ROOT.TCanvas("cv","",1900,1150); cv.Divide(3,2,0.004,0.006)
keep=[]
def h2(name,title,x,y,ww,nx,x0,x1,ny,y0,y1):
    h=ROOT.TH2D(name,title,nx,x0,x1,ny,y0,y1); ok=np.isfinite(x)&np.isfinite(y)
    H=np.histogram2d(x[ok],y[ok],[np.linspace(x0,x1,nx+1),np.linspace(y0,y1,ny+1)],weights=ww[ok])[0]
    for i in range(nx):
        for j in range(ny): h.SetBinContent(i+1,j+1,max(H[i,j],0))
    keep.append(h); return h
# (1) all AstroPix clusters
allc=np.arange(G['n']); Xa=np.where(G['col']>=0,G['chip']*20.0+G['col']*0.5+0.25,np.nan); Ya=np.where(G['row']>=0,G['row']*0.5+0.25,np.nan)
cv.cd(1); colz(); h=h2("ha","All AstroPix clusters (row+col);x along module [mm], chip 0 = North;y [mm]",Xa,Ya,np.ones(G['n']),180,0,180,35,0,17.5); h.Draw("colz")
# (2) true KCU0 partners
cv.cd(2); colz(); h=h2("ht","Clusters with a KCU0 partner (bkg-subtracted);x [mm];y [mm]",X,Y,w,90,0,180,35,0,17.5); h.Draw("colz")
# (3) paddle fire fraction vs x
cv.cd(3); ROOT.gPad.SetGrid()
edges=np.arange(0,181,10.0); frame=ROOT.TH1D("f3","Paddle hit fraction vs AstroPix x;x [mm];P(paddle hit)",1,0,180); frame.SetMaximum(1.5); frame.SetMinimum(0); frame.Draw(); keep.append(frame)
lg=legend(0.22,0.74,0.9,0.86,0.042,4); cols=(ROOT.kRed+1,ROOT.kBlue+1,ROOT.kGreen+2,ROOT.kMagenta+1)
for p,name in enumerate(("P3","P4","P5","P6")):
    m,e,sw=wmean(X[hp],w[hp],hit(p)[hp],edges); okb=(sw>400)&np.isfinite(m)
    g=ROOT.TGraphErrors(int(okb.sum()),(0.5*(edges[1:]+edges[:-1]))[okb].astype(float),m[okb].astype(float),np.zeros(okb.sum()),e[okb].astype(float))
    g.SetMarkerStyle(20); g.SetMarkerColor(cols[p]); g.SetLineColor(cols[p]); g.Draw("P same"); keep.append(g); lg.AddEntry(g,name,"p")
lg.Draw(); keep.append(lg)
for c0 in range(0,181,20):
    l=ROOT.TLine(c0,0,c0,1.05); l.SetLineStyle(3); l.Draw(); keep.append(l)
# (4) SFIL North-South asymmetry of the beam row vs x, per layer
cv.cd(4); ROOT.gPad.SetGrid()
frame=ROOT.TH1D("f4","Beam-row (N-S)/(N+S) minus layer mean;x [mm];#Delta (N-S)/(N+S)",1,0,180); frame.SetMinimum(-0.25); frame.SetMaximum(0.42); frame.Draw(); keep.append(frame)
lg=legend(0.22,0.74,0.9,0.86,0.042,4)
for L in range(4):
    S=a[:,2*4+L]; N=a[:,16+2*4+L]; asym=(N-S)/(N+S+1.0); okk=h0&((N+S)>60)
    m0,_,_=wmean(np.zeros(okk.sum()),w[okk],asym[okk],np.array([-1,1]))
    m,e,sw=wmean(X[okk],w[okk],asym[okk]-m0[0],edges); okb=(sw>300)&np.isfinite(m)
    g=ROOT.TGraphErrors(int(okb.sum()),(0.5*(edges[1:]+edges[:-1]))[okb].astype(float),m[okb].astype(float),np.zeros(okb.sum()),e[okb].astype(float))
    g.SetMarkerStyle(20+L); g.SetMarkerColor(cols[L]); g.SetLineColor(cols[L]); g.Draw("P same"); keep.append(g); lg.AddEntry(g,f"SFIL{L+1}","p")
lg.Draw(); keep.append(lg)
# (5) SFIL2 vertical profile vs AstroPix y: energy fraction in cell 3 (beam row), cell 4, cells 1+2
cv.cd(5); ROOT.gPad.SetGrid()
frame=ROOT.TH1D("f5","SFIL2 vertical profile vs AstroPix y;y [mm];fraction of SFIL2 (S+N) energy",1,0,17.5); frame.SetMinimum(0); frame.SetMaximum(1.45); frame.Draw(); keep.append(frame)
E=[a[:,r*4+1]+a[:,16+r*4+1] for r in range(4)]; Et=E[0]+E[1]+E[2]+E[3]; okk=h0&(Et>100)
ye=np.arange(0,17.6,1.25); lg=legend(0.22,0.68,0.9,0.86,0.04)
for k,(name,frac) in enumerate((("cell 3 (SiPM3, beam row)",E[2]/np.maximum(Et,1)),("cell 4 (SiPM4)",E[3]/np.maximum(Et,1)),("cells 1+2 (SiPM1+2)",(E[0]+E[1])/np.maximum(Et,1)))):
    m,e,sw=wmean(Y[okk],w[okk],frac[okk],ye); okb=(sw>300)&np.isfinite(m)
    g=ROOT.TGraphErrors(int(okb.sum()),(0.5*(ye[1:]+ye[:-1]))[okb].astype(float),m[okb].astype(float),np.zeros(okb.sum()),e[okb].astype(float))
    g.SetMarkerStyle(20+k); g.SetMarkerColor(cols[k]); g.SetLineColor(cols[k]); g.Draw("P same"); keep.append(g); lg.AddEntry(g,name,"p")
lg.Draw(); keep.append(lg)
# (6) paddle fire fraction vs AstroPix y, in the two best-populated paddle regions
cv.cd(6); ROOT.gPad.SetGrid()
frame=ROOT.TH1D("f6","Paddle response vs AstroPix y;y [mm];P(paddle hit)",1,0,17.5); frame.SetMinimum(0); frame.SetMaximum(1.05); frame.Draw(); keep.append(frame)
lg=legend(0.3,0.2,0.85,0.36,0.045); ye=np.arange(0,17.6,1.25)
for k,(name,sel,p) in enumerate((("P3 | chip 3",(chip==3)&hp,0),("P6 | chips 5-6",((chip==5)|(chip==6))&hp,3))):
    m,e,sw=wmean(Y[sel],w[sel],hit(p)[sel],ye); okb=(sw>300)&np.isfinite(m)
    g=ROOT.TGraphErrors(int(okb.sum()),(0.5*(ye[1:]+ye[:-1]))[okb].astype(float),m[okb].astype(float),np.zeros(okb.sum()),e[okb].astype(float))
    g.SetMarkerStyle(20+k); g.SetMarkerColor(cols[0 if p==0 else 3]); g.SetLineColor(cols[0 if p==0 else 3]); g.Draw("P same"); keep.append(g); lg.AddEntry(g,name,"p")
lg.Draw(); keep.append(lg)
finish(cv); cv.SaveAs("out/astro_geo_kcu0.png")
