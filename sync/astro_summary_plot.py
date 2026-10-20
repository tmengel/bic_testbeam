"""AstroPix summary: 10 MHz phase lock, readout latency, paddle hit fraction per chip for log10LR > 2.5
matches, and the final matches by latency, chip and run time -> out/astro_summary.png"""
import numpy as np, ROOT
from plotstyle import setup, finish, legend, headroom, colz
setup()
H=np.load('out/astro_ts10_H.npy')                        # [true/shuffled][sub][256]
lat=np.load('out/astro_latency.npz'); M=np.load('out/astro_kcu0_matches.npz')
Kd=np.load('out/kcu0_amp.npz')
cv=ROOT.TCanvas("cv","",2100,1300); cv.Divide(3,2,0.002,0.002)
# 1 phase
h1=ROOT.TH1D("h1","AstroPix ts vs KCU0 trigger at 10 MHz;((K+4)/16 - ts) mod 256;10^{3} KCU0 triggers in window",256,-0.5,255.5)
h2=ROOT.TH1D("h2","",256,-0.5,255.5); h2.SetLineColor(2)
for i in range(256): h1.SetBinContent(i+1,H[0,4,i]/1e3); h2.SetBinContent(i+1,H[1,4,i]/1e3)
cv.cd(1); h1.SetLineColor(ROOT.kBlue+1); lo1=min(H[0,4].min(),H[1,4].min())*0.998/1e3; headroom(h1,lo1,H[0,4].max()/1e3,0.4); h1.Draw("hist"); h2.Draw("hist same")
lg=legend(0.21,0.72,0.75,0.85); lg.AddEntry(h1,"true pairing (19.3#sigma)","l"); lg.AddEntry(h2,"ts shuffled (3.0#sigma)","l"); lg.Draw()
# 2 latency
e=np.arange(-130,6,1.0); hm=np.histogram(lat['dm'],e)[0]-np.histogram(lat['dx'],e)[0]
h3=ROOT.TH1D("h3","Readout latency (phase-matched minus control);K - AstroPix stamp [ms];excess pairs",len(e)-1,e[0],e[-1])
for i,v in enumerate(hm): h3.SetBinContent(i+1,v)
cv.cd(2); h3.SetLineColor(ROOT.kBlue+1); h3.Draw("hist")
# 3 purity vs cut
pass
# 4 paddle vs chip for high-purity matches
sel=(M['log10LR']>2.5)&(M['pad']>=0)
h4=ROOT.TH2D("h4","Paddle hit fraction, log_{10}LR > 2.5 matches;AstroPix chip;paddle",9,-0.5,8.5,4,0.5,4.5)
for ch in range(9):
    s=sel&(M['chip']==ch)
    if s.sum()==0: continue
    for p in range(4): h4.SetBinContent(ch+1,p+1,((M['pad'][s]>>p)&1).mean())
for p,n in enumerate(("P3","P4","P5","P6")): h4.GetYaxis().SetBinLabel(p+1,n)
cv.cd(3); colz(); h4.SetMarkerSize(1.3); h4.SetMinimum(0); h4.SetMaximum(1); h4.Draw("colz text")
sel=M['log10LR']>2.5
# 5 latency of final matches by chip-independent pops
h5=ROOT.TH1D("h5","Final matches (log_{10}LR > 2.5);K - AstroPix stamp [ms];matches",130,-125,5)
for v in -M['latency_ms'][sel]: h5.Fill(v)
cv.cd(4); h5.SetLineColor(ROOT.kBlue+1); h5.Draw("hist")
# 6 per-chip matches
h6=ROOT.TH1D("h6","Final matches per chip (log_{10}LR > 2.5);AstroPix chip;10^{3} matches",9,-0.5,8.5)
for v in M['chip'][sel]: h6.Fill(v)
cv.cd(5); h6.Scale(1e-3); h6.SetLineColor(ROOT.kBlue+1); h6.Draw("hist")
# 7 matched-trigger time spread over run
h7=ROOT.TH1D("h7","Final matches vs run time;KCU0 time [s];matches / 10 s",150,0,1500)
t0=Kd['ts'][0]
for v in (M['trig_ts'][sel]-t0)*6.25e-9: h7.Fill(v)
cv.cd(6); h7.SetLineColor(ROOT.kBlue+1); h7.Draw("hist")
finish(cv); cv.SaveAs("out/astro_summary.png")
