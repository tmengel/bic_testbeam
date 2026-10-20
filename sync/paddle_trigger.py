"""Common paddle trigger: raw KCU0 paddle waveforms and the KCU0 vs CAEN paddle pattern -> out/paddle_trigger.png
 1 KCU0 raw paddle waveforms (segment 0001, ts_dump kcuwf): P3 when P3 fired (saturates at ADC 0)
   and when another paddle fired (crosstalk spike only), mean and three single triggers
 2 P(CAEN paddle | KCU0 paddle) on all triggers both recorded (trigger_table.py)
 3 paddle multiplicity per trigger, KCU0 and CAEN"""
import numpy as np, uproot, ROOT
from plotstyle import setup, finish, legend, colz
setup(); ROOT.gStyle.SetPaintTextFormat(".3f")
d = uproot.open('out/kcuwf_0001.root')['kcuwf'].arrays(['wf'], library='np'); w = d['wf'].astype(float).reshape(-1, 8, 20)
hit = (w[:, :4, 8:17] <= 5).sum(2) >= 3
keep = []; cv = ROOT.TCanvas("cv", "", 2100, 760); cv.Divide(3, 1, 0.004, 0.004)
cv.cd(1); ROOT.gPad.SetGrid()
fr = ROOT.TH1D("f1", "KCU0 paddle P3, raw ADC;sample (25 ns);ADC", 1, -0.5, 19.5); fr.SetMinimum(-5); fr.SetMaximum(290); fr.Draw(); keep.append(fr)
lg = legend(0.22, 0.72, 0.9, 0.86); keep.append(lg)
for sel, col, name in ((hit[:, 0], ROOT.kRed + 1, "P3 fired: saturates at 0"), (hit[:, 1:].any(1) & ~hit[:, 0], ROOT.kBlue + 1, "another paddle fired")):
    g = ROOT.TGraph(20, np.arange(20.0), w[sel, 0].mean(0)); g.SetLineColor(col); g.SetLineWidth(4); g.Draw("L same"); keep.append(g); lg.AddEntry(g, name + " (mean)", "l")
    for j in np.nonzero(sel)[0][:3]:
        g = ROOT.TGraph(20, np.arange(20.0), w[j, 0]); g.SetLineColor(col); g.SetLineStyle(2); g.Draw("L same"); keep.append(g)
lg.Draw()
TT = np.load('out/triggers_3062.npz'); sh = (TT['pad_kcu0'] >= 0) & (TT['pad_caen'] >= 0); p0, pc = TT['pad_kcu0'][sh], TT['pad_caen'][sh]
cv.cd(2); colz(); ROOT.gPad.SetLeftMargin(0.2)
hm = ROOT.TH2D("hm", "P(CAEN paddle | KCU0 paddle), %.0fk triggers;KCU0 paddle;CAEN paddle" % (sh.sum() / 1e3), 4, 0.5, 4.5, 4, 0.5, 4.5); keep.append(hm)
for k, nm in enumerate(("P3", "P4", "P5", "P6")):
    hm.GetXaxis().SetBinLabel(k + 1, nm); hm.GetYaxis().SetBinLabel(k + 1, nm)
    for b in range(4): hm.SetBinContent(k + 1, b + 1, ((p0 >> k) & (pc >> b) & 1).sum() / max(1, ((p0 >> k) & 1).sum()))
hm.SetMinimum(0); hm.SetMaximum(1); hm.SetMarkerSize(1.8); hm.GetXaxis().SetLabelSize(0.06); hm.GetYaxis().SetLabelSize(0.06); hm.Draw("colz text")
cv.cd(3); ROOT.gPad.SetGrid()
nb = np.array([bin(x).count('1') for x in range(16)])
h0 = ROOT.TH1D("h0", "Paddles hit per trigger;paddles hit;fraction of triggers", 5, -0.5, 4.5)
for i in range(5): h0.GetXaxis().SetBinLabel(i + 1, str(i))
h0.GetXaxis().SetLabelSize(0.065); hC = ROOT.TH1D("hC", "", 5, -0.5, 4.5); keep += [h0, hC]
for h, p in ((h0, TT['pad_kcu0'][TT['pad_kcu0'] >= 0]), (hC, TT['pad_caen'][TT['pad_caen'] >= 0])):
    for i, v in enumerate(np.bincount(nb[p], minlength=5) / len(p)): h.SetBinContent(i + 1, v)
h0.SetLineColor(ROOT.kBlack); h0.SetLineWidth(3); hC.SetLineColor(ROOT.kRed + 1); hC.SetLineWidth(3); hC.SetLineStyle(2); h0.SetMaximum(1.3); h0.Draw("hist"); hC.Draw("hist same")
lg = legend(0.45, 0.72, 0.93, 0.86); lg.AddEntry(h0, "KCU0 (dip > 65 ADC)", "l"); lg.AddEntry(hC, "CAEN (> 200 ADC)", "l"); lg.Draw(); keep.append(lg)
finish(cv); cv.SaveAs("out/paddle_trigger.png")
