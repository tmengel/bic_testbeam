"""ADC sum check revisited: matched KCU0+KCU1 triggers with an empty KCU1 (SFIL5 South) -> out/k1empty_wf.png

Classes (event tree, matched KCU0+KCU1 triggers; KCU1 sum = sum of max - ped over 72 ch):
  0  empty/empty:      KCU1 < 300 and KCU0 S+N sum < 1500 ADC (20k random of 188k)
  1  KCU1 empty only:  KCU1 < 300 and KCU0 >= 1500 (all 20k)
  2  both showering:   KCU1 >= 300 and KCU0 >= 1500 (20k random)
Raw waveforms (wf_acc with WFRAW=1, out/k1e/*.root; stock decoders), classes 100-105 / 106-111 = single
triggers of class 0 / 1.
 1 KCU0: mean over all 32 SFIL1-4 channels (S+N), minus samples 1-3
 2 KCU0 SFIL3 beam row (cell 3, S and N), minus samples 1-3
 3 KCU1: mean of the 72-channel sum, minus samples 1-3
 4 CAEN: Baby BCAL South sum, minus samples 0-19
 5 single class-1 triggers: KCU0 mean over the 32 SFIL channels, raw minus samples 1-3
 6 single class-1 triggers: KCU1 72-channel sum, raw minus samples 1-3"""
import numpy as np, uproot, glob, ROOT
import sys; sys.path.insert(0, '/sphenix/user/tmengel/eic/bic_testbeam/common')
from plotstyle import setup, finish, legend, headroom, colz
setup()
def sums(pattern):
    S, W = {}, {}
    for f in sorted(glob.glob(pattern)):
        F = uproot.open(f)
        for k in F.keys(cycle=False):
            c = int(k[1:])
            if k[0] == 's': S[c] = S.get(c, 0) + F[k].values()
            if k[0] == 'w': W[c] = W.get(c, 0) + F[k].values()
    return S, W
K0, K1, KC = sums('out/k1e/kcu0_*.root'), sums('out/k1e/kcu1_*.root'), sums('out/k1e/caen.root')
mean = lambda A, c: A[0][c] / A[1][c][0]
bsub = lambda y, n=(1, 4): y - y[n[0]:n[1]].mean()
NAMES = {0: "empty / empty", 1: "KCU1 empty, KCU0 #geq 1500", 2: "both showering"}
COLS = {0: ROOT.kGray + 2, 1: ROOT.kRed + 1, 2: ROOT.kBlue + 1}
keep = []; cv = ROOT.TCanvas("cv", "", 2100, 1250); cv.Divide(3, 2, 0.004, 0.006)
def frame(pad, title, nx, ymin, ymax, logy=False):
    cv.cd(pad); ROOT.gPad.SetGrid(); ROOT.gPad.SetLogy(logy); fr = ROOT.TH1D(f"f{pad}", title, 1, -0.5, nx - 0.5); fr.SetMinimum(ymin); fr.SetMaximum(ymax); fr.Draw(); keep.append(fr)
def line(y, col, w=3, sty=1):
    g = ROOT.TGraph(len(y), np.arange(len(y), dtype=float), y.astype(float)); g.SetLineColor(col); g.SetLineWidth(w); g.SetLineStyle(sty); g.Draw("L same"); keep.append(g); return g
print("class weights (KCU0, KCU1, CAEN):", {c: (K0[1][c][0], K1[1][c][0], KC[1].get(c, [0])[0]) for c in (0, 1, 2)})
# 1 KCU0 all 32 SFIL channels
frame(1, "KCU0: mean over all 32 SFIL1-4 channels (S+N);sample (25 ns);ADC - samples 1-3", 20, -10, 200)
lg = legend(0.2, 0.66, 0.8, 0.86, 0.042); keep.append(lg)
for c in (0, 1, 2): lg.AddEntry(line(bsub(mean(K0, c)[:32].mean(0)), COLS[c]), NAMES[c], "l")
lg.Draw()
# 2 KCU0 beam row
frame(2, "KCU0 SFIL3 beam row (cell 3): S solid, N dashed;sample (25 ns);ADC - samples 1-3", 20, -20, 700)
for c in (0, 1, 2): line(bsub(mean(K0, c)[10]), COLS[c]); line(bsub(mean(K0, c)[26]), COLS[c], 3, 2)
# 3 KCU1
frame(3, "KCU1 SFIL5 South: 72-channel sum;sample (25 ns);ADC - samples 1-3", 20, -200, 9000)
for c in (0, 1, 2): line(bsub(mean(K1, c).sum(0)), COLS[c])
t = ROOT.TLatex(1, 7600, "raw pedestal, 72-ch sum: " + " / ".join(f"{mean(K1, c).sum(0)[1:4].mean():.0f}" for c in (0, 1, 2))); t.SetTextFont(42); t.SetTextSize(0.04); t.Draw(); keep.append(t)
# 4 CAEN BBCAL
nsC = int(np.nonzero(np.abs(KC[0][0]).sum(0) > 0)[0].max()) + 1
frame(4, "CAEN Baby BCAL South sum;sample (1 ns);ADC - samples 0-19", nsC, -1400, 150)
for c in (0, 1, 2):
    if c in KC[0]: line(bsub(mean(KC, c)[:12].sum(0)[:nsC], (0, 20)), COLS[c])
lg4 = legend(0.36, 0.2, 0.93, 0.42, 0.042); keep.append(lg4)
for c in (0, 1, 2): lg4.AddEntry(keep[-4 + c], NAMES[c], "l")
lg4.Draw()
# 5, 6 single class-1 triggers
frame(5, "Six KCU1-empty, KCU0 #geq 1500 triggers: KCU0 32-ch mean;sample (25 ns);ADC - samples 1-3", 20, -20, 250)
for k, c in enumerate(range(106, 112)):
    if c in K0[0]: line(bsub(K0[0][c][:32].mean(0)), [ROOT.kRed + 1, ROOT.kOrange + 7, ROOT.kMagenta + 1, ROOT.kGreen + 2, ROOT.kAzure + 1, ROOT.kBlack][k], 2)
frame(6, "Same six triggers: KCU1 72-channel sum;sample (25 ns);ADC - samples 1-3", 20, -400, 1500)
for k, c in enumerate(range(106, 112)):
    if c in K1[0]: line(bsub(K1[0][c].sum(0)), [ROOT.kRed + 1, ROOT.kOrange + 7, ROOT.kMagenta + 1, ROOT.kGreen + 2, ROOT.kAzure + 1, ROOT.kBlack][k], 2)
finish(cv); cv.SaveAs("out/k1empty_wf.png")
for c in (0, 1, 2):
    y = bsub(mean(K0, c)[:32].mean(0)); print(f"class {c} {NAMES[c]:28s} KCU0 32-ch mean: samples 0-4 {np.round(y[:5], 1)}  6-9 {np.round(y[6:10], 1)}  15-19 {np.round(y[15:], 1)}")
