"""Waveform-level background subtraction for the AstroPix -> trigger association -> out/wf_bkgsub.png

From wf_acc sums (out/wfacc/*.root, classes from wf_keys.py). For each AstroPix chip region the
candidate triggers in the 10 MHz signal phase band (phase 189) contain the true partner plus
random triggers; the two control bands (phases 60, 120) contain random triggers only. Per
sample:  mean(signal), mean(control), and the true-partner mean
         (S_sig - 0.5 S_ctl) / (W_sig - 0.5 W_ctl).
Columns: chip 0 (P4), chips 2-3 (P3), chips 5-6 (P6), chip 8 (P5).
Rows: KCU0 paddles (true partner), KCU0 SFIL3 beam row S and N, KCU1 SFIL5 arrays B and D (sums),
      CAEN paddle channels (true partner), CAEN Baby BCAL South sum."""
import numpy as np, uproot, glob, ROOT
from plotstyle import setup, finish, legend, headroom, colz
setup()

def sums(pattern):
    S, Q, W = {}, {}, {}
    for f in sorted(glob.glob(pattern)):
        F = uproot.open(f)
        for k in F.keys(cycle=False):
            c = int(k[1:]); v = F[k].values()
            D = {'s': S, 'q': Q, 'w': W}[k[0]]; D[c] = D.get(c, 0) + v
    return S, Q, W
K0, K1, KC = sums('out/wfacc/kcu0_*.root'), sums('out/wfacc/kcu1_*.root'), sums('out/wfacc/caen.root')

def true_mean(A, g, chans):
    S, Q, W = A; s, c = 2 * g, 2 * g + 1
    ss, sc = S[s][chans].sum(0), S[c][chans].sum(0); ws, wc = W[s][0], W[c][0]
    tm = (ss - 0.5 * sc) / (ws - 0.5 * wc)
    err = np.sqrt(Q[s][chans].sum(0) + 0.25 * Q[c][chans].sum(0)) / (ws - 0.5 * wc)
    return ss / ws, sc / wc, tm, err, ws - 0.5 * wc

GROUPS = ("chip 0 (over P4)", "chips 2-3 (over P3)", "chips 5-6 (over P6)", "chip 8 (over P5)")
ARR = {'B': list(range(18, 26)) + list(range(27, 35)), 'D': list(range(54, 62)) + list(range(63, 71))}
pcol = (ROOT.kRed + 1, ROOT.kBlue + 1, ROOT.kGreen + 2, ROOT.kMagenta + 1)      # P3 P4 P5 P6 (kcuadc order 32..35)
cv = ROOT.TCanvas("cv", "", 2200, 2000); cv.Divide(4, 5, 0.003, 0.004); keep = []
def frame(pad, title, n, ymin, ymax):
    cv.cd(pad); ROOT.gPad.SetGrid(); fr = ROOT.TH1D(f"fr{pad}", title, 1, -0.5, n - 0.5); fr.SetMinimum(ymin); fr.SetMaximum(ymax); fr.Draw(); keep.append(fr)
def line(x, y, col, width=2, style=1, e=None):
    g = ROOT.TGraphErrors(len(x), x.astype(float), y.astype(float), np.zeros(len(x)), (e if e is not None else np.zeros(len(x))).astype(float))
    g.SetLineColor(col); g.SetMarkerColor(col); g.SetLineWidth(width); g.SetLineStyle(style); g.Draw("L same"); keep.append(g); return g
xs = np.arange(20.0)
for g in range(4):
    # row 1: KCU0 paddles, true partner
    frame(1 + g, f"{GROUPS[g]}: KCU0 paddles, true partner;sample (25 ns);ADC - ped", 20, -95, 110)
    lg = legend(0.62, 0.72, 0.94, 0.85, 0.04, 2); keep.append(lg)
    for k, nm in enumerate(("P3", "P4", "P5", "P6")):
        _, _, tm, er, nt = true_mean(K0, g, [32 + k]); lg.AddEntry(line(xs, tm, pcol[k], e=er), nm, "l")
    lg.Draw(); t = ROOT.TLatex(0.5, 88, f"{nt:.0f} true partners"); t.SetTextFont(42); t.SetTextSize(0.042); t.Draw(); keep.append(t)
    # row 2: SFIL3 beam row, the subtraction itself
    frame(5 + g, "KCU0 SFIL3 beam row (cell 3), S solid / N dashed;sample (25 ns);ADC - ped", 20, -20, 1250)
    lg = legend(0.2, 0.68, 0.93, 0.86, 0.04); keep.append(lg)
    for side, ch, sty in (("S", 10, 1), ("N", 26, 2)):
        ms, mc, tm, er, _ = true_mean(K0, g, [ch])
        a = line(xs, ms, ROOT.kGray + 2, 2, sty); b = line(xs, mc, ROOT.kOrange + 7, 2, sty); c = line(xs, tm, ROOT.kBlack, 3, sty, er)
        if side == "S": lg.AddEntry(a, "signal band (true + random)", "l"); lg.AddEntry(b, "control bands (random)", "l"); lg.AddEntry(c, "true partner = sig - #frac{1}{2} ctl", "l")
    lg.Draw()
    # row 3: SFIL5 arrays B, D
    frame(9 + g, "KCU1 SFIL5 South arrays B and D, sums;sample (25 ns);ADC - ped", 20, -100, 8500)
    lg = legend(0.45, 0.66, 0.93, 0.86, 0.04, 2); keep.append(lg)
    for nm, col in (("B", ROOT.kRed + 1), ("D", ROOT.kBlue + 1)):
        ms, mc, tm, er, _ = true_mean(K1, g, ARR[nm]); c_ = line(xs, mc, col, 1, 2); t_ = line(xs, tm, col, 3, 1, er)
        lg.AddEntry(t_, f"{nm} true", "l"); lg.AddEntry(c_, f"{nm} control", "l")
    lg.Draw()
    # rows 4, 5: CAEN
    nsC = int(np.nonzero(np.abs(KC[0][0]).sum(0) > 0)[0].max()) + 1; xc = np.arange(nsC, dtype=float)   # recorded CAEN samples
    lo, hi = 0, nsC
    frame(13 + g, "CAEN paddle channels 12-15, true partner;sample (ns);ADC - ped", nsC, -1600, 200)
    lg = legend(0.62, 0.18, 0.93, 0.45, 0.042); keep.append(lg)
    for k, nm in enumerate(("P3", "P4", "P5", "P6")):
        _, _, tm, er, _ = true_mean(KC, g, [12 + k]); tm = tm[:nsC]; lg.AddEntry(line(xc, tm, pcol[k], 2), nm, "l")
    lg.Draw()
    frame(17 + g, "CAEN Baby BCAL South, sum of ch 0-11;sample (ns);ADC - ped", nsC, -1500, 200)
    ms, mc, tm, er, _ = true_mean(KC, g, list(range(12))); a_ = line(xc, ms[:nsC], ROOT.kGray + 2); b_ = line(xc, mc[:nsC], ROOT.kOrange + 7); c_ = line(xc, tm[:nsC], ROOT.kBlack, 3)
    lg = legend(0.5, 0.18, 0.93, 0.4, 0.04); keep.append(lg); lg.AddEntry(a_, "signal band", "l"); lg.AddEntry(b_, "control bands", "l"); lg.AddEntry(c_, "true partner", "l"); lg.Draw()
finish(cv); cv.SaveAs("out/wf_bkgsub.png")
for g in range(4):
    _, _, _, _, n0 = true_mean(K0, g, [0]); _, _, _, _, n1 = true_mean(K1, g, [0]); _, _, _, _, nc = true_mean(KC, g, [0])
    print(f"{GROUPS[g]:22s} true partners: KCU0 {n0:7.0f}  KCU1 {n1:7.0f}  CAEN {nc:7.0f}   signal/control weight KCU0 {K0[2][2*g][0]:.0f}/{K0[2][2*g+1][0]:.0f}")
