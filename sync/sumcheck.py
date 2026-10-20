"""ADC sum check, KCU0 vs KCU1 on matched triggers, before and after the amplitude fix -> out/sumcheck.png

Old amplitude: max(samples) - pedestal (mean of samples 1-3). On a noise-only waveform the maximum of
20 samples sits several sigma above the pedestal, so every channel contributes +10-20 ADC with no signal.
New amplitude: mean of the peak samples (KCU0 6-9, KCU1 5-8) - pedestal - common mode, where the common
mode is the same quantity for the median unconnected channel of the ASIC (removes event-wide baseline
shifts). Noise then averages to zero.
Sums: KCU0 = 32 SFIL1-4 channels (S+N); KCU1 = 64 connected SFIL5 South channels.
 row 1 (old) / row 2 (new): matched pairs, decoy pairs (each KCU0 trigger with a later KCU1 trigger),
 col 3: KCU0 sum for triggers whose KCU1 sum is empty (the band at KCU1 = 0), and r vs pairing shift."""
import numpy as np, uproot, glob, ROOT
from plotstyle import setup, finish, legend, headroom, colz
setup()
OFF1 = 1143665946
def load(pat, tree, br):
    A = [uproot.open(f)[tree].arrays(br, library='np') for f in sorted(glob.glob(pat))]
    return {k: np.concatenate([a[k] for a in A]) for k in br}
K0 = load('out/adc_beam_ROC_0-00003062-*.root', 'kcuadc', ['ts', 'amp', 'win', 'cm'])
K1 = load('out/adc_beam_ROC_1-00003062-*.root', 'kcu1adc', ['ts', 'a1', 'w1', 'cm1'])
for D in (K0, K1):
    o = np.argsort(D['ts'], kind='stable'); u = np.r_[True, np.diff(D['ts'][o]) > 0]
    for k in D: D[k] = D[k][o][u]
t0 = K0['ts'].astype(np.int64); t1 = K1['ts'].astype(np.int64) + OFF1
side = np.r_[np.zeros(16, int), np.full(16, 2)]                      # ASIC of each SFIL channel: A0 South, A2 North
old0 = K0['amp'][:, :32].sum(1); new0 = (K0['win'][:, :32] - K0['cm'][:, side]).sum(1)
conn = np.array([(c % 9) != 8 for c in range(72)])
old1 = K1['a1'][:, conn].astype(float).sum(1); new1 = (K1['w1'][:, conn] - K1['cm1'][:, None]).sum(1)
i0 = np.full(len(t1), -1)
for d in (0, -4):
    j = np.clip(np.searchsorted(t0, t1 + d), 0, len(t0) - 1); hit = (t0[j] == t1 + d) & (i0 < 0); i0[hit] = j[hit]
m = i0 >= 0; a, b = i0[m], np.nonzero(m)[0]
print("matched KCU0-KCU1 pairs:", m.sum())
keep = []; cv = ROOT.TCanvas("cv", "", 2100, 1400); cv.Divide(3, 2, 0.002, 0.002)
def h2(pad, name, title, x, y, rng):
    cv.cd(pad); colz(); ROOT.gPad.SetLogz()
    h = ROOT.TH2D(name, title, 80, rng[0] / 1e3, rng[1] / 1e3, 80, rng[2] / 1e3, rng[3] / 1e3); keep.append(h)
    h.FillN(len(x), x.astype(float) / 1e3, y.astype(float) / 1e3, np.ones(len(x)))
    h.Draw("colz"); return h
r = lambda x, y: np.corrcoef(x, y)[0, 1]
for row, (nm, s0, s1, rng) in enumerate((("old: max - ped", old0, old1, (0, 5500, 0, 13000)),
                                          ("new: window - ped - CM", new0, new1, (-400, 5500, -1000, 13000)))):
    x, y = s0[a], s1[b]; yd = s1[np.roll(b, 7)]
    h2(1 + 3 * row if row == 0 else 4, f"hm{row}", f"Matched, {nm}: r = {r(x, y):.3f};KCU0 #Sigma amp [10^{{3}} ADC];KCU1 #Sigma amp [10^{{3}} ADC]", x, y, rng)
    h2(2 + 3 * row if row == 0 else 5, f"hd{row}", f"Decoys, {nm}: r = {r(x, yd):.3f};KCU0 #Sigma amp [10^{{3}} ADC];KCU1 #Sigma amp [10^{{3}} ADC]", x, yd, rng)
    print(f"{nm:40s} matched r {r(x, y):.3f}  decoy r {r(x, yd):.3f}  KCU0 sum for KCU1-empty (old KCU1 < 300): median {np.median(x[old1[b] < 300]):.0f}")
# col 3, row 1: KCU0 sum of the KCU1-empty triggers, old vs new
cv.cd(3)
e_old = old1[b] < 300; e_new = e_old                  # the same triggers in both definitions
ho = ROOT.TH1D("ho", "KCU0 #Sigma amp, triggers with KCU1 empty;KCU0 #Sigma amp [ADC];10^{3} triggers / 30 ADC", 120, -600, 3000); hn = ROOT.TH1D("hn", "", 120, -600, 3000); keep += [ho, hn]
ho.FillN(int(e_old.sum()), old0[a][e_old].astype(float), np.ones(int(e_old.sum()))); hn.FillN(int(e_new.sum()), new0[a][e_new].astype(float), np.ones(int(e_new.sum())))
ho.Scale(1e-3); hn.Scale(1e-3); ho.SetLineColor(ROOT.kGray + 2); hn.SetLineColor(ROOT.kRed + 1); headroom(ho, 0, max(ho.GetMaximum(), hn.GetMaximum()), 0.35); ho.Draw("hist"); hn.Draw("hist same")
ho.GetXaxis().SetNdivisions(505)
lg = legend(0.45, 0.62, 0.92, 0.85, 0.04); lg.SetHeader("median KCU0 #Sigma amp"); lg.AddEntry(ho, "max - ped: %.0f" % np.median(old0[a][e_old]), "l"); lg.AddEntry(hn, "window - ped - CM: %.0f" % np.median(new0[a][e_new]), "l"); lg.Draw(); keep.append(lg)
# col 3, row 2: r vs pairing shift
cv.cd(6)
ks = np.arange(-20, 21); rr = {}
for nm, s0, s1 in (("old", old0, old1), ("new", new0, new1)):
    rr[nm] = np.array([r(s0[a], s1[np.roll(b, k)]) for k in ks])
fr = ROOT.TH1D("fs", "Summed-amplitude correlation vs pairing shift;shift k (matched triggers);Pearson r", 1, -21, 21); headroom(fr, -0.1, max(rr['old'].max(), rr['new'].max()), 0.3); fr.Draw(); keep.append(fr)
lg = legend(0.21, 0.7, 0.7, 0.85); keep.append(lg)
for nm, col, ms in (("old", ROOT.kGray + 2, 24), ("new", ROOT.kRed + 1, 20)):
    g = ROOT.TGraph(len(ks), ks.astype(float), rr[nm].astype(float)); g.SetMarkerStyle(ms); g.SetMarkerColor(col); g.SetLineColor(col); g.Draw("LP same"); keep.append(g)
    lg.AddEntry(g, ("max - pedestal" if nm == "old" else "window - ped - CM") + ": r(0) = %.3f" % rr[nm][20], "lp")
lg.Draw()
finish(cv); cv.SaveAs("out/sumcheck.png")
