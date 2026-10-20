"""KCU1 -> KCU0 alignment, run 3062 -> out/kcu1_align.png (stock-decoder timestamps).

Both KCUs count the same 160 MHz clock (6.25 ns ticks) but start at different resets.
 1 coarse: each board's reset time on the host clock, L = RTCLOCK - 6.25 ns x ts (median over the run);
   L(KCU1) - L(KCU0) gives the counter offset to about 10 ms
 2 fine: for every KCU1 trigger (100k sampled), the difference ts(KCU0) - ts(KCU1) to every KCU0 trigger
   within +-20 ms of the coarse offset; random pairs are flat, the shared triggers pile up at one value
 3 exact: the most common difference is the offset; residual after it, in ticks
Panels: (1) differences around the coarse offset, 10 us bins; (2) residual of the matched pairs after
        the offset, 1-tick bins."""
import numpy as np, uproot, glob, ROOT
from plotstyle import setup, finish, legend, headroom
setup()
def kcu(pat):
    A = [uproot.open(f)['kcu'].arrays(['ts', 'nsamp', 'rt'], library='np') for f in sorted(glob.glob(pat))]
    ts = np.concatenate([a['ts'][a['nsamp'] == 20] for a in A]).astype(np.int64)
    rt = np.concatenate([a['rt'][a['nsamp'] == 20] for a in A]).astype(np.float64)
    ok = ts > 1e8
    return ts[ok], rt[ok]
t0, r0 = kcu('out/beam_ROC_0-00003062-*.root'); t1, r1 = kcu('out/beam_ROC_1-00003062-*.root')
L0, L1 = np.median(r0 - 6.25 * t0), np.median(r1 - 6.25 * t1)
coarse = (L1 - L0) / 6.25                                   # ts(KCU0) - ts(KCU1) in ticks
t0 = np.unique(t0); t1 = np.unique(t1)
rng = np.random.default_rng(5); s1 = np.sort(rng.choice(t1, 100000, replace=False))
W = int(20e6 / 6.25)                                       # 20 ms in ticks
lo = np.searchsorted(t0, s1 + int(coarse) - W); hi = np.searchsorted(t0, s1 + int(coarse) + W)
d = np.concatenate([t0[a:b] - s for a, b, s in zip(lo, hi, s1)])
vals, cnt = np.unique(d, return_counts=True); off = vals[cnt.argmax()]
print(f"coarse offset {coarse:.0f} ticks ({coarse*6.25e-9:.4f} s); exact offset {off} ticks; coarse - exact = {(coarse-off)*6.25e-6:.2f} ms")
res = d - off; m = np.abs(res) <= 8
print("residual shares (0, -4, other within +-8):", np.round([np.mean(res[m] == 0), np.mean(res[m] == -4), np.mean(~np.isin(res[m], [0, -4]))], 4))
keep = []; cv = ROOT.TCanvas("cv", "", 1500, 700); cv.Divide(2, 1, 0.002, 0.002)
cv.cd(1); ROOT.gPad.SetLogy()
x = (d - coarse) * 6.25e-6                                 # ms from the coarse offset
h1 = ROOT.TH1D("h1", "ts(KCU0) - ts(KCU1), all pairs within #pm20 ms;difference - coarse offset [ms];pairs / 10 #mus", 4000, -20, 20); keep.append(h1)
h1.FillN(len(x), x.astype(float), np.ones(len(x))); headroom(h1, 0.5, h1.GetMaximum(), 0.35, log=True); h1.SetLineColor(ROOT.kBlue + 1); h1.Draw("hist")
t = ROOT.TLatex(); t.SetNDC(); t.SetTextFont(42); t.SetTextSize(0.042); t.DrawLatex(0.21, 0.8, "spike = shared triggers"); keep.append(t)
cv.cd(2)
h3 = ROOT.TH1D("h3", "Matched pairs after the offset;ts(KCU0) - ts(KCU1) - offset [6.25 ns ticks];fraction of matched pairs", 17, -8.5, 8.5); keep.append(h3)
h3.FillN(int(m.sum()), res[m].astype(float), np.ones(int(m.sum()))); h3.Scale(1 / h3.Integral()); headroom(h3, 0, 1, 0.3)
h3.SetFillColorAlpha(ROOT.kAzure + 1, 0.6); h3.SetLineColor(ROOT.kAzure + 2); h3.Draw("hist")
t3 = ROOT.TLatex(); t3.SetNDC(); t3.SetTextFont(42); t3.SetTextSize(0.042)
t3.DrawLatex(0.21, 0.8, "0: %.1f%%    -4 (KCU1 +25 ns): %.1f%%" % (100 * np.mean(res[m] == 0), 100 * np.mean(res[m] == -4))); keep.append(t3)
finish(cv); cv.SaveAs("out/kcu1_align.png")
