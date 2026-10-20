"""CAEN -> KCU0 trigger matching, run 3062, stock-decoder timestamps.

Model: t_KCU0 = t_CAEN + coarse + off(t), off(t) found as the coincidence peak in
2 s slices every 25 s, fitted with a line, then refined with a running median of
matched-pair residuals. Shift control: same matching with CAEN moved by +10 ms.
"""
import uproot, numpy as np
import ROOT
from plotstyle import setup, finish, legend, headroom
setup()

tk = np.load('out/kcu0_ts.npy') * 6.25
c = uproot.open('out/caen_3062.root')['caen'].arrays(library='np')
tc = c['gtt'].astype(np.int64) * (1e3 / 117.1875); rtc = c['rt']
Lk = 804631.17857e9
coarse = np.percentile(rtc - tc, 1) - Lk

def peak(x, guess, win, binw):
    lo, hi = np.searchsorted(tk, [x[0] + guess - win, x[-1] + guess + win]); ref = tk[lo:hi]
    h = np.zeros(int(2 * win / binw), np.int64)
    for xi in x + guess:
        i0, i1 = np.searchsorted(ref, [xi - win, xi + win])
        h += np.bincount(((ref[i0:i1] - xi + win) / binw).astype(np.int64), minlength=len(h))[:len(h)]
    k = h.argmax(); return guess + (k + 0.5) * binw - win, h[k] / max(np.median(h), 1)

# 1) slice peaks
x0 = tc + coarse; T, O = [], []
guess = -0.322e9
for t0 in np.arange(tk[0] + 5e9, tk[-1] - 5e9, 25e9):
    sel = (x0 > t0) & (x0 < t0 + 0.5e9)
    if sel.sum() < 100: continue
    o, snr = peak(x0[sel], guess, 0.2e6, 1e3)            # +-200 us around running guess, 1 us bins
    if snr > 10: T.append(t0); O.append(o); guess = o
T, O = np.array(T), np.array(O)
p = np.polyfit(T - tk[0], O, 1)
print(f"slices with a peak: {len(T)};  drift {p[0]*1e6:+.3f} ppm,  "
      f"ratio t_KCU0/t_CAEN = {(1+p[0])*(1e3/117.1875)/6.25:.9f} KCU ticks per CAEN tick")
print("slice residual to the line [us] p05/p50/p95:", np.round(np.percentile((O - np.polyval(p, T - tk[0])) / 1e3, [5, 50, 95]), 2))

# 2) nearest-neighbour matching with the linear model, then running-median refinement
def nearest(x):
    i = np.clip(np.searchsorted(tk, x), 1, len(tk) - 1)
    l, r = tk[i - 1], tk[i]; j = np.where(x - l < r - x, i - 1, i)
    return tk[j] - x, j
pred = x0 + np.polyval(p, x0 - tk[0])
inspan = (pred > tk[0]) & (pred < tk[-1])
r, j = nearest(pred)
m = inspan & (np.abs(r) < 50e3)
# running median over 201 matched pairs
idx = np.nonzero(m)[0]; rm = np.array([np.median(r[idx[max(0,k-100):k+101]]) for k in range(0, len(idx), 50)])
corr = np.interp(pred, pred[idx[::50]], rm)
r2, j2 = nearest(pred + corr)
for name, sh in (("matched", 0.0), ("shifted +10 ms", 10e6)):
    rr, _ = nearest(pred + corr + sh)
    print(f"  {name:15s}: of {inspan.sum()} CAEN in KCU0 span, |dt|<1 us: {100*(np.abs(rr[inspan])<1e3).mean():5.1f}%  "
          f"<100 ns: {100*(np.abs(rr[inspan])<100).mean():5.1f}%  <50 us: {100*(np.abs(rr[inspan])<50e3).mean():5.1f}%")
good = inspan & (np.abs(r2) < 1e3)
print("matched residual [ns] p05/p50/p95:", np.round(np.percentile(r2[good], [5, 50, 95]), 1), " MAD %.1f ns" % np.median(np.abs(r2[good] - np.median(r2[good]))))
np.savez('out/caen_kcu0_match.npz', icaen=np.nonzero(good)[0], ikcu0=j2[good], res=r2[good],
         tpred=pred + corr, inspan=inspan)                  # every CAEN trigger on the KCU0 clock [ns]

cv = ROOT.TCanvas("cv", "", 2100, 700); cv.Divide(3, 1, 0.002, 0.002)
g = ROOT.TGraph(len(T), ((T - tk[0]) / 1e9).astype(float), ((O - np.polyval(p, T - tk[0])) / 1e3).astype(float))
g.SetTitle(f"Slice offset - linear fit, drift {p[0]*1e6:+.2f} ppm;KCU0 time [s];offset - line [#mus]"); g.SetMarkerStyle(20); g.SetMarkerSize(1.2)
cv.cd(1); g.Draw("AP")
h1 = ROOT.TH1D("h1", "Nearest-KCU0 residual, all CAEN triggers;KCU0 - predicted [#mus];CAEN triggers", 200, -5, 5)
h2 = ROOT.TH1D("h2", "", 200, -5, 5); h2.SetLineColor(2)
for v in r2[inspan] / 1e3: h1.Fill(v)
rr, _ = nearest(pred + corr + 10e6)
for v in rr[inspan] / 1e3: h2.Fill(v)
cv.cd(2); ROOT.gPad.SetLogy(); headroom(h1, 0.5, max(h1.GetMaximum(), h2.GetMaximum()), 0.33, log=True); h1.SetLineColor(ROOT.kBlue + 1); h1.Draw("hist"); h2.Draw("hist same")
leg = legend(0.21, 0.72, 0.9, 0.85); leg.AddEntry(h1, "nearest KCU0 trigger", "l"); leg.AddEntry(h2, "CAEN shifted +10 ms (control)", "l"); leg.Draw()
h3 = ROOT.TH1D("h3", "Zoom on matched pairs;KCU0 - predicted [ns];10^{3} CAEN triggers", 100, -50, 50)
for v in r2[good]: h3.Fill(v)
cv.cd(3); h3.GetXaxis().SetNdivisions(505); h3.Scale(1e-3); h3.SetLineColor(ROOT.kBlue + 1); h3.Draw("hist")
finish(cv); cv.SaveAs("out/caen_kcu0.png")
