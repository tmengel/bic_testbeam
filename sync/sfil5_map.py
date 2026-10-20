"""Data-driven vertical map of SFIL5 South (KCU1 ASIC2, individual SiPMs) from AstroPix position.

For each channel: background-subtracted true-partner response (max - pedestal, clipped at 0)
as a function of AstroPix y, using the unbiased weights of astro_geo_common (phase-matched
candidates +1, control bands -1/2; no physics selection). Two height estimators:
  dy    = response-weighted mean y minus the plain mean y of the same events  [mm]
  slope = d ln(response) / dy from a weighted straight-line fit over 1.25 mm bins [1/mm]
For a lateral shower profile ~ exp(-(y - y_ch)^2 / 2 s^2), slope = (y_ch - y)/s^2, so both
order the channels by height. SFIL1-4 cells (KCU0, known order) are measured the same way
as the reference.
Cross-check without AstroPix: on every trigger with both KCUs, the SFIL3-4 upper-cell fraction
U = (cells 1+2) / (cells 1+2+3), S+N, tracks the shower's vertical position (and partly its
width); each array's share of the SFIL5 South energy is followed across quintiles of U.
-> out/sfil5_map.png, out/sfil5_map.npz, table on stdout."""
import numpy as np, uproot, sys, ROOT
sys.path.insert(0, '.')
from astro_geo_common import load, cluster_xy, wmean
from plotstyle import setup, finish, legend, headroom, colz
setup()

G = load(); _, _, _, X, Y = cluster_xy()
ev = uproot.open('out/events_3062.root')['events'].arrays(['has_kcu0', 'has_kcu1', 'sfil_s', 'sfil_n', 'sfil5_s'], library='np')
ki, w = G['ki'], G['w']; y = Y[G['cl']]; x = X[G['cl']]
ARR = {'A': [0, 1, 2, 3, 4, 5, 6, 7, 9, 10, 11, 12, 13, 14, 15, 16], 'B': list(range(18, 26)) + list(range(27, 35)),
       'C': list(range(36, 44)) + list(range(45, 53)), 'D': list(range(54, 62)) + list(range(63, 71))}
arr_of = {c: k for k, v in ARR.items() for c in v}
YE = np.arange(0.0, 17.51, 1.25); YC = 0.5 * (YE[1:] + YE[:-1])

def measure(a, sel):
    """a: (n, nch) amplitudes for rows sel; returns per-channel mean, dy, edy, slope, eslope, rel profile"""
    yy, ww = y[sel], w[sel]; ok = np.isfinite(yy) & (yy > 0.5) & (yy < 17.0); yy, ww, a = yy[ok], ww[ok], a[ok].clip(0, None)
    ybar = (ww * yy).sum() / ww.sum(); out = []
    for c in range(a.shape[1]):
        o = a[:, c]; S = (ww * o).sum(); m = S / ww.sum()
        if m < 3: out.append((m, np.nan, np.nan, np.nan, np.nan, np.full(len(YC), np.nan))); continue
        cy = (ww * o * yy).sum() / S; ecy = np.sqrt((ww ** 2 * o ** 2 * (yy - cy) ** 2).sum()) / S
        pm, pe, sw = wmean(yy, ww, o, YE); good = (sw > 300) & (pm > 0) & np.isfinite(pe) & (pe > 0)
        rel = pm / m
        if good.sum() >= 4:
            L, eL = np.log(pm[good]), pe[good] / pm[good]
            (s, b), cov = np.polyfit(YC[good], L, 1, w=1 / eL, cov='unscaled'); es = np.sqrt(cov[0, 0])
        else: s = es = np.nan
        out.append((m, cy - ybar, ecy, s, es, rel))
    return out

s1 = ev['has_kcu1'][ki] == 1; r5 = measure(ev['sfil5_s'][ki][s1].astype(np.float32), s1)
s0 = ev['has_kcu0'][ki] == 1
rS = measure(ev['sfil_s'][ki][s0].reshape(-1, 16), s0); rN = measure(ev['sfil_n'][ki][s0].reshape(-1, 16), s0)

print("Reference, SFIL1-4 cells (KCU0; cell 3 = beam row, cells 1+2 = SiPM1+2 on v2 boards):")
print("  layer cell |  South: mean   dy[mm]       slope[1/mm]   |  North: mean   dy[mm]       slope[1/mm]")
for L in range(4):
    for r in range(4):
        k = r * 4 + L; a, b = rS[k], rN[k]
        print("  SFIL%d  c%d  | %7.1f  %+5.2f+-%.2f  %+6.3f+-%.3f | %7.1f  %+5.2f+-%.2f  %+6.3f+-%.3f" % (L + 1, r + 1, a[0], a[1], a[2], a[3], a[4], b[0], b[1], b[2], b[3], b[4]))
print("\nSFIL5 South (KCU1 ch = 144 + index), sorted by dy within each array:")
rows = []
for k, v in ARR.items():
    live = [c for c in v if np.isfinite(r5[c][1])]
    print(f"  array {k}: {len(live)}/16 live" + ("" if live else ""))
    for c in sorted(live, key=lambda c: r5[c][1]):
        m, dy, edy, s, es, _ = r5[c]; rows.append((c, k, m, dy, edy, s, es))
        print("     ch %3d  mean %6.1f  dy %+5.2f+-%.2f  slope %+6.3f+-%.3f" % (144 + c, m, dy, edy, s, es))
R = np.array([(r[0], r[2], r[3], r[4], r[5], r[6]) for r in rows])
for k in ARR:
    sel = np.array([r[1] == k for r in rows])
    if sel.any():
        wt = 1 / R[sel, 3] ** 2
        print(f"  array {k}: weighted mean dy {np.sum(wt * R[sel, 2]) / wt.sum():+.2f} mm, spread (rms) {R[sel, 2].std():.2f} mm, n={sel.sum()}")
np.savez('out/sfil5_map.npz', ch=144 + R[:, 0].astype(int), array=np.array([r[1] for r in rows]), mean=R[:, 1], dy=R[:, 2], edy=R[:, 3],
         slope=R[:, 4], eslope=R[:, 5], ref_south=np.array([r[1:5] for r in rS]), ref_north=np.array([r[1:5] for r in rN]))

# along the fibres: response vs x should be common to all channels (attenuation only)
xs = x[s1]; okx = np.isfinite(xs) & (xs > 40) & (xs < 180); a5 = ev['sfil5_s'][ki][s1].astype(np.float32).clip(0, None)
xbar = (w[s1][okx] * xs[okx]).sum() / w[s1][okx].sum()
dx = [((w[s1][okx] * a5[okx, c] * xs[okx]).sum() / (w[s1][okx] * a5[okx, c]).sum() - xbar) for c in R[:, 0].astype(int)]
print("\nalong the fibres (x): response-weighted mean x minus mean x, live SFIL5 channels: median %+.2f mm, rms %.2f mm" % (np.median(dx), np.std(dx)))

# ---- plots
# ---- cross-check without AstroPix: array shares vs the SFIL3-4 vertical proxy
E = uproot.open('out/events_3062.root')['events'].arrays(['has_kcu0', 'has_kcu1', 'sfil_s', 'sfil_n', 'sfil5_s'], library='np')
m = (E['has_kcu0'] == 1) & (E['has_kcu1'] == 1); S_, N_ = E['sfil_s'][m].astype(float), E['sfil_n'][m].astype(float)
a5a = E['sfil5_s'][m].astype(float).clip(0, None); del E
up = (S_[:, 0, 2:] + S_[:, 1, 2:] + N_[:, 0, 2:] + N_[:, 1, 2:]).sum(1); lo = (S_[:, 2, 2:] + N_[:, 2, 2:]).sum(1)
U = up / np.maximum(up + lo, 1); tot = a5a.sum(1); px = (up + lo > 800) & (tot > 600) & ((a5a > 900).sum(1) == 0)
qb = np.digitize(U[px], np.quantile(U[px], [0.2, 0.4, 0.6, 0.8])); fr5 = a5a[px] / tot[px, None]
share = {k: np.array([fr5[qb == i][:, v].sum(1).mean() for i in range(5)]) for k, v in ARR.items()}
print("\nSFIL3-4 proxy, %d triggers: array share of SFIL5 S energy in U quintiles (low -> high), log(top/bottom):" % px.sum())
for k in ARR: print("  %s %s  %+.2f" % (k, np.round(share[k], 3), np.log(share[k][4] / share[k][0])))
# array energy share vs AstroPix y (background-subtracted)
a5t = ev['sfil5_s'][ki][s1].astype(np.float32).clip(0, None); t5 = a5t.sum(1); okE = (t5 > 200) & np.isfinite(y[s1])
ashare = {k: wmean(y[s1][okE], w[s1][okE], a5t[okE][:, v].sum(1) / t5[okE], YE) for k, v in ARR.items()}

cv = ROOT.TCanvas("cv", "", 1900, 1100); cv.Divide(3, 2, 0.004, 0.006); keep = []
cols = {'A': ROOT.kRed + 1, 'B': ROOT.kBlue + 1, 'C': ROOT.kGreen + 2, 'D': ROOT.kMagenta + 1}
cv.cd(1); ROOT.gPad.SetGrid()
fr = ROOT.TH1D("f1", "SFIL5 S: response-weighted y offset per channel;KCU1 channel;#Deltay = <y>_{resp} - <y> [mm]", 1, 143, 216)
fr.SetMinimum(-1.2); fr.SetMaximum(1.6); fr.Draw(); keep.append(fr); lg = legend(0.22, 0.74, 0.93, 0.86, 0.042, 4); keep.append(lg)
for k in ARR:
    sel = np.array([r[1] == k for r in rows])
    if not sel.any(): continue
    g = ROOT.TGraphErrors(int(sel.sum()), (144 + R[sel, 0]).astype(float), R[sel, 2].astype(float), np.zeros(sel.sum()), R[sel, 3].astype(float))
    g.SetMarkerStyle(20); g.SetMarkerColor(cols[k]); g.SetLineColor(cols[k]); g.Draw("P same"); keep.append(g); lg.AddEntry(g, f"array {k}", "p")
lg.Draw()
cv.cd(2); ROOT.gPad.SetGrid()
fr = ROOT.TH1D("f2", "SFIL1-4 cells, same estimator (S filled, N open);SFIL layer;#Deltay [mm]", 1, 0.5, 4.5); fr.SetMinimum(-1.5); fr.SetMaximum(4.6); fr.Draw(); keep.append(fr)
lg = legend(0.2, 0.68, 0.93, 0.86, 0.036, 2); keep.append(lg); ccols = (ROOT.kOrange + 7, ROOT.kAzure + 1, ROOT.kBlack, ROOT.kGray + 1)
for r_, nm in enumerate(("cell 1 (SiPM1+2)", "cell 2 (SiPM1+2)", "cell 3 (SiPM3, beam row)", "cell 4 (SiPM4)")):
    for side, rr, off, ms in (("S", rS, -0.08, 20), ("N", rN, 0.08, 24)):
        xv = np.array([L + 1 + off + (r_ - 1.5) * 0.05 for L in range(4)]); yv = np.array([rr[r_ * 4 + L][1] for L in range(4)]); ev_ = np.array([rr[r_ * 4 + L][2] for L in range(4)])
        ok = np.isfinite(yv)
        g = ROOT.TGraphErrors(int(ok.sum()), xv[ok], yv[ok], np.zeros(ok.sum()), ev_[ok]); g.SetMarkerStyle(ms); g.SetMarkerColor(ccols[r_]); g.SetLineColor(ccols[r_])
        g.Draw("P same"); keep.append(g)
        if side == "S": lg.AddEntry(g, nm, "p")
lg.Draw()
cv.cd(3); colz(); ROOT.gStyle.SetPalette(ROOT.kTemperatureMap)
order = R[np.argsort(R[:, 2]), 0].astype(int)
h = ROOT.TH2D("h3", "SFIL5 S response / mean - 1, ordered by #Deltay;AstroPix y [mm];channel rank (low #rightarrow high)", len(YC), YE[0], YE[-1], len(order), -0.5, len(order) - 0.5); keep.append(h)
for i, c in enumerate(order):
    for b, v in enumerate(r5[c][5]):
        if np.isfinite(v): h.SetBinContent(b + 1, i + 1, max(-0.6, min(0.6, v - 1)) + 1e-6)
h.SetMinimum(-0.6); h.SetMaximum(0.6); h.Draw("colz")
cv.cd(4); ROOT.gPad.SetGrid()
fr = ROOT.TH1D("f4", "#Deltay vs ln-slope estimator;#Deltay [mm];d ln(response)/dy [1/mm]", 1, -1.2, 1.2); fr.SetMinimum(-0.15); fr.SetMaximum(0.15); fr.Draw(); keep.append(fr)
for k in ARR:
    sel = np.array([r[1] == k for r in rows])
    if not sel.any(): continue
    g = ROOT.TGraphErrors(int(sel.sum()), R[sel, 2].astype(float), R[sel, 4].astype(float), R[sel, 3].astype(float), R[sel, 5].astype(float))
    g.SetMarkerStyle(20); g.SetMarkerColor(cols[k]); g.SetLineColor(cols[k]); g.Draw("P same"); keep.append(g)
cv.cd(5); ROOT.gPad.SetGrid()
fr = ROOT.TH1D("f5", "SFIL5 S array share vs AstroPix y (bkg-sub.);AstroPix y [mm];array share of SFIL5 S energy", 1, YE[0], YE[-1]); fr.SetMinimum(0); fr.SetMaximum(1.1); fr.Draw(); keep.append(fr)
lg = legend(0.22, 0.74, 0.93, 0.86, 0.042, 4); keep.append(lg)
for k in ARR:
    mm, ee, sw = ashare[k]; ok = (sw > 300) & np.isfinite(mm)
    g = ROOT.TGraphErrors(int(ok.sum()), YC[ok], mm[ok], np.zeros(ok.sum()), ee[ok]); g.SetMarkerStyle(20); g.SetMarkerColor(cols[k]); g.SetLineColor(cols[k]); g.Draw("P same"); keep.append(g); lg.AddEntry(g, f"array {k}", "p")
lg.Draw()
cv.cd(6); ROOT.gPad.SetGrid()
fr = ROOT.TH1D("f6", "No AstroPix: share vs SFIL3-4 upper fraction U;quintile of U (shower low #rightarrow high);array share of SFIL5 S energy", 1, 0.5, 5.5); fr.SetMinimum(0); fr.SetMaximum(0.8); fr.Draw(); keep.append(fr)
lg = legend(0.22, 0.74, 0.93, 0.86, 0.042, 4); keep.append(lg)
for k in ARR:
    g = ROOT.TGraph(5, np.arange(1, 6, dtype=float), share[k].astype(float)); g.SetMarkerStyle(20); g.SetLineWidth(2); g.SetMarkerColor(cols[k]); g.SetLineColor(cols[k]); g.Draw("LP same"); keep.append(g); lg.AddEntry(g, f"array {k}", "lp")
lg.Draw()
finish(cv); cv.SaveAs("out/sfil5_map.png")
