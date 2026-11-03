"""Shower shape and relative energy response, run 3062 -> out/shower.png, numbers on stdout.

There is no per-event beam energy: the August paddles were never calibrated against the PS.
Each paddle sees a narrow slice of the dispersed PS positron beam; in May the same paddles sat
at 0.48 (P4) < 0.55 (P3) < 0.57 (P6) < 0.68 GeV (P5), i.e. energy rises from the North end
(chip 0) to the South end (chip 8) of the module. Here the paddle (single-paddle triggers only)
is used as an UNCALIBRATED, ORDERED energy tag.

Energy proxies (uncalibrated), all pedestal subtracted so noise averages to zero:
  SFIL layer L: sum over cells 1, 3, 4 of sqrt(South * North) per cell (cell 2 repeats cell 1 on
                v2 boards), cell amplitude = mean(samples 6-9) - pedestal - common mode (sfil_*_w).
                The geometric mean cancels light attenuation along the fibres. Saturated channels
                (max - ped > 900 ADC) are kept and counted.
  SFIL5 South:  sum of the 64 connected KCU1 channels (sfil5_s_w, single-ended)
  Baby BCAL:    sum of the CAEN ch 0-11 pulse integrals (caen_q, South end only; North was dead)
Single-ended sums carry an attenuation bias that grows toward the South end (P5).

 1 longitudinal profile: median per layer, per paddle, normalised to P4
 2 deep / shallow ratio vs paddle (shower maximum moving deeper with energy)
 3 energy-proxy distributions per paddle (SFIL1-4 sum)
 4 relative width sigma/mu of a Gaussian fit to the peak, per paddle and proxy
 5 SFIL upper-cell share vs AstroPix y per layer (background subtracted) with erf fits
 6 fitted lateral width vs layer"""
import numpy as np, uproot, sys, ROOT
sys.path.insert(0, '.')
from astro_geo_common import load, cluster_xy, wmean
from plotstyle import setup, finish, legend, headroom, colz
setup()
ev = uproot.open('out/events_3062.root')['events'].arrays(['has_kcu0', 'has_kcu1', 'has_caen', 'pad', 'sfil_s', 'sfil_n', 'sfil_s_w', 'sfil_n_w',
                                                         'sfil5_s_w', 'caen_q'], library='np')
S, N = ev['sfil_s_w'].astype(np.float32).clip(1, None), ev['sfil_n_w'].astype(np.float32).clip(1, None)
G2 = np.sqrt(S * N)                                                  # [event, cell, layer]
EL = G2[:, [0, 2, 3], :].sum(1)                                      # per layer, cells 1, 3, 4
sat = ((ev['sfil_s'] > 900) | (ev['sfil_n'] > 900)).any(1)           # [event, layer]
conn = np.array([(c % 9) != 8 for c in range(72)])
E5 = ev['sfil5_s_w'][:, conn].sum(1); EB = ev['caen_q'][:, :12].sum(1) / 1e3   # BBCAL in 10^3 ADC x ns
PADS = (("P4", 2), ("P3", 1), ("P6", 8), ("P5", 4)); cols = (ROOT.kBlue + 1, ROOT.kGreen + 2, ROOT.kOrange + 7, ROOT.kRed + 1)
h0, h1, hc = ev['has_kcu0'] == 1, ev['has_kcu1'] == 1, ev['has_caen'] == 1
ESF = EL.sum(1); ECUT = 400.0                                       # valley between the no-shower peak (at 0) and the shower peak
shw = h0 & (ESF > ECUT)                                              # showering triggers (needs KCU0)
sel = {nm: (ev['pad'] == p) for nm, p in PADS}
print("triggers with KCU0, share below the SFIL1-4 valley (%d ADC; little or no shower in SFIL1-4):" % ECUT)
print("  " + "  ".join(f"{nm} {np.mean(ESF[sel[nm] & h0] <= ECUT):.3f}" for nm, _ in PADS))
keep = []; cv = ROOT.TCanvas("cv", "", 2100, 1250); cv.Divide(3, 2, 0.004, 0.006)

# 1 longitudinal profile
labels = ["SFIL1", "SFIL2", "SFIL3", "SFIL4", "SFIL5 S", "BBCAL S"]
med = {}
print("showering triggers: median energy proxy [ADC] per layer (SFIL1-4: sqrt(S*N), cells 1,3,4) and saturated fraction:")
for nm, _ in PADS:
    m0 = sel[nm] & shw
    med[nm] = np.r_[np.median(EL[m0], 0), np.median(E5[m0 & h1]), np.median(EB[m0 & hc])]
    print(f"  {nm}: " + "  ".join(f"{l} {v:6.0f}" for l, v in zip(labels, med[nm])) + "   sat " + " ".join(f"{x:.2f}" for x in sat[m0].mean(0)))
cv.cd(1); ROOT.gPad.SetGrid()
fr = ROOT.TH1D("f1", "Longitudinal profile relative to P4;;median / median(P4)", 6, 0.5, 6.5); fr.SetMinimum(0.4); fr.SetMaximum(2.9); keep.append(fr)
for i, l in enumerate(labels): fr.GetXaxis().SetBinLabel(i + 1, l)
fr.GetXaxis().SetLabelSize(0.05); fr.Draw(); lg = legend(0.2, 0.7, 0.93, 0.86, 0.038, 2); keep.append(lg)
for k, (nm, _) in enumerate(PADS):
    g = ROOT.TGraph(6, np.arange(1, 7, dtype=float), (med[nm] / med["P4"]).astype(float)); g.SetMarkerStyle(20 + k); g.SetMarkerColor(cols[k]); g.SetLineColor(cols[k]); g.SetLineWidth(2)
    g.Draw("LP same"); keep.append(g); lg.AddEntry(g, nm + (" (May: %s GeV)" % {"P4": "0.48", "P3": "0.55", "P6": "0.57", "P5": "0.68"}[nm]), "lp")
lg.Draw(); l = ROOT.TLine(4.5, 0.8, 4.5, 1.5); l.SetLineStyle(2); l.Draw(); keep.append(l)
t = ROOT.TLatex(4.6, 0.83, "single-ended"); t.SetTextSize(0.035); t.Draw(); keep.append(t)

# 2 deep/shallow
cv.cd(2); ROOT.gPad.SetGrid()
fr = ROOT.TH1D("f2", "Shower depth vs paddle (energy order);paddle;median ratio", 4, 0.5, 4.5); fr.SetMinimum(0.3); fr.SetMaximum(2.8); keep.append(fr)
for i, (nm, _) in enumerate(PADS): fr.GetXaxis().SetBinLabel(i + 1, nm)
fr.GetXaxis().SetLabelSize(0.06); fr.Draw(); lg = legend(0.2, 0.74, 0.93, 0.86, 0.042, 2); keep.append(lg)
for k, (name, num, den) in enumerate((("SFIL3+4 / SFIL1+2", [2, 3], [0, 1]), ("SFIL4 / SFIL2", [3], [1]))):
    v, e = [], []
    for nm, _ in PADS:
        m = sel[nm] & shw; r = EL[m][:, num].sum(1) / np.maximum(EL[m][:, den].sum(1), 1)
        v.append(float(np.median(r))); e.append(float(1.2533 * r.std() / np.sqrt(len(r))))
    print(f"  {name}: " + "  ".join(f"{nm} {x:.3f}+-{y:.3f}" for (nm, _), x, y in zip(PADS, v, e)))
    g = ROOT.TGraphErrors(4, np.arange(1, 5, dtype=float), np.array(v, float), np.zeros(4), np.array(e, float)); g.SetMarkerStyle(20 + k); g.SetLineWidth(2)
    g.SetMarkerColor(cols[3 - k]); g.SetLineColor(cols[3 - k]); g.Draw("LP same"); keep.append(g); lg.AddEntry(g, name, "lp")
lg.Draw()

# 3 + 4 energy-proxy distributions and relative widths
def peakfit(x, name, lo=0.0):
    """Gaussian fit to the shower peak: start at the histogram mode above lo, refit in [mu - 1 sigma, mu + 1.5 sigma]"""
    h = ROOT.TH1D(name, "", 120, 0, np.percentile(x, 99.5)); keep.append(h)
    for v in x: h.Fill(v)
    xu = x[x > lo]; mu0, sg0 = float(np.median(xu)), float(0.5 * (np.percentile(xu, 84) - np.percentile(xu, 16)))   # fallback: robust width
    hb = np.histogram(xu, 60)[0]; e_ = np.histogram(xu, 60)[1]; mu0 = float(0.5 * (e_[hb.argmax()] + e_[hb.argmax() + 1]))
    mu, sg, er = mu0, sg0, sg0 / np.sqrt(len(x)) / mu0
    for _ in range(3):
        f = ROOT.TF1(name + "f", "gaus", mu - 1.0 * sg, mu + 1.5 * sg); f.SetParameters(h.GetMaximum(), mu, sg)
        if int(h.Fit(f, "QRN0")) != 0 or not (0 < f.GetParameter(1) < 2 * mu0) or not (0 < abs(f.GetParameter(2)) < 2 * sg0):
            mu, sg, er = mu0, sg0, sg0 / np.sqrt(len(x)) / mu0; break
        mu, sg, er = f.GetParameter(1), abs(f.GetParameter(2)), f.GetParError(2) / f.GetParameter(1)
    return h, mu, sg, er
cv.cd(3); ROOT.gPad.SetGrid()
res = {}
for k, (nm, _) in enumerate(PADS):
    m = sel[nm] & h0 & hc; ms = m & (ESF > ECUT)
    comb = ESF[ms] / np.median(ESF[ms]) + EB[ms] / np.median(EB[ms])
    for proxy, x, lo in (("SFIL1-4", ESF[m], ECUT), ("BBCAL S", EB[ms], 0.0), ("SFIL1-4 + BBCAL", comb, 0.0)):
        h, mu, sg, er = peakfit(x, f"h{nm}{proxy}".replace(" ", "").replace("+", "p").replace(",", "").replace(".", "").replace("-", ""), lo)
        res[(nm, proxy)] = (float(sg / mu), float(er))
        if proxy == "SFIL1-4":
            h.Scale(1 / h.Integral()); h.SetLineColor(cols[k]); h.SetLineWidth(2); h.SetTitle("SFIL1-4 energy proxy per paddle;#Sigma_{layers} #sqrt{S#timesN} [ADC];fraction")
            h.SetMinimum(2e-4); h.SetMaximum(30); h.GetXaxis().SetNdivisions(505); h.Draw("hist" if k == 0 else "hist same")
            if k == 0:
                ROOT.gPad.SetLogy(); lc = ROOT.TLine(ECUT, 2e-4, ECUT, 0.3); lc.SetLineStyle(2); lc.Draw(); keep.append(lc)
lg = legend(0.45, 0.74, 0.93, 0.86, 0.042, 4); keep.append(lg)
for k, (nm, _) in enumerate(PADS): lg.AddEntry(keep[[i for i, o in enumerate(keep) if isinstance(o, ROOT.TH1D) and o.GetName() == f"h{nm}SFIL14"][0]], nm, "l")
lg.Draw()
print("\nrelative width sigma/mu of a Gaussian fit to the shower peak (showering triggers; uncalibrated; includes beam spread,")
print("leakage, saturation and the single-ended BBCAL):")
for proxy in ("SFIL1-4", "BBCAL S", "SFIL1-4 + BBCAL"):
    print(f"  {proxy:18s} " + "  ".join(f"{nm} {res[(nm, proxy)][0]:.3f}" for nm, _ in PADS))
cv.cd(4); ROOT.gPad.SetGrid()
fr = ROOT.TH1D("f4", "Shower-peak width (upper bound on resolution);paddle (energy order);#sigma / #mu", 4, 0.5, 4.5); fr.SetMinimum(0); fr.SetMaximum(0.9); keep.append(fr)
for i, (nm, _) in enumerate(PADS): fr.GetXaxis().SetBinLabel(i + 1, nm)
fr.GetXaxis().SetLabelSize(0.06); fr.Draw(); lg = legend(0.2, 0.7, 0.93, 0.86, 0.042); keep.append(lg)
for k, proxy in enumerate(("SFIL1-4", "BBCAL S", "SFIL1-4 + BBCAL")):
    g = ROOT.TGraphErrors(4, np.arange(1, 5, dtype=float), np.array([res[(nm, proxy)][0] for nm, _ in PADS]), np.zeros(4), np.array([res[(nm, proxy)][1] for nm, _ in PADS]))
    g.SetMarkerStyle(20 + k); g.SetMarkerColor(cols[k]); g.SetLineColor(cols[k]); g.SetLineWidth(2); g.Draw("LP same"); keep.append(g); lg.AddEntry(g, proxy, "lp")
lg.Draw()

# 5 + 6 lateral spread per layer from the cell boundary vs AstroPix y
G = load(); _, _, _, X, Y = cluster_xy(); ki, w = G['ki'], G['w']; y = Y[G['cl']]
okA = (ev['has_kcu0'][ki] == 1) & np.isfinite(y)
ye = np.arange(4.0, 16.01, 0.75); yc = 0.5 * (ye[1:] + ye[:-1])
cv.cd(5); ROOT.gPad.SetGrid()
fr = ROOT.TH1D("f5", "Upper-cell share, cell 1 / (cell 1 + cell 3);AstroPix y [mm];upper share", 1, 4, 16); fr.SetMinimum(0); fr.SetMaximum(1.05); fr.Draw(); keep.append(fr)
lg = legend(0.2, 0.68, 0.93, 0.86, 0.04); keep.append(lg); widths = []
for L in range(1, 4):                                       # SFIL1 carries almost no shower signal
    up, lo = G2[ki][:, 0, L], G2[ki][:, 2, L]; u = up / np.maximum(up + lo, 1); ok = okA & ((up + lo) > 50)
    m, e, sw = wmean(y[ok], w[ok], u[ok], ye); good = (sw > 300) & np.isfinite(m) & (e > 0)
    g = ROOT.TGraphErrors(int(good.sum()), yc[good].astype(float), m[good].astype(float), np.zeros(good.sum()), e[good].astype(float)); g.SetMarkerStyle(20 + L); g.SetMarkerColor(cols[L]); g.SetLineColor(cols[L])
    f = ROOT.TF1(f"erf{L}", "[0]+[1]*0.5*(1+TMath::Erf((x-[2])/(sqrt(2)*[3])))", 4, 15); f.SetParameters(0.05, 0.5, 12, 1.5); f.SetParLimits(3, 0.2, 8); f.SetLineColor(cols[L])
    g.Fit(f, "QR"); g.Draw("P same"); keep += [g, f]
    widths.append((f.GetParameter(3), f.GetParError(3), f.GetParameter(2), f.GetParError(2)))
    lg.AddEntry(g, f"SFIL{L+1}: #sigma = {f.GetParameter(3):.2f} #pm {f.GetParError(3):.2f} mm", "p")
lg.Draw()
print("\nlateral width from the cell-boundary step (erf sigma) and boundary position:")
for L, (s_, es, y0, ey0) in enumerate(widths, 1): print(f"  SFIL{L+1}: sigma {s_:.2f} +- {es:.2f} mm, boundary at AstroPix y = {y0:.2f} +- {ey0:.2f} mm")
cv.cd(6); ROOT.gPad.SetGrid()
fr = ROOT.TH1D("f6", "Lateral spread vs depth (erf width of cell boundary);SFIL layer;#sigma_{y} [mm]", 1, 1.5, 4.5); fr.SetMinimum(0); fr.SetMaximum(3.5); fr.GetXaxis().SetNdivisions(3); fr.Draw(); keep.append(fr)
g = ROOT.TGraphErrors(3, np.arange(2, 5, dtype=float), np.array([x[0] for x in widths]), np.zeros(3), np.array([x[1] for x in widths])); g.SetMarkerStyle(20); g.SetMarkerSize(1.5); g.SetLineWidth(2); g.Draw("LP same"); keep.append(g)
finish(cv); cv.SaveAs("out/shower.png")
