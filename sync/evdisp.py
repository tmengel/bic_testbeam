"""Event displays of synchronized triggers -> out/evdisp/ev<row>.png and out/evdisp/events.json

Events: out/evdisp_events.npz (wf_keys.py: two per paddle region, most confident AstroPix match
with KCU0, KCU1 and CAEN all present). Amplitudes from out/events_3062.root, raw waveforms from
the wf_acc classes 100+i (out/wfacc/*.root). Panels:
 1 AstroPix module (9 chips, nominal 20 mm pitch, 0.5 mm pixels) with the matched cluster
   (first row/column hit; a missing row is drawn as a band) and the four paddles (lit if hit)
 2 South end cell amplitudes, drawn in AstroPix-y order (cell 4 bottom, cells 1+2 top);
   SFIL5 South arrays in their measured order (A, B, D, C bottom to top)
 3 North end cell amplitudes: SFIL1-4, same order
 4 event summary
 5 longitudinal profile: SFIL1-4 sqrt(S*N) (cells 1,3,4), SFIL5 South sum / 10, Baby BCAL South sum
 6 KCU0 waveforms: SFIL beam row (cell 3) South and North, paddles
 7 KCU1 SFIL5 South array sums
 8 CAEN waveforms: Baby BCAL South sum, paddles"""
import numpy as np, uproot, glob, json, os, ROOT
from plotstyle import setup, finish, legend, headroom, colz
setup(); ROOT.gStyle.SetPaintTextFormat(".0f")
os.makedirs('out/evdisp', exist_ok=True)
rows = np.load('out/evdisp_events.npz')['rows']
br = ['trig_ts', 't_s', 'pad', 'pad_kcu0', 'pad_caen', 'sfil_s_w', 'sfil_n_w', 'sfil5_s_w', 'caen', 'caen_q', 'astro_x', 'astro_y', 'astro_chip', 'astro_purity', 'astro_log10LR', 'astro_latency_ms']
E = uproot.open('out/events_3062.root')['events'].arrays(br, library='np', entry_start=0)
def wf(pattern, c):
    s = None
    for f in sorted(glob.glob(pattern)):
        F = uproot.open(f)
        if f"s{c}" in F: s = F[f"s{c}"].values() if s is None else s + F[f"s{c}"].values()
    return s
ARR = {'A': [0, 1, 2, 3, 4, 5, 6, 7, 9, 10, 11, 12, 13, 14, 15, 16], 'B': list(range(18, 26)) + list(range(27, 35)),
       'C': list(range(36, 44)) + list(range(45, 53)), 'D': list(range(54, 62)) + list(range(63, 71))}
PADX = {"P4": (0, 20), "P3": (40, 80), "P6": (100, 140), "P5": (160, 180)}; PBIT = {"P3": 1, "P4": 2, "P5": 4, "P6": 8}
pcol = {"P3": ROOT.kRed + 1, "P4": ROOT.kBlue + 1, "P5": ROOT.kGreen + 2, "P6": ROOT.kMagenta + 1}
out = []
for i, r in enumerate(rows):
    w0, w1, wc = wf('out/wfacc/kcu0_*.root', 100 + i), wf('out/wfacc/kcu1_*.root', 100 + i), wf('out/wfacc/caen.root', 100 + i)
    s, n, s5, ca = E['sfil_s_w'][r].astype(float), E['sfil_n_w'][r].astype(float), E['sfil5_s_w'][r].astype(float), E['caen'][r].astype(float)
    cq = E['caen_q'][r].astype(float)
    pad = int(E['pad'][r]); ax, ay, chip = float(E['astro_x'][r]), float(E['astro_y'][r]), int(E['astro_chip'][r])
    a5 = {k: float(s5[v].clip(0, None).sum()) for k, v in ARR.items()}
    lon = [float(np.sqrt(s.clip(1, None) * n.clip(1, None))[[0, 2, 3], L].sum()) for L in range(4)] + [sum(a5.values()) / 10, float(cq[:12].sum()) / 50]
    keep = []; cv = ROOT.TCanvas(f"cv{i}", "", 2400, 1150)
    body = ROOT.TPad(f"body{i}", "", 0, 0, 1, 0.955); body.Draw(); keep.append(body); body.Divide(4, 2, 0.004, 0.008)
    cv.cd = body.cd                                                                   # panels live in the body pad
    head = (f"Run 3062, trigger {r}, t = {E['t_s'][r]:.3f} s, KCU0+KCU1+CAEN, paddle "
            + ",".join(nm for nm, b in PBIT.items() if pad >= 0 and pad & b) + f";  AstroPix chip {chip}, purity {E['astro_purity'][r]:.2f}, latency {E['astro_latency_ms'][r]:.1f} ms")
    # 1 AstroPix + paddles
    cv.cd(1); ROOT.gPad.SetGrid(0, 0)
    fr = ROOT.TH2D(f"fa{i}", "AstroPix module + trigger paddles;x along module [mm] (chip 0 = North);y [mm]", 1, -5, 185, 1, -9, 22); fr.Draw(); keep.append(fr)
    for c in range(9):
        b = ROOT.TBox(c * 20, 0, c * 20 + 17.5, 17.5); b.SetFillColor(ROOT.kGray if c != chip else ROOT.kYellow - 9); b.SetLineColor(ROOT.kBlack); b.Draw("l"); keep.append(b)
        t = ROOT.TLatex(c * 20 + 6, 18.5, str(c)); t.SetTextSize(0.04); t.Draw(); keep.append(t)
    for nm, (x0, x1) in PADX.items():
        on = pad >= 0 and bool(pad & PBIT[nm]); b = ROOT.TBox(x0, -8, x1, -2); b.SetFillColor(pcol[nm] if on else ROOT.kGray + 1); b.SetFillStyle(1001 if on else 3004)
        b.Draw(); keep.append(b); t = ROOT.TLatex(x0 + 1, -7, nm); t.SetTextSize(0.045); t.SetTextColor(ROOT.kWhite if on else ROOT.kBlack); t.Draw(); keep.append(t)
    if np.isfinite(ay):
        m = ROOT.TMarker(ax, ay, 29); m.SetMarkerSize(3.5); m.SetMarkerColor(ROOT.kRed); m.Draw(); keep.append(m)
    else:
        b = ROOT.TBox(ax - 0.5, 0, ax + 0.5, 17.5); b.SetFillColor(ROOT.kRed); b.Draw(); keep.append(b)
    # 2, 3 cell maps
    for pad_i, (side, A_, with5) in ((2, ("South", s, True)), (3, ("North", n, False))):
        cv.cd(pad_i); colz(); ROOT.gPad.SetLeftMargin(0.24)
        nl = 5 if with5 else 4
        h = ROOT.TH2D(f"h{side}{i}", f"{side} end: cell amplitude, window - ped - CM [ADC];layer;cell", nl, 0.5, nl + 0.5, 4, 0.5, 4.5); keep.append(h)
        for L in range(4):
            h.GetXaxis().SetBinLabel(L + 1, f"SFIL{L+1}")
            for c in range(4): h.SetBinContent(L + 1, 4 - c, max(A_[c, L], 0.1))   # cell 1 on top, cell 4 at the bottom
        if with5:
            h.GetXaxis().SetBinLabel(5, "SFIL5 / 8")
            for k, nm in enumerate("ABDC"): h.SetBinContent(5, k + 1, a5[nm] / 8)
        for c, l in enumerate(("1 (SiPM1+2)", "2 (SiPM1+2)", "3 (SiPM3)", "4 (SiPM4)")): h.GetYaxis().SetBinLabel(4 - c, l)
        h.SetMinimum(0); h.SetMaximum(1000); h.SetMarkerSize(1.6); h.Draw("colz text")
        if with5:
            for k, nm in enumerate("ABDC"):
                t = ROOT.TLatex(5.25, k + 0.6, nm); t.SetTextSize(0.04); t.SetTextColor(ROOT.kWhite); t.Draw(); keep.append(t)
    # 4 summary
    cv.cd(4); pv = ROOT.TPaveText(0.03, 0.03, 0.97, 0.97, "NDC"); pv.SetFillColor(0); pv.SetTextAlign(12); pv.SetTextSize(0.045); pv.SetTextFont(42); keep.append(pv)
    for l in (f"Run 3062, trigger row {r}", f"t = {E['t_s'][r]:.3f} s into the run", "streams: KCU0, KCU1, CAEN",
              "paddle: " + ",".join(nm for nm, b in PBIT.items() if pad >= 0 and pad & b) + f"   (KCU0 {int(E['pad_kcu0'][r])}, CAEN {int(E['pad_caen'][r])})",
              f"AstroPix chip {chip}, x = {ax:.1f} mm, y = " + (f"{ay:.1f} mm" if np.isfinite(ay) else "n/a (no row hit)"),
              f"match: log10 LR {E['astro_log10LR'][r]:.2f}, purity {E['astro_purity'][r]:.2f}", f"readout latency {E['astro_latency_ms'][r]:.1f} ms",
              f"SFIL1-4 #Sigma#sqrt{{S#timesN}} = {sum(lon[:4]):.0f} ADC", f"SFIL5 S {sum(a5.values()):.0f} ADC, BBCAL S {cq[:12].sum()/1e3:.1f}k ADC ns"):
        pv.AddText(l)
    pv.Draw()
    # 5 longitudinal
    cv.cd(5); ROOT.gPad.SetGrid()
    h = ROOT.TH1D(f"hl{i}", "Longitudinal profile (pedestal-subtracted);;ADC", 6, 0.5, 6.5); keep.append(h)
    for k, (l, v) in enumerate(zip(["SFIL1", "SFIL2", "SFIL3", "SFIL4", "SFIL5 S/10", "BBCAL S q/50"], lon)): h.GetXaxis().SetBinLabel(k + 1, l); h.SetBinContent(k + 1, v)
    h.SetFillColor(ROOT.kAzure - 4); h.GetXaxis().SetLabelSize(0.05); h.Draw("hist bar")
    # 6 KCU0 waveforms
    cv.cd(6); ROOT.gPad.SetGrid()
    fr = ROOT.TH1D(f"f5{i}", "KCU0 waveforms: beam row S (solid) / N (dashed), paddles (thin);sample (25 ns);ADC - pedestal", 1, -0.5, 19.5); fr.SetMinimum(-100); fr.SetMaximum(1000); fr.Draw(); keep.append(fr)
    lc = (ROOT.kBlack, ROOT.kRed + 1, ROOT.kBlue + 1, ROOT.kGreen + 2)
    lg = legend(0.7, 0.58, 0.93, 0.86, 0.042); keep.append(lg)
    for L in range(4):
        for k, sty in ((8 + L, 1), (24 + L, 2)):
            g = ROOT.TGraph(20, np.arange(20.0), w0[k].astype(float)); g.SetLineColor(lc[L]); g.SetLineWidth(2); g.SetLineStyle(sty); g.Draw("L same"); keep.append(g)
            if sty == 1: lg.AddEntry(g, f"SFIL{L+1}", "l")
    for k, nm in enumerate(("P3", "P4", "P5", "P6")):
        g = ROOT.TGraph(20, np.arange(20.0), w0[32 + k].astype(float)); g.SetLineColor(pcol[nm]); g.SetLineWidth(1); g.Draw("L same"); keep.append(g)
    lg.Draw()
    # 7 KCU1 arrays
    cv.cd(7); ROOT.gPad.SetGrid()
    fr = ROOT.TH1D(f"f7{i}", "KCU1 SFIL5 South array sums;sample (25 ns);ADC - pedestal", 1, -0.5, 19.5); fr.SetMinimum(-200); fr.SetMaximum(max(1000, 1.6 * max(w1[v].sum(0).max() for v in ARR.values()))); fr.Draw(); keep.append(fr)
    lg = legend(0.7, 0.6, 0.93, 0.86, 0.042); keep.append(lg)
    for nm, col in (("C", ROOT.kGreen + 2), ("D", ROOT.kBlue + 1), ("B", ROOT.kRed + 1), ("A", ROOT.kGray + 1)):
        g = ROOT.TGraph(20, np.arange(20.0), w1[ARR[nm]].sum(0).astype(float)); g.SetLineColor(col); g.SetLineWidth(2); g.Draw("L same"); keep.append(g); lg.AddEntry(g, f"array {nm}", "l")
    lg.Draw()
    # 8 CAEN
    cv.cd(8); ROOT.gPad.SetGrid()
    ns = int(np.nonzero(np.abs(wc).sum(0) > 0)[0].max()) + 1; wc = wc[:, :ns]; xc = np.arange(ns, dtype=float)   # recorded CAEN samples
    fr = ROOT.TH1D(f"f6{i}", "CAEN: Baby BCAL South sum (black), paddles;sample (1 ns);ADC - pedestal", 1, 0, ns); fr.SetMinimum(-2200); fr.SetMaximum(600); fr.Draw(); keep.append(fr)
    g = ROOT.TGraph(ns, xc, wc[:12].sum(0).astype(float)); g.SetLineColor(ROOT.kBlack); g.SetLineWidth(2); g.Draw("L same"); keep.append(g)
    for k, nm in enumerate(("P3", "P4", "P5", "P6")):
        g = ROOT.TGraph(ns, xc, wc[12 + k].astype(float)); g.SetLineColor(pcol[nm]); g.Draw("L same"); keep.append(g)
    ROOT.TCanvas.cd(cv); tt = ROOT.TLatex(0.01, 0.975, head); tt.SetNDC(); tt.SetTextSize(0.02); tt.Draw(); keep.append(tt)
    finish(body); cv.SaveAs(f"out/evdisp/ev{r}.png")
    out.append(dict(row=int(r), t_s=float(E['t_s'][r]), pad=pad, astro=dict(x=ax if np.isfinite(ax) else None, y=ay if np.isfinite(ay) else None, chip=chip,
                purity=float(E['astro_purity'][r]), latency_ms=float(E['astro_latency_ms'][r])),
                sfil_s=s.tolist(), sfil_n=n.tolist(), sfil5_arrays=a5, sfil5=s5.tolist(), bbcal=ca[:12].tolist(), caen_pad=ca[12:].tolist(),
                wf_kcu0_beamrow_s=w0[8:12].tolist(), wf_kcu0_beamrow_n=w0[24:28].tolist(), wf_kcu0_pad=w0[32:36].tolist(),
                wf_kcu1_arrays={k: w1[v].sum(0).tolist() for k, v in ARR.items()}, wf_caen_bbcal=wc[:12].sum(0)[::2].tolist(), wf_caen_pad=wc[12:16][:, ::2].tolist()))
json.dump(out, open('out/evdisp/events.json', 'w'))
print("wrote", len(out), "event displays to out/evdisp/")
