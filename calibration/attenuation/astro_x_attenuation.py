"""Model-free attenuation and North/South gain ratios from the AstroPix x position, and the beam energy versus x from
the geometric mean.  Per layer, beam-row cell (cell 3): ln(A_N / A_S) = ln(g_N/g_S) + 2 x / lambda_eff for an exponential
attenuation with x measured from the fibre centre (x > 0 towards South; AstroPix x = 20 mm x chip + 0.5 mm x column - 80).
sqrt(A_N A_S) cancels the attenuation to first order, so its x dependence is the pair-spectrometer energy curve E(x) times a
constant.  Input: fitted beam-row amplitudes S3/N3 (template fits, KCU0) and layer sums S/N, joined to the AstroPix clusters.
    python3 astro_x_attenuation.py  -> astro_x_attenuation.txt/.png
"""
import numpy as np, uproot, sys, glob, ROOT
sys.path.insert(0, '/sphenix/user/tmengel/eic/bic_testbeam/common'); from plotstyle import setup, legend, label
setup(); ROOT.gROOT.SetBatch(True); out = open('astro_x_attenuation.txt', 'w')
def say(*a): print(*a, flush=True); print(*a, file=out)
def at(xv, xcf, pf): return np.nan if (xv < xcf.min() - 5 or xv > xcf.max() + 5) else np.interp(xv, xcf, pf)   # no extrapolation outside the populated bins
D = np.load('/sphenix/user/tmengel/eic/bic_testbeam/sim/analysis/data_fit_layers.npz'); S3, N3, S, N, sat, satN, pad = D['S3'].astype(float), D['N3'].astype(float), D['S'].astype(float), D['N'].astype(float), D['satS'], D['satN'], D['pad']
ts0 = uproot.open('/sphenix/user/tmengel/eic/bic_testbeam/data/kcu0_3062.root')['header']['ts0'].array(library='np')[D['trig']]
ev = uproot.open('/sphenix/user/tmengel/eic/bic_testbeam/data/events_3062.root')['events'].arrays(['trig_ts', 'astro_x', 'astro_purity', 'astro_chip'], library='np')
ets = ev['trig_ts']; order = np.argsort(ets); ets_s = ets[order]; pos = np.clip(np.searchsorted(ets_s, ts0), 0, len(ets_s) - 1); hit = ets_s[pos] == ts0; row = np.where(hit, order[pos], -1)
x = np.full(len(ts0), np.nan); x[hit] = ev['astro_x'][row[hit]] - 88.75;   pur = np.zeros(len(ts0)); pur[hit] = ev['astro_purity'][row[hit]]   # x: AstroPix frame (20 mm x chip + 0.5 mm x column) -> stack frame, chip 4 centre at 0 (P4 chip 0 = -80, P5 chip 8 = +80)
shower = (S.sum(axis=1) > 600)
m0 = shower & np.isfinite(x) & (pur >= 0.3)
say('fit-cache showers %d; with an AstroPix x of purity >= 0.3 (chips 3, 5, 6, 8): %d; x range %.0f to %.0f mm (stack frame, North negative)' % (shower.sum(), m0.sum(), np.nanmin(x[m0]), np.nanmax(x[m0])))
# model attenuation for reference: alpha exp(-d/l1) + (1-alpha) exp(-d/l2), fibre 510 mm
alpha, l1, l2, L = 0.181, 608.7, 4182., 510.
att = lambda d: alpha * np.exp(-d / l1) + (1 - alpha) * np.exp(-d / l2)
xs = np.linspace(-100, 100, 201); model_slope = np.polyfit(xs, np.log(att(L / 2 + xs) / att(L / 2 - xs)), 1)[0]
say('model (epic two-component attenuation) slope of ln(N/S) vs x over +-100 mm: %.5f per mm -> lambda_eff = -2/slope = %.0f mm' % (model_slope, -2 / model_slope))
edges = np.arange(-100, 101, 10); xc = 0.5 * (edges[1:] + edges[:-1])
c = ROOT.TCanvas('c', '', 2400, 1600); c.Divide(3, 2); keep = []
cols = [ROOT.kBlue + 1, ROOT.kGreen + 2, ROOT.kOrange + 7, ROOT.kRed + 1]
c.cd(1); fr = ROOT.gPad.DrawFrame(-100, -1.5, 100, 1.5, ';AstroPix x [mm];ln(A_{N} / A_{S}), beam cell, median per 10 mm'); keep.append(fr); lg = legend(0.15, 0.65, 0.6, 0.88)
say('\nattenuation per layer from the beam-row cell (both ends unsaturated, A_S and A_N > 50 ADC):')
for Lr in range(4):
    m = m0 & (S3[:, Lr] > 50) & (N3[:, Lr] > 50) & (sat[:, Lr] == 0) & (satN[:, Lr] == 0)
    if Lr == 0: m &= True
    r = np.log(N3[m, Lr] / S3[m, Lr]); xx = x[m]
    prof = np.array([np.median(r[(xx >= lo) & (xx < hi)]) if ((xx >= lo) & (xx < hi)).sum() > 30 else np.nan for lo, hi in zip(edges[:-1], edges[1:])])
    fin = np.isfinite(prof)
    if fin.sum() < 3: say('  SFIL%d: too few populated x bins (%d events)' % (Lr + 1, m.sum())); continue
    p = np.polyfit(xc[fin], prof[fin], 1); res = prof[fin] - np.polyval(p, xc[fin])
    pe = np.polyfit(xx, r, 1)
    say('  SFIL%d: n %6d; profile fit slope %.5f /mm -> lambda_eff %.0f mm; gN/gS at x = 0: %.3f (event fit: slope %.5f, ratio %.3f); profile residual rms %.3f; ratio of the slope to the model %.2f' % (Lr + 1, m.sum(), p[0], -2 / p[0] if p[0] else np.nan, np.exp(p[1]), pe[0], np.exp(pe[1]), res.std(), p[0] / model_slope))
    g = ROOT.TGraph(int(fin.sum()), xc[fin].astype(float), prof[fin]); g.SetMarkerStyle(20); g.SetMarkerColor(cols[Lr]); g.SetLineColor(cols[Lr]); g.Draw('lp same'); keep.append(g); lg.AddEntry(g, 'SFIL%d: #lambda_{eff} %.0f mm, g_{N}/g_{S} %.2f' % (Lr + 1, -2 / p[0], np.exp(p[1])), 'lp')
lg.Draw(); keep.append(lg); keep.append(label('data only: slope = -2 / #lambda_{eff}, intercept = ln(g_{N}/g_{S}); model slope gives %.0f mm' % (-2 / model_slope), x=0.12, y=0.93, size=0.03))
# ---- geometric mean vs x: energy curve
c.cd(2); fr2 = ROOT.gPad.DrawFrame(-100, 0.6, 100, 1.6, ';AstroPix x [mm];#sqrt{A_{N} A_{S}} / value at x = 0'); keep.append(fr2); lg2 = legend(0.15, 0.65, 0.6, 0.88)
say('\ngeometric mean sqrt(A_N A_S) of the beam cell versus x (median per 10 mm, normalised at x = 0): the energy curve E(x) up to attenuation residuals')
Ecurves = {}
for Lr in range(1, 4):   # SFIL2-4 (SFIL1 North column noisy)
    m = m0 & (S3[:, Lr] > 50) & (N3[:, Lr] > 50) & (sat[:, Lr] == 0) & (satN[:, Lr] == 0); gm = np.sqrt(S3[m, Lr] * N3[m, Lr]); xx = x[m]
    prof = np.array([np.median(gm[(xx >= lo) & (xx < hi)]) if ((xx >= lo) & (xx < hi)).sum() > 30 else np.nan for lo, hi in zip(edges[:-1], edges[1:])])
    fin = np.isfinite(prof); norm = np.interp(0., xc[fin], prof[fin]); prof /= norm; Ecurves[Lr] = prof
    p = np.polyfit(xc[fin], prof[fin], 1); say('  SFIL%d (unsaturated both ends, n %d): linear slope %.4f per mm -> %.0f %% change over 160 mm; values at x = -80/-30/+30/+80: %s' % (Lr + 1, m.sum(), p[0], 100 * p[0] * 160, ' '.join('%.3f' % at(v, xc[fin], prof[fin]) for v in (-80, -30, 30, 80))))
    g = ROOT.TGraph(int(fin.sum()), xc[fin].astype(float), prof[fin]); g.SetMarkerStyle(20); g.SetMarkerColor(cols[Lr]); g.SetLineColor(cols[Lr]); g.Draw('lp same'); keep.append(g); lg2.AddEntry(g, 'SFIL%d beam cell' % (Lr + 1), 'lp')
# layer sums over SFIL2-4 (fitted, saturation-recovered) geometric mean per layer, then summed
m = m0 & (N[:, 1:4] > 50).all(axis=1) & (S[:, 1:4] > 50).all(axis=1); tot = np.sqrt(S[m, 1:4] * N[m, 1:4]).sum(axis=1); xx = x[m]
prof = np.array([np.median(tot[(xx >= lo) & (xx < hi)]) if ((xx >= lo) & (xx < hi)).sum() > 30 else np.nan for lo, hi in zip(edges[:-1], edges[1:])]); fin = np.isfinite(prof); prof /= np.interp(0., xc[fin], prof[fin])
p = np.polyfit(xc[fin], prof[fin], 1); say('  sum of SFIL2-4 geometric-mean layer sums (fits, all showers, n %d): slope %.4f per mm -> %.0f %% over 160 mm; at -80/-30/+30/+80: %s' % (m.sum(), p[0], 100 * p[0] * 160, ' '.join('%.3f' % at(v, xc[fin], prof[fin]) for v in (-80, -30, 30, 80))))
g = ROOT.TGraph(int(fin.sum()), xc[fin].astype(float), prof[fin]); g.SetMarkerStyle(21); g.SetMarkerColor(ROOT.kBlack); g.SetLineColor(ROOT.kBlack); g.SetLineWidth(2); g.Draw('lp same'); keep.append(g); lg2.AddEntry(g, 'SFIL2-4 layer sums (fits)', 'lp'); lg2.Draw(); keep.append(lg2)
keep.append(label('geometric mean cancels an exponential attenuation: E(x) / E(0) of the pair-spectrometer arm', x=0.12, y=0.93, size=0.03))
# ---- comparison with the paddle energy ordering from the share fits (P4 3.37, P3 3.53, P6 3.73, P5 3.91 GeV, earlier model fit)
say('\nearlier model-anchored paddle energies (share fits): P4 3.37, P3 3.53, P6 3.73, P5 3.91 GeV at x = -80/-30/+30/+80 -> relative to x = 0 (interp 3.63): 0.93 / 0.97 / 1.03 / 1.08')
# ---- x distribution per paddle (AstroPix) and the x-coverage
c.cd(3); fr3 = ROOT.gPad.DrawFrame(-100, 0, 100, 0.2, ';AstroPix x [mm];fraction of matched showers'); keep.append(fr3); lg3 = legend(0.4, 0.65, 0.88, 0.88)
PB = {4: 2, 3: 1, 6: 8, 5: 4}
for p_, colr in zip((4, 3, 6, 5), cols):
    mm = m0 & (pad == p_)   # fit-cache pad is the paddle number
    if mm.sum() < 20: say('  paddle P%d: %d matched showers (chip not in the matched set)' % (p_, mm.sum())); continue
    h = ROOT.TH1F('hx%d' % p_, '', 40, -100, 100); [h.Fill(v) for v in x[mm]]; h.Scale(1. / max(1, h.GetEntries())); h.SetLineColor(colr); h.SetLineWidth(2); h.Draw('hist same'); keep.append(h); lg3.AddEntry(h, 'P%d: median x %.0f mm' % (p_, np.median(x[mm])), 'l')
    say('  paddle P%d: matched showers %d, AstroPix x median %.0f, 16-84 %% %.0f-%.0f mm' % (p_, mm.sum(), np.median(x[mm]), np.percentile(x[mm], 16), np.percentile(x[mm], 84)))
lg3.Draw(); keep.append(lg3)
# ---- N/S ratio per layer vs x as 2D for SFIL2 (quality look)
c.cd(4); Lr = 1; m = m0 & (S3[:, Lr] > 50) & (N3[:, Lr] > 50) & (sat[:, Lr] == 0) & (satN[:, Lr] == 0)
h2 = ROOT.TH2F('h2', 'SFIL2 beam cell;AstroPix x [mm];ln(A_{N}/A_{S})', 40, -100, 100, 60, -1.5, 1.5); [h2.Fill(a, b) for a, b in zip(x[m], np.log(N3[m, Lr] / S3[m, Lr]))]; h2.Draw('col'); keep.append(h2)
c.cd(5); m = m0 & (S3[:, 2] > 50) & (N3[:, 2] > 50) & (sat[:, 2] == 0) & (satN[:, 2] == 0)
h3 = ROOT.TH2F('h3', 'SFIL3 beam cell;AstroPix x [mm];ln(A_{N}/A_{S})', 40, -100, 100, 60, -1.5, 1.5); [h3.Fill(a, b) for a, b in zip(x[m], np.log(N3[m, 2] / S3[m, 2]))]; h3.Draw('col'); keep.append(h3)
c.cd(6); m = m0 & (S3[:, 1] > 50) & (N3[:, 1] > 50) & (sat[:, 1] == 0) & (satN[:, 1] == 0)
h4 = ROOT.TH2F('h4', 'SFIL2 beam cell;AstroPix x [mm];#sqrt{A_{N} A_{S}} [ADC]', 40, -100, 100, 60, 0, 1500); [h4.Fill(a, b) for a, b in zip(x[m], np.sqrt(S3[m, 1] * N3[m, 1]))]; h4.Draw('col'); keep.append(h4)
c.SaveAs('astro_x_attenuation.png'); c.SaveAs('astro_x_attenuation.pdf'); say('wrote astro_x_attenuation.png')
