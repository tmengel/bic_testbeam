"""Response versus time in run 3064 from the hadd-ed gain_stability.C output (12.3 h, 600 s bins).

    python3 gain_stability_figs.py [kcu0_file] [kcu1_file]  -> gain_stability_{rate,kcu0,ns,ped,kcu1}.png, gain_stability.txt

Per 600 s bin the median of the fitted-amplitude distribution (column of the TH2) is used, divided by the median over the whole
run of the same channel / sum; the 10-min medians have 0.2-0.5 % statistical error at 3064's rates. The TAC period at the end
of the run (low current) is where the single-paddle rate falls below 20 % of its run median.
"""
import sys, os, numpy as np, ROOT
sys.path.insert(0, '/sphenix/user/tmengel/eic/bic_testbeam/common'); from plotstyle import setup, finish, legend, headroom
setup()
T = '/sphenix/tg/tg01/jets/tmengel/eic'
F0 = sys.argv[1] if len(sys.argv) > 1 else T + '/kcu0/stability/stability_3064.root'
F1 = sys.argv[2] if len(sys.argv) > 2 else T + '/kcu1/stability/stability_3064.root'
f0, f1 = ROOT.TFile.Open(F0), ROOT.TFile.Open(F1)
out = open('gain_stability.txt', 'w')
def say(*a): s = ' '.join(str(x) for x in a); print(s); out.write(s + '\n')
keep = []

def h2arr(h):
    nx, ny = h.GetNbinsX(), h.GetNbinsY(); a = np.zeros((nx, ny)); yc = np.array([h.GetYaxis().GetBinCenter(j + 1) for j in range(ny)])
    for i in range(nx):
        for j in range(ny): a[i, j] = h.GetBinContent(i + 1, j + 1)
    return a, yc, np.array([h.GetXaxis().GetBinCenter(i + 1) for i in range(nx)])
def qmed(col, yc):
    """median of a binned distribution, interpolated linearly inside the bin that holds the 50 % point"""
    n = col.sum(); c = np.cumsum(col); j = int(np.searchsorted(c, 0.5 * n)); w = yc[1] - yc[0]
    prev = c[j - 1] if j > 0 else 0.; return yc[j] - w / 2 + w * (0.5 * n - prev) / max(col[j], 1e-9)
def medians(h, nmin=200):
    """per time bin the interpolated median of the column, nan where fewer than nmin entries; plus the run median"""
    a, yc, xc = h2arr(h); med = np.full(len(xc), np.nan)
    for i in range(len(xc)):
        if a[i].sum() >= nmin: med[i] = qmed(a[i], yc)
    tot = a.sum(axis=0); run = qmed(tot, yc) if tot.sum() else np.nan
    return xc, med, run, a.sum(axis=1)
def graph(x, y, color, marker=20, size=0.8):
    m = np.isfinite(y); g = ROOT.TGraph(int(m.sum()), np.ascontiguousarray(x[m], float), np.ascontiguousarray(y[m], float)); g.SetMarkerColor(color); g.SetLineColor(color); g.SetMarkerStyle(marker); g.SetMarkerSize(size); keep.append(g); return g
def frame(name, title, xmax, ymin, ymax):
    fr = ROOT.TH2D(name, title, 10, 0, xmax, 10, ymin, ymax); fr.SetStats(0); keep.append(fr); fr.Draw(); return fr
def mkleg(*a, **k):
    lg = legend(*a, **k); keep.append(lg); return lg
COL = [ROOT.kBlack, ROOT.kRed + 1, ROOT.kBlue + 1, ROOT.kGreen + 2, ROOT.kMagenta + 1, ROOT.kOrange + 7]
H2H = 3600.

# ---------------------------------------------------------------- rates and the TAC period
hN = {p: f0.Get('hN_%d' % p) for p in (0, 3, 4, 5, 6)}
nb = hN[0].GetNbinsX(); tb = np.array([hN[0].GetXaxis().GetBinCenter(i + 1) for i in range(nb)]) / H2H
rate = {p: np.array([hN[p].GetBinContent(i + 1) for i in range(nb)]) / 600. for p in hN}
live = rate[0] > 0; medr = np.median(rate[0][live]); high = np.flatnonzero(rate[0] > 0.2 * medr); tac = tb[high[-1] + 1:]   # the TAC period starts after the last high-rate bin; beam-off gaps inside the run stay nan in the medians
say('rate per 600 s bin [Hz]:', ' '.join('%.0f' % r for r in rate[0]))
tmax_h = tb[live].max() + 0.1
say('run 3064: %d bins of 600 s with fitted showers; single-paddle shower rate median %.0f Hz; beam-off gaps (rate < 20 %% of it) inside the beam period: %d bins; TAC period (after the last high-rate bin): %d bins from %.2f h' % (live.sum(), medr, int((rate[0][:high[-1] + 1] < 0.2 * medr).sum()), len(tac), tac.min() if len(tac) else -1))
tac0 = tac.min() if len(tac) else tmax_h
cv = ROOT.TCanvas('c0', '', 900, 650); keep.append(cv); cv.SetLogy()
fr = frame('fr0', 'run 3064: triggers with a fitted beam-cell pulse per 600 s;time in run [h];rate [Hz]', tmax_h, 0.5, 5000)
lg = mkleg(0.20, 0.72, 0.60, 0.87, ncol=2)
for i, p in enumerate((0, 3, 4, 5, 6)):
    g = graph(tb, np.where(rate[p] > 0, rate[p], np.nan), COL[i]); g.Draw('P'); lg.AddEntry(g, 'any single paddle' if p == 0 else 'P%d' % p, 'p')
lg.Draw(); ln = ROOT.TLine(tac0, 0.5, tac0, 5000); ln.SetLineStyle(2); ln.Draw(); keep.append(ln)
tx = ROOT.TLatex(); tx.SetNDC(); tx.SetTextSize(0.035); tx.DrawLatex(0.62, 0.80, 'TAC period from %.2f h' % tac0); keep.append(tx)
finish(cv); cv.SaveAs('gain_stability_rate.png')

# ---------------------------------------------------------------- KCU0 beam-cell medians per layer (paddle: any) and the South total
cv = ROOT.TCanvas('c1', '', 1500, 650); cv.Divide(2, 1, 0.002, 0.002); keep.append(cv)
cv.cd(1); fr = frame('fr1', 'SFIL1-4 South beam cell, all single-paddle showers;time in run [h];median A_{fit} / run median', tmax_h, 0.90, 1.14)
lg = mkleg(0.19, 0.60, 0.62, 0.87, size=0.032); say('\nKCU0 South beam cell: 10-min medians relative to the run median (before the TAC period)')
for L in range(4):
    xc, med, run, nn = medians(f0.Get('h2_%d_0' % (8 * L + 2))); r = med / run; g = graph(xc / H2H, r, COL[L]); g.Draw('P'); lg.AddEntry(g, 'SFIL%d S cell 3 (%.0f ADC)' % (L + 1, run), 'p')
    pre = np.isfinite(r) & (xc / H2H < tac0); p = np.polyfit(xc[pre] / H2H, r[pre], 1) if pre.sum() > 2 else (np.nan, np.nan)
    say('  SFIL%d S cell 3: run median %.0f ADC, rms of the bin medians %.2f %%, slope %+.3f %%/h, TAC-period median %.3f' % (L + 1, run, 100 * np.nanstd(r[pre]), 100 * p[0], np.nanmedian(r[np.isfinite(r) & (xc / H2H >= tac0)])))
for L in range(4):
    xc, med, run, nn = medians(f0.Get('hS_%d_0' % (L + 1))); r = med / run; pre = np.isfinite(r) & (xc / H2H < tac0)
    say('  SFIL%d S layer sum: run median %.0f ADC, rms %.2f %%' % (L + 1, run, 100 * np.nanstd(r[pre])))
xc, med, run, nn = medians(f0.Get('hT_0')); r = med / run; g = graph(xc / H2H, r, ROOT.kGray + 2, 24); g.Draw('P'); lg.AddEntry(g, 'SFIL1-4 South total (%.0f ADC)' % run, 'p'); lg.Draw()
pre = np.isfinite(r) & (xc / H2H < tac0); say('  SFIL1-4 South total: run median %.0f ADC, rms %.2f %%, slope %+.3f %%/h' % (run, 100 * np.nanstd(r[pre]), 100 * np.polyfit(xc[pre] / H2H, r[pre], 1)[0]))
ln = ROOT.TLine(tac0, 0.90, tac0, 1.14); ln.SetLineStyle(2); ln.Draw(); keep.append(ln)
cv.cd(2); fr = frame('fr2', 'SFIL3 South beam cell per paddle;time in run [h];median A_{fit} / run median', tmax_h, 0.90, 1.14)
lg = mkleg(0.20, 0.66, 0.55, 0.87); say('\nSFIL3 S cell 3 per paddle:')
for i, p in enumerate((3, 4, 5, 6)):
    xc, med, run, nn = medians(f0.Get('h2_18_%d' % p)); r = med / run; g = graph(xc / H2H, r, COL[i + 1]); g.Draw('P'); lg.AddEntry(g, 'P%d (%.0f ADC)' % (p, run), 'p')
    pre = np.isfinite(r) & (xc / H2H < tac0); say('  P%d: run median %.0f ADC, rms %.2f %%, slope %+.3f %%/h' % (p, run, 100 * np.nanstd(r[pre]), 100 * np.polyfit(xc[pre] / H2H, r[pre], 1)[0]))
lg.Draw(); ln2 = ROOT.TLine(tac0, 0.90, tac0, 1.14); ln2.SetLineStyle(2); ln2.Draw(); keep.append(ln2)
finish(cv); cv.SaveAs('gain_stability_kcu0.png')

# ---------------------------------------------------------------- North / South and pedestals
cv = ROOT.TCanvas('c2', '', 1500, 650); cv.Divide(2, 1, 0.002, 0.002); keep.append(cv)
cv.cd(1); fr = frame('fr3', 'beam cell ln(A_{N}/A_{S}) - run median, all single-paddle showers;time in run [h];#Delta ln(A_{N}/A_{S})', tmax_h, -0.06, 0.10)
lg = mkleg(0.20, 0.74, 0.60, 0.87); say('\nln(N/S) of the beam cell: 10-min medians minus the run median')
for L in range(1, 4):
    xc, med, run, nn = medians(f0.Get('hNS_%d_0' % (L + 1))); r = med - run; g = graph(xc / H2H, r, COL[L]); g.Draw('P'); lg.AddEntry(g, 'SFIL%d (run median %.3f)' % (L + 1, run), 'p')
    pre = np.isfinite(r) & (xc / H2H < tac0); say('  SFIL%d: run median %.3f, rms %.4f, slope %+.4f /h' % (L + 1, run, np.nanstd(r[pre]), np.polyfit(xc[pre] / H2H, r[pre], 1)[0]))
lg.Draw(); ln3 = ROOT.TLine(tac0, -0.06, tac0, 0.10); ln3.SetLineStyle(2); ln3.Draw(); keep.append(ln3)
cv.cd(2); fr = frame('fr4', 'pedestal of the beam cells - run mean;time in run [h];#Delta pedestal [ADC]', tmax_h, -6, 10)
lg = mkleg(0.19, 0.60, 0.62, 0.87, size=0.030); say('\npedestals (samples 1-3) of the beam cells: 10-min means minus the run mean')
j = 0
for L in range(4):
    for e, en in ((0, 'S'), (1, 'N')):
        k = 8 * L + 4 * e + 2; pp = f0.Get('pPed_%d' % k); x = np.array([pp.GetBinCenter(i + 1) for i in range(pp.GetNbinsX())]); y = np.array([pp.GetBinContent(i + 1) if pp.GetBinEntries(i + 1) > 100 else np.nan for i in range(pp.GetNbinsX())])
        if not np.isfinite(y).any(): continue
        run = np.nanmean(y); g = graph(x / H2H, y - run, COL[L], 20 if e == 0 else 24); g.Draw('P'); lg.AddEntry(g, 'SFIL%d %s (%.0f ADC)' % (L + 1, en, run), 'p'); j += 1
        say('  SFIL%d %s cell 3: run mean %.1f ADC, range %+.1f / %+.1f' % (L + 1, en, run, np.nanmin(y - run), np.nanmax(y - run)))
lg.Draw(); ln4 = ROOT.TLine(tac0, -6, tac0, 10); ln4.SetLineStyle(2); ln4.Draw(); keep.append(ln4)
finish(cv); cv.SaveAs('gain_stability_ns_ped.png')

# ---------------------------------------------------------------- clipped fraction
cv = ROOT.TCanvas('c3', '', 900, 650); keep.append(cv)
fr = frame('fr5', 'South beam cells: fraction of fitted pulses at the 1023 ceiling;time in run [h];clipped fraction', tmax_h, 0, 0.9)
lg = mkleg(0.20, 0.74, 0.60, 0.87, ncol=2); say('\nclipped fraction of the South beam cells (run mean)')
for L in range(4):
    pp = f0.Get('pSat_%d' % (8 * L + 2)); x = np.array([pp.GetBinCenter(i + 1) for i in range(pp.GetNbinsX())]); y = np.array([pp.GetBinContent(i + 1) if pp.GetBinEntries(i + 1) > 100 else np.nan for i in range(pp.GetNbinsX())])
    g = graph(x / H2H, y, COL[L]); g.Draw('P'); lg.AddEntry(g, 'SFIL%d S cell 3' % (L + 1), 'p'); say('  SFIL%d: %.3f' % (L + 1, np.nanmean(y[x / H2H < tac0])))
lg.Draw(); finish(cv); cv.SaveAs('gain_stability_clip.png')

# ---------------------------------------------------------------- KCU1: 64 single SiPMs, unclipped pulses
cv = ROOT.TCanvas('c4', '', 1500, 650); cv.Divide(2, 1, 0.002, 0.002); keep.append(cv)
cv.cd(1); fr = frame('fr6', 'KCU1 SFIL5 South, unclipped pulses, mean over the channels of each array;time in run [h];median A_{fit} / run median', tmax_h, 0.90, 1.14)
lg = mkleg(0.19, 0.64, 0.60, 0.87, size=0.034); say('\nKCU1 SFIL5 South: per channel 10-min medians / run median, averaged over the array')
R = {}
for k in range(64):
    xc, med, run, nn = medians(f1.Get('h2_%d' % k), nmin=100)
    if np.isfinite(run) and run > 0: R[k] = med / run
for a, an in enumerate('ABCD'):
    ks = [k for k in R if k // 16 == a]
    if not ks: continue
    M = np.nanmean(np.stack([R[k] for k in ks]), axis=0); g = graph(xc / H2H, M, COL[a]); g.Draw('P'); lg.AddEntry(g, 'array %s (%d channels)' % (an, len(ks)), 'p')
    pre = np.isfinite(M) & (xc / H2H < tac0); say('  array %s: %d channels, rms of the mean curve %.2f %%, slope %+.3f %%/h' % (an, len(ks), 100 * np.nanstd(M[pre]), 100 * np.polyfit(xc[pre] / H2H, M[pre], 1)[0]))
lg.Draw(); ln5 = ROOT.TLine(tac0, 0.90, tac0, 1.14); ln5.SetLineStyle(2); ln5.Draw(); keep.append(ln5)
cv.cd(2); fr = frame('fr7', 'KCU1: clipped fraction (nsat > 0), mean over the array;time in run [h];clipped fraction', tmax_h, 0, 0.5)
lg = mkleg(0.20, 0.74, 0.60, 0.87, ncol=2)
for a, an in enumerate('ABCD'):
    ys = []
    for k in range(16 * a, 16 * a + 16):
        pp = f1.Get('pClip_%d' % k); y = np.array([pp.GetBinContent(i + 1) if pp.GetBinEntries(i + 1) > 100 else np.nan for i in range(pp.GetNbinsX())])
        if np.isfinite(y).any(): ys.append(y)
    if not ys: continue
    x = np.array([pp.GetBinCenter(i + 1) for i in range(pp.GetNbinsX())]); M = np.nanmean(np.stack(ys), axis=0); g = graph(x / H2H, M, COL[a]); g.Draw('P'); lg.AddEntry(g, 'array %s' % an, 'p')
lg.Draw(); finish(cv); cv.SaveAs('gain_stability_kcu1.png')
out.close()
