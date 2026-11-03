"""MIP peaks of the SFIL cells from the cosmic_mip.C output (hadd of the chunks of run 3088, or 3091).

    python3 cosmic_mip_fit.py [cosmic.root] [tag]   -> cosmic_mip_<tag>.txt, cosmic_mip_gm_<tag>.png, cosmic_mip_ends_<tag>.png, cosmic_mip_shape_<tag>.png

Muon tag: both ends of a cell fire (A_S and A_N above the threshold; the two ends sit on different ASICs, so their noise is
independent), both peak at the trigger time (samples 7-9, within one sample of each other; the baseline ramps peak at the
window edge) and at least one end has the pulse shape (from the stored waveform: the sample 4 after the peak is below 0.8 of the peak,
the sample 2 before it below 0.5). Per cell the spectrum of sqrt(A_S A_N) (attenuation-free, the light at the
middle of the fibre) and of each end are fitted with a Landau convolved with a Gaussian; the MPV is the ADC per MIP at the
fibre middle. The noise-only expectation of the tag is measured from the late window (Al, samples 13-19, no pulses).
"""
import sys, os, numpy as np, uproot, ROOT
sys.path.insert(0, '/sphenix/user/tmengel/eic/bic_testbeam/common'); from plotstyle import setup, finish, legend, headroom
setup()
F = sys.argv[1] if len(sys.argv) > 1 else '/sphenix/tg/tg01/jets/tmengel/eic/kcu0/cosmics/cosmic_3088.root'
TAG = sys.argv[2] if len(sys.argv) > 2 else '3088'
THR = 8.; TAIL = 0.8
out = open('cosmic_mip_%s.txt' % TAG, 'w')
def say(*a): s = ' '.join(str(x) for x in a); print(s); out.write(s + '\n')
t = uproot.open(F)['cos']; n = t.num_entries
A = np.zeros((n, 36), np.float32); Al = np.zeros((n, 36), np.float32); PK = np.zeros((n, 36), np.int8); shape = np.zeros((n, 36), bool); nstored = 0
wsum = np.zeros(20); wn = 0                                    # mean waveform of the tagged SFIL3 cell 3 South muons (illustration)
i0 = 0
for arr in t.iterate(['A', 'Al', 'pk', 'nwf', 'wf'], step_size=40000, library='ak'):
    import awkward as ak
    m = len(arr); a = ak.to_numpy(arr['A'])[:, :36]; al = ak.to_numpy(arr['Al'])[:, :36]; pkk = ak.to_numpy(arr['pk'])[:, :36].astype(int); nw = ak.to_numpy(arr['nwf'])
    A[i0:i0 + m] = a; Al[i0:i0 + m] = al; PK[i0:i0 + m] = pkk
    has = nw > 0; nstored += int(has.sum())
    if has.any():
        wf = ak.to_numpy(arr['wf'][has]).reshape(int(has.sum()), 36, 20) / 10.
        pc = np.clip(pkk[has], 5, 11); tail = np.take_along_axis(wf, np.clip(pc + 4, 0, 19)[:, :, None], axis=2)[:, :, 0]; rise = np.take_along_axis(wf, np.clip(pc - 2, 0, 19)[:, :, None], axis=2)[:, :, 0]
        sh = (tail < TAIL * a[has]) & (rise < 0.5 * a[has]); shape[i0:i0 + m][has] = sh
        mu = (a[has][:, 18] > THR) & (a[has][:, 22] > THR) & (sh[:, 18] | sh[:, 22]) & (pkk[has][:, 18] >= 7) & (pkk[has][:, 18] <= 9) & (pkk[has][:, 22] >= 7) & (pkk[has][:, 22] <= 9)
        if mu.any(): wsum += wf[mu][:, 18, :].sum(axis=0); wn += int(mu.sum())
    i0 += m
say('%s: %d triggers, %d with stored waveforms' % (F, n, nstored))
def kS(L, c): return 8 * L + c
def kN(L, c): return 8 * L + 4 + c
keep = []
# Landau x Gauss via numerical convolution (the classic langaufun)
ROOT.gROOT.ProcessLine('''
double langaufun(double *x, double *par) {
  // par: width (Landau scale), MP, area, gaussian sigma
  double invsq2pi = 0.3989422804014; double mpshift = -0.22278298; int np = 100; double sc = 5.0;
  double mpc = par[1] - mpshift * par[0]; double xlow = x[0] - sc * par[3], xupp = x[0] + sc * par[3]; double step = (xupp - xlow) / np; double sum = 0;
  for (int i = 1; i <= np / 2; i++) { double xx = xlow + (i - .5) * step; sum += TMath::Landau(xx, mpc, par[0]) / par[0] * TMath::Gaus(x[0], xx, par[3]); xx = xupp - (i - .5) * step; sum += TMath::Landau(xx, mpc, par[0]) / par[0] * TMath::Gaus(x[0], xx, par[3]); }
  return par[2] * step * sum * invsq2pi / par[3]; }
''')
def fit_langaus(h, lo, hi, name):
    f = ROOT.TF1(name, ROOT.langaufun, lo, hi, 4); keep.append(f)
    mx = h.GetBinCenter(h.GetMaximumBin()); f.SetParameters(0.12 * mx, mx, h.Integral() * h.GetBinWidth(1), 0.15 * mx)
    f.SetParLimits(0, 0.3, 50); f.SetParLimits(1, 3, 200); f.SetParLimits(3, 0.5, 50); f.SetLineColor(ROOT.kRed + 1); f.SetLineWidth(2)
    r = h.Fit(f, 'RQLS0'); ok = int(r) == 0 and f.GetParError(1) > 0
    return f, (f.GetParameter(1), f.GetParError(1), f.GetParameter(0), f.GetParameter(3), f.GetChisquare() / max(1, f.GetNDF())) if ok else None
def H(name, title, vals, nb=60, lo=0, hi=120):
    h = ROOT.TH1D(name, title, nb, lo, hi); keep.append(h)
    for v in vals: h.Fill(float(v))
    return h
say('\nper cell: muons (both ends > %.0f, peak at samples 7-9 within one sample of each other, pulse shape at one end at least), accidental expectation for independent ends, MPV of sqrt(A_S A_N), A_S, A_N [ADC] (Landau width, Gauss sigma, chi2/ndf)' % THR)
say('%-6s %7s %7s | %-22s | %-22s | %-22s' % ('cell', 'muons', 'accid.', 'GM: MPV +- err (w, s, chi2)', 'South', 'North'))
res = {}
cv = ROOT.TCanvas('cgm', '', 2000, 1500); cv.Divide(4, 4, 0.002, 0.002); keep.append(cv)
cv2 = ROOT.TCanvas('cends', '', 2000, 1500); cv2.Divide(4, 4, 0.002, 0.002); keep.append(cv2)
for L in range(4):
    for c in range(4):
        s, nn = kS(L, c), kN(L, c)
        # muon tag: both ends above threshold, both peaking at the trigger time (samples 7-9) within one sample of each other, both with the pulse shape
        pS = (A[:, s] > THR) & (PK[:, s] >= 7) & (PK[:, s] <= 9); pN = (A[:, nn] > THR) & (PK[:, nn] >= 7) & (PK[:, nn] <= 9)
        m = pS & pN & (np.abs(PK[:, s].astype(int) - PK[:, nn].astype(int)) <= 1) & (shape[:, s] | shape[:, nn])
        mn = np.array([pS.mean() * pN.mean() * n])                      # accidental expectation if the two ends were independent (upper limit of the noise in the tag)
        gm = np.sqrt(A[m, s] * A[m, nn])
        pad = (L, c)
        hg = H('hgm_%d_%d' % (L, c), 'SFIL%d cell %d: #sqrt{A_{S}A_{N}}, %d muons;#sqrt{A_{S}A_{N}} [ADC];muons' % (L + 1, c + 1, m.sum()), gm)
        hs = H('hs_%d_%d' % (L, c), 'SFIL%d cell %d South;A_{S} [ADC];muons' % (L + 1, c + 1), A[m, s]); hn = H('hn_%d_%d' % (L, c), 'SFIL%d cell %d North;A_{N} [ADC];muons' % (L + 1, c + 1), A[m, nn])
        rg = rs = rn = None
        if m.sum() >= 40:
            fg, rg = fit_langaus(hg, 12, 120, 'fg_%d_%d' % (L, c))
            fs, rs = fit_langaus(hs, 12, 120, 'fs_%d_%d' % (L, c)); fn, rn = fit_langaus(hn, 12, 120, 'fn_%d_%d' % (L, c))
        fmt = lambda r: '%5.1f +- %4.1f (%4.1f %4.1f %4.1f)' % r if r and r[1] < 0.25 * r[0] else '      -               '
        q = np.percentile(gm, [25, 50, 75]) if m.sum() else (0, 0, 0)
        say('S%dc%d %7d %7.1f | %s | %s | %s | GM median %5.1f [%4.1f, %4.1f]' % (L + 1, c + 1, m.sum(), float(mn[0]), fmt(rg), fmt(rs), fmt(rn), q[1], q[0], q[2]))
        res[pad] = (m.sum(), float(mn[0]), rg, rs, rn)
        cv.cd(4 * L + c + 1); hg.SetLineColor(ROOT.kBlack); hg.SetLineWidth(2); headroom(hg, 0, hg.GetMaximum(), 0.3); hg.Draw('hist')
        tx = ROOT.TLatex(); tx.SetNDC(); tx.SetTextSize(0.06); keep.append(tx)
        if rg and rg[1] < 0.25 * rg[0]: fg.Draw('same'); tx.DrawLatex(0.40, 0.78, 'MPV %.1f #pm %.1f ADC' % rg[:2])
        elif m.sum(): tx.DrawLatex(0.40, 0.78, 'median %.1f ADC' % q[1])
        cv2.cd(4 * L + c + 1); hs.SetLineColor(ROOT.kBlue + 1); hn.SetLineColor(ROOT.kRed + 1); hs.SetLineWidth(2); hn.SetLineWidth(2); headroom(hs, 0, max(hs.GetMaximum(), hn.GetMaximum()), 0.3); hs.Draw('hist'); hn.Draw('hist same')
        lg = legend(0.45, 0.70, 0.92, 0.86); keep.append(lg); lg.AddEntry(hs, 'South' + (' MPV %.1f' % rs[0] if rs else ''), 'l'); lg.AddEntry(hn, 'North' + (' MPV %.1f' % rn[0] if rn else ''), 'l'); lg.Draw()
finish(cv); cv.SaveAs('cosmic_mip_gm_%s.png' % TAG); finish(cv2); cv2.SaveAs('cosmic_mip_ends_%s.png' % TAG)
# in-board tag (the other two of cells 1-3 of the same end and layer; cell 4 by cells 2-3). WARNING (2026-10-30): on SFIL1 S, SFIL2 S
# and SFIL5 N this selects within-board bursts, not muons (corr(A_c1, A_c2) = 0.97-0.99, ratio 0.91, thousands of events where SFIL3
# has 4-26); the table is kept as a monitor of that noise, its MPVs are not MIP values
say('\nin-board tag (other two cells of 1-3 > %.0f at samples 7-9 with the shape; cell 4 by cells 2-3): tagged, accidental expectation, MPV or median [ADC]. NOT muons on SFIL1 S, SFIL2 S, SFIL5 N: within-board bursts (corr 0.97-0.99)' % THR)
good = lambda k: (A[:, k] > THR) & (PK[:, k] >= 7) & (PK[:, k] <= 9) & shape[:, k]
cv5 = ROOT.TCanvas('cvert', '', 2000, 1500); cv5.Divide(4, 4, 0.002, 0.002); keep.append(cv5)
for L in range(4):
    for e, en in ((0, 'S'), (1, 'N')):
        line = '  SFIL%d %s:' % (L + 1, en)
        for c in range(4):
            k = 8 * L + 4 * e + c; others = [8 * L + 4 * e + cc for cc in (0, 1, 2) if cc != c] if c < 3 else [8 * L + 4 * e + 1, 8 * L + 4 * e + 2]
            tag = good(others[0]) & good(others[1]); m = tag & (PK[:, k] >= 7) & (PK[:, k] <= 9)
            acc = good(others[0]).mean() * good(others[1]).mean() * n
            vals = A[m, k]; q = np.percentile(vals, [25, 50, 75]) if m.sum() else (0, 0, 0)
            r = None
            if e == 0:
                h = H('hv_%d_%d' % (L, c), 'SFIL%d cell %d South, vertical tag (%d);A_{S} [ADC];muons' % (L + 1, c + 1, m.sum()), vals, 60, -20, 100)
                if m.sum() >= 40: f, r = fit_langaus(h, 6, 100, 'fv_%d_%d' % (L, c))
                cv5.cd(4 * L + c + 1); h.SetLineWidth(2); headroom(h, 0, h.GetMaximum(), 0.3); h.Draw('hist')
                tx = ROOT.TLatex(); tx.SetNDC(); tx.SetTextSize(0.06); keep.append(tx)
                if r and r[1] < 0.25 * r[0]: f.Draw('same'); tx.DrawLatex(0.40, 0.78, 'MPV %.1f #pm %.1f' % r[:2])
                elif m.sum(): tx.DrawLatex(0.40, 0.78, 'median %.1f' % q[1])
            line += '  c%d %5d (%5.0f) %s med %5.1f |' % (c + 1, m.sum(), acc, ('MPV %5.1f+-%4.1f' % r[:2]) if r and r[1] < 0.25 * r[0] else '              ', q[1])
        say(line)
finish(cv5); cv5.SaveAs('cosmic_mip_vertical_%s.png' % TAG)
# SFIL5 North (no South partner): vertical tag from the other SFIL5 N cells, A with the shape cut
say('\nSFIL5 North (KCU0 ch 153/156/163/167), in-board tag by the other two of cells 1-3 > %.0f with the shape cut - within-board bursts dominate (corr 0.70), not a MIP value:' % THR)
cv3 = ROOT.TCanvas('c5', '', 2000, 500); cv3.Divide(4, 1, 0.002, 0.002); keep.append(cv3)
for c in range(4):
    k = 32 + c; others = [32 + cc for cc in (0, 1, 2) if cc != c] if c < 3 else [33, 34]
    m = np.all([(A[:, o] > THR) & shape[:, o] & (PK[:, o] >= 7) & (PK[:, o] <= 9) for o in others], axis=0) & shape[:, k] & (PK[:, k] >= 7) & (PK[:, k] <= 9)
    h = H('h5_%d' % c, 'SFIL5 North cell %d, %d tagged;A [ADC];muons' % (c + 1, m.sum()), A[m, k], 60, -20, 100)
    r = None
    if m.sum() >= 40: f, r = fit_langaus(h, 5, 100, 'f5_%d' % c)
    say('  N5c%d: %d tagged, MPV %s' % (c + 1, m.sum(), '%.1f +- %.1f' % r[:2] if r else '-'))
    cv3.cd(c + 1); h.SetLineWidth(2); headroom(h, 0, h.GetMaximum(), 0.3); h.Draw('hist')
    if r: f.Draw('same')
finish(cv3); cv3.SaveAs('cosmic_mip_sfil5n_%s.png' % TAG)
# shape-cut illustration: mean corrected waveform of the tagged SFIL3 cell 3 South muons
cv4 = ROOT.TCanvas('c4', '', 900, 650); keep.append(cv4)
hw = ROOT.TH1D('hw', 'SFIL3 cell 3 South: mean corrected waveform of the tagged muons (%d);sample (25 ns);ADC' % wn, 20, -0.5, 19.5); keep.append(hw)
if wn:
    for i in range(20): hw.SetBinContent(i + 1, wsum[i] / wn)
hw.SetLineWidth(2); hw.Draw('hist'); finish(cv4); cv4.SaveAs('cosmic_mip_shape_%s.png' % TAG)
out.close()
