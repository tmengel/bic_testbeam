"""AstroPix ToT calibration of run 3062 along the lines of the FTBF paper (arXiv:2602.06084):
  1. pixel hits = row record + column record of the same chip with |dToA| <= 1 (8-bit 10 MHz on-chip counter) and
     |ToT_col - ToT_row| / ToT_col < 10 %, within the same RCDAQ event and a few records of each other;
  2. per-pixel ToT spectra of the beam hits, Landau (x Gauss at chip level) fits -> ToT most probable value per pixel;
  3. the paper converted ToT to keV with per-pixel 241Am / 133Ba photopeak curves and found a MIP MPV of 34.0 +- 1.9 keV
     (128 um effective depletion) for 120 GeV protons.  No source data exist for this module here, so the inverse is done:
     the per-pixel MPV map is the relative calibration, and the absolute scale is anchored at the expected MIP MPV, which for
     3-6 GeV positrons in 128 um of silicon is the same 34 keV (Fermi plateau; computed below from the Landau-Vavilov MPV with
     the Sternheimer density effect).
    python3 astro_tot_calib.py  -> astro_tot_calib.txt, astro_tot_calib.png, astro_tot_pixels.npz
"""
import numpy as np, uproot, sys, ROOT
sys.path.insert(0, '/sphenix/user/tmengel/eic/bic_testbeam/common'); from plotstyle import setup, legend, label
setup(); ROOT.gROOT.SetBatch(True); out = open('astro_tot_calib.txt', 'w')
def say(*a): print(*a, flush=True); print(*a, file=out)
# ---- expected MIP MPV in silicon (Landau-Vavilov most probable value, PDG eq., Sternheimer delta for Si)
def mpv_si(bg, x_um):
    K, ZA, rho, I, me = 0.307075, 0.49848, 2.329, 173e-6, 0.510999   # MeV cm2/g, -, g/cm3, MeV, MeV
    b2 = bg ** 2 / (1 + bg ** 2); x = x_um * 1e-4; xi = K / 2 * ZA * rho * x / b2
    X = np.log10(bg); X0, X1, a, k, C = 0.2015, 2.8716, 0.1492, 3.2546, 4.4355
    delta = 0 if X < X0 else (2 * np.log(10) * X - C + (a * (X1 - X) ** k if X < X1 else 0))
    return 1e3 * xi * (np.log(2 * me * bg ** 2 / I) + np.log(xi / I) + 0.2 - b2 - delta)
say('expected MIP most-probable energy loss in 128 um Si: 120 GeV protons (beta gamma 128) %.1f keV; 4 GeV positrons (beta gamma 7830) %.1f keV; 3 / 6 GeV %.1f / %.1f keV' % (mpv_si(128., 128), mpv_si(4000. / 0.511, 128), mpv_si(3000. / 0.511, 128), mpv_si(6000. / 0.511, 128)))
# ---- records
t = uproot.open('/sphenix/tg/tg01/jets/tmengel/eic/astropix/astropix_3062.root')['hits']
a = t.arrays(['id', 'iscol', 'location', 'ts', 'tot', 'tot_nibble', 'valid', 'evt', 'layer', 'payload'], library='np')
v = a['valid'] & (a['id'] <= 8) & (a['location'] <= 34) & (a['layer'] == 1) & (a['payload'] == 4) & (a['tot_nibble'] == 0)
ids, isc, loc, ts, tot, ev = a['id'][v].astype(int), a['iscol'][v].astype(bool), a['location'][v].astype(int), a['ts'][v].astype(int), a['tot'][v].astype(int), a['evt'][v]
say('records: %d total, %d clean (valid, chip 0-8, location 0-34, upper nibble 0)' % (len(v), v.sum()))
# ---- pairing: row record i with the nearest column record j (same chip, same event, within 8 records) with |dts| <= 1 and relative ToT difference < 10 %
n = len(ids); idx = np.arange(n); rows = []; cols = []
used_c = np.zeros(n, bool)
for d in list(range(1, 9)) + list(range(-1, -9, -1)):
    i = idx[max(0, -d): n - max(0, d)]; j = i + d
    m = (~isc[i]) & isc[j] & (ids[i] == ids[j]) & (ev[i] == ev[j]) & ((np.abs(ts[i] - ts[j]) <= 1) | (np.abs(ts[i] - ts[j]) >= 255)) & (tot[j] > 0) & (np.abs(tot[j] - tot[i]) < 0.1 * tot[j]) & ~used_c[j]
    rows.append(i[m]); cols.append(j[m]); used_c[j[m]] = True
ri = np.concatenate(rows); ci = np.concatenate(cols)
# one column per row record: keep the first assignment per row
_, first = np.unique(ri, return_index=True); ri, ci = ri[first], ci[first]
chip, row, col = ids[ri], loc[ri], loc[ci]; totc = tot[ci]; totr = tot[ri]
say('pixel hits (row + column records paired with |dToA| <= 1 and |dToT|/ToT_col < 10 %%): %d; row-only records with ToT > 0 left unpaired: %d' % (len(ri), ((~isc) & (tot > 0)).sum() - len(ri)))
say('per chip: hits, median column ToT, fraction of column records with ToT = 0 (chips 5-8 carry no ToT):')
for c in range(9): say('  chip %d: hits %7d  median ToT %5.0f   zero-ToT column records %.3f' % (c, (chip == c).sum(), np.median(totc[chip == c]) if (chip == c).sum() else -1, (tot[isc & (ids == c)] == 0).mean()))
good_chips = [c for c in range(9) if (tot[isc & (ids == c)] == 0).mean() < 0.1]
say('chips with a working ToT: %s' % good_chips)
# ---- langaus (Landau x Gauss) for chip-level spectra
ROOT.gInterpreter.Declare('''
double langaufun(double *x, double *par) {
  // par: 0 Landau width, 1 MPV, 2 area, 3 Gauss sigma  (ROOT tutorial langaus.C)
  double invsq2pi = 0.3989422804014; double mpshift = -0.22278298; double np = 100.0; double sc = 5.0;
  double xx, mpc, fland, sum = 0.0, xlow, xupp, step, i;
  mpc = par[1] - mpshift * par[0]; xlow = x[0] - sc * par[3]; xupp = x[0] + sc * par[3]; step = (xupp - xlow) / np;
  for (i = 1.0; i <= np / 2; i++) { xx = xlow + (i - .5) * step; fland = TMath::Landau(xx, mpc, par[0]) / par[0]; sum += fland * TMath::Gaus(x[0], xx, par[3]);
    xx = xupp - (i - .5) * step; fland = TMath::Landau(xx, mpc, par[0]) / par[0]; sum += fland * TMath::Gaus(x[0], xx, par[3]); }
  return (par[2] * step * sum * invsq2pi / par[3]); }
''')
c1 = ROOT.TCanvas('c1', '', 2400, 1600); c1.Divide(3, 2); keep = []
chipmpv = {}
say('\nchip-level ToT spectra (beam band rows 10-24, columns >= 3), Landau x Gauss fits: MPV, Landau width, Gauss sigma [ToT units; x 5 ns if the ToT clock is 200 MHz]')
for p, c in enumerate(good_chips[:5]):
    m = (chip == c) & (row >= 10) & (row <= 24) & (col >= 3); v_ = totc[m]
    c1.cd(p + 1); h = ROOT.TH1F('hc%d' % c, 'chip %d;column ToT [units];pixel hits' % c, 100, 0, 4000); [h.Fill(x) for x in v_]; h.SetLineWidth(2); h.Draw('hist'); keep.append(h)
    f = ROOT.TF1('lg%d' % c, ROOT.langaufun, 100, 4000, 4); x0 = h.GetBinCenter(h.GetMaximumBin()); f.SetParameters(0.12 * x0, x0, h.GetEntries() * h.GetBinWidth(1), 0.1 * x0); f.SetParLimits(0, 10, 2000); f.SetParLimits(3, 5, 2000)
    h.Fit(f, 'QRN', '', 0.4 * x0, 3.0 * x0); f.SetLineColor(ROOT.kRed + 1); f.Draw('same'); keep.append(f); chipmpv[c] = f.GetParameter(1)
    say('  chip %d: hits %6d, mode %4.0f, MPV %6.1f +- %4.1f, width %5.1f, sigma %5.1f, chi2/ndf %.1f' % (c, m.sum(), x0, f.GetParameter(1), f.GetParError(1), f.GetParameter(0), f.GetParameter(3), f.GetChisquare() / max(1, f.GetNDF())))
    keep.append(label('MPV %.0f' % f.GetParameter(1), x=0.55, y=0.8, size=0.05))
# ---- per-pixel MPV (Landau fit around the mode) for pixels with >= 150 hits
pix = chip * 10000 + row * 100 + col; up, cnt = np.unique(pix, return_counts=True); sel = up[cnt >= 150]
say('\npixels with >= 150 paired hits: %d (of %d hit pixels); per-pixel Landau fits' % (len(sel), len(up)))
fl = ROOT.TF1('fl', 'landau', 0, 4000); PM = []; order = np.argsort(pix); ps = pix[order]; ts_ = totc[order]
starts = np.searchsorted(ps, sel); ends = np.searchsorted(ps, sel, side='right')
for k, (s0, s1) in enumerate(zip(starts, ends)):
    vv = ts_[s0:s1]; h = ROOT.TH1F('hp', '', 60, 0, 3000); [h.Fill(x) for x in vv]; x0 = h.GetBinCenter(h.GetMaximumBin()); fl.SetParameters(h.GetMaximum(), x0, 0.15 * x0)
    r = h.Fit(fl, 'QRNS', '', 0.5 * x0, 2.0 * x0); PM.append((sel[k] // 10000, (sel[k] // 100) % 100, sel[k] % 100, s1 - s0, fl.GetParameter(1), fl.GetParError(1), fl.GetParameter(2))); h.Delete()
PM = np.array(PM); np.savez('astro_tot_pixels.npz', pm=PM, chipmpv=np.array([(c, chipmpv[c]) for c in chipmpv]))
for c in good_chips:
    m = PM[:, 0] == c
    if m.sum() < 3: continue
    say('  chip %d: %3d pixels, per-pixel MPV median %5.0f, 16-84 %% %5.0f-%5.0f, rms/median %.3f, median fit error %.1f %%' % (c, m.sum(), np.median(PM[m, 4]), np.percentile(PM[m, 4], 16), np.percentile(PM[m, 4], 84), PM[m, 4].std() / np.median(PM[m, 4]), 100 * np.median(PM[m, 5] / PM[m, 4])))
allm = np.median(PM[:, 4]); say('all fitted pixels: %d, MPV median %.0f, rms/median %.3f (the paper: 4.7 %% pixel-to-pixel spread in keV after per-pixel calibration)' % (len(PM), allm, PM[:, 4].std() / allm))
say('\none-point energy scale: E = 34.0 keV x ToT / MPV_pixel (relative per pixel; absolute anchored to the FTBF MIP value, +-6 %%; the ToT-energy curve is not linear through zero, threshold ~20 keV, so this holds near the MIP only). In ToT units per keV at the chip level: ' + ', '.join('chip %d %.1f' % (c, chipmpv[c] / 34.0) for c in chipmpv))
# ---- figure part 2: pixel MPV map and distribution
c1.cd(6); fr = ROOT.gPad.DrawFrame(0, 0, 3000, max(10, np.histogram(PM[:, 4], bins=60, range=(0, 3000))[0].max() * 1.3), ';per-pixel ToT MPV [units];pixels'); keep.append(fr)
hm = ROOT.TH1F('hm', '', 60, 0, 3000); [hm.Fill(x) for x in PM[:, 4]]; hm.SetLineWidth(2); hm.Draw('hist same'); keep.append(hm); keep.append(label('pixels with >= 150 hits, chips %s' % good_chips, x=0.15, y=0.93, size=0.035))
c1.SaveAs('astro_tot_calib.png'); c1.SaveAs('astro_tot_calib.pdf')
c2 = ROOT.TCanvas('c2', '', 2400, 600); c2.Divide(len(good_chips), 1); keep2 = []
for p, c in enumerate(good_chips):
    c2.cd(p + 1); h2 = ROOT.TH2F('map%d' % c, 'chip %d: ToT MPV / chip median;column;row' % c, 35, -0.5, 34.5, 35, -0.5, 34.5); m = PM[:, 0] == c; med = np.median(PM[m, 4]) if m.sum() else 1
    for r_ in PM[m]: h2.SetBinContent(int(r_[2]) + 1, int(r_[1]) + 1, r_[4] / med)
    h2.SetMinimum(0.6); h2.SetMaximum(1.4); h2.Draw('colz'); keep2.append(h2)
c2.SaveAs('astro_tot_mpv_map.png'); say('wrote astro_tot_calib.png, astro_tot_mpv_map.png, astro_tot_pixels.npz')
