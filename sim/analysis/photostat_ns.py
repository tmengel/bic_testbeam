"""Route 3: ADC per photoelectron of the SFIL2-4 channels from the North/South photostatistics of run 3062, no simulation.

For a cell read at both ends, S and N see the same scintillation light L up to the (fixed) attenuation; after
equalising N to S by the median ratio r, x = S and y = N/r are two measurements of L with independent fluctuations:
    var(x) - cov(x, y) = sigma_S^2,   var(y) - cov(x, y) = sigma_N^2 / r^2        (per bin of s = (x + y) / 2)
and sigma^2 = c + a A + b A^2 with a = g F (ADC per photoelectron times the excess-noise factor), c the electronic
noise + fit floor, b the amplitude-proportional terms (template shape error, beam-position attenuation spread).
Equal photoelectron yields at both ends would make the amplitude ratio r equal the gain ratio a_N / a_S; a ratio r
smaller than a_N / a_S means fewer photoelectrons at the North end (light), not a smaller gain.
Pairs: template-fit amplitudes of both ends accepted and unsaturated (peak sample below 1020 ADC).
    python3 photostat_ns.py  -> photostat_pairs.npz (cache), photostat_ns.txt, photostat_ns.png
"""
import os, sys, numpy as np, uproot, ROOT
sys.path.insert(0, '/sphenix/user/tmengel/eic/bic_testbeam/common'); from plotstyle import setup, legend, label
F = '/sphenix/tg/tg01/jets/tmengel/eic/kcu0/fit_full/kcu0_sfil_fit_full.root'
PAIRS = [(L, c) for L in (1, 2, 3) for c in range(4)]      # SFIL2-4 (SFIL1 North is masked), cells 1-4
if not os.path.exists('photostat_pairs.npz'):
    t = uproot.open(F)['fits']; acc = {p: [] for p in PAIRS}; pads = {p: [] for p in PAIRS}
    for a in t.iterate(['pad', 'ok', 'amp', 'afit', 'ped'], step_size=200000, library='np'):
        ok = a['ok']; raw = np.nan_to_num(a['amp']); fit = np.nan_to_num(a['afit']); ped = np.nan_to_num(a['ped']); unsat = raw + ped < 1020.
        for (L, c) in PAIRS:
            kS, kN = 8 * L + c, 8 * L + 4 + c; m = ok[:, kS] & ok[:, kN] & unsat[:, kS] & unsat[:, kN]
            acc[(L, c)].append(np.stack([fit[m, kS], fit[m, kN]], axis=1)); pads[(L, c)].append(a['pad'][m])
        print('.', end='', flush=True)
    np.savez('photostat_pairs.npz', **{'p%d_%d' % p: np.concatenate(acc[p]) for p in PAIRS}, **{'pad%d_%d' % p: np.concatenate(pads[p]) for p in PAIRS}); print()
P = np.load('photostat_pairs.npz')
out = open('photostat_ns.txt', 'w')
def say(*a): print(*a); print(*a, file=out)
say('pairs with both ends fitted and unsaturated:', {'SFIL%d c%d' % (L + 1, c + 1): len(P['p%d_%d' % (L, c)]) for (L, c) in PAIRS})

def analyse(xy, edges, label_=''):
    x, N = xy[:, 0], xy[:, 1]; big = (x + N) > 200
    r = np.median(N[big] / x[big]) if big.sum() > 50 else np.median(N / np.maximum(x, 1)); y = N / r; s = 0.5 * (x + y)
    rows = []
    for lo, hi in zip(edges[:-1], edges[1:]):
        m = (s >= lo) & (s < hi)
        if m.sum() < 200: continue
        xs, ys = x[m], y[m]; d = xs - ys; med = np.median(d); mad = 1.4826 * np.median(np.abs(d - med)); keep = np.abs(d - med) < 5 * max(mad, 1.)   # pile-up outliers
        xs, ys = xs[keep], ys[keep]; n = len(xs)
        vx, vy, cxy = xs.var(ddof=1), ys.var(ddof=1), np.cov(xs, ys)[0, 1]
        rows.append((xs.mean(), ys.mean(), n, vx - cxy, vy - cxy, cxy))
    rows = np.array(rows)
    def fit(A, var, n):
        w = (n - 1) / (2 * np.maximum(var, 1.) ** 2)      # 1/variance of a variance estimate
        M = np.stack([np.ones_like(A), A, A ** 2], axis=1); W = np.sqrt(w)[:, None]
        coef, *_ = np.linalg.lstsq(M * W, var * W[:, 0], rcond=None); cov = np.linalg.inv((M * W).T @ (M * W))
        return coef, np.sqrt(np.diag(cov))
    use = rows[:, 0] < (620. if edges[-1] < 1000 else 900.)   # below the saturation edge (the unsaturated requirement truncates the top bins)
    cS, eS = fit(rows[use, 0], rows[use, 3], rows[use, 2]); cN, eN = fit(rows[use, 1], rows[use, 4] * r * r, rows[use, 2])   # sigma_N^2 in N's own ADC, vs N amplitude
    return dict(r=r, rows=rows, aS=cS[1], aS_err=eS[1], aN=cN[1], aN_err=eN[1], cS=cS, cN=cN, n=len(x))

edges = np.array([60, 90, 120, 160, 200, 250, 300, 360, 430, 500, 580, 670, 770, 880])
res = {}
say('\nper channel pair: r = median N/S amplitude ratio; a = linear coefficient of sigma^2 vs amplitude = ADC per photoelectron x excess noise (South, North); a_N/a_S')
say('%-12s %8s %6s %10s %10s %10s %12s %14s' % ('pair', 'pairs', 'r', 'a_S', 'a_N', 'a_N/a_S', 'r/(a_N/a_S)', 'c_S, b_S'))
for (L, c) in PAIRS:
    xy = P['p%d_%d' % (L, c)]
    if len(xy) < 2000: say('%-12s %8d  too few pairs' % ('SFIL%d c%d' % (L + 1, c + 1), len(xy))); continue
    R = analyse(xy, edges); res[(L, c)] = R; ratio = R['aN'] / R['aS']; eratio = ratio * np.hypot(R['aS_err'] / R['aS'], R['aN_err'] / R['aN'])
    say('%-12s %8d %6.3f %5.2f+-%.2f %5.2f+-%.2f %6.3f+-%.3f %12.3f %8.1f, %.5f' % ('SFIL%d c%d' % (L + 1, c + 1), R['n'], R['r'], R['aS'], R['aS_err'], R['aN'], R['aN_err'], ratio, eratio, R['r'] / ratio, R['cS'][0], R['cS'][2]))
say('\n(r / (a_N/a_S) = ratio of photoelectron yields North/South: 1 = the amplitude difference is all gain; r itself = all light)')
say('\nbeam-row (cell 3) sigma^2 table per bin: A_S, n, sigma_S^2, sigma_N^2 (in N ADC):')
for (L, c) in PAIRS:
    if c != 2 or (L, c) not in res: continue
    R = res[(L, c)]; say(' SFIL%d:' % (L + 1) + ' '.join('[%.0f n%d %.0f %.0f]' % (row[0], row[2], row[3], row[4] * R['r'] ** 2) for row in R['rows']))
say('\nLayer sums (all four cells per end, cache data_fit_layers.npz; triggers with no saturated channel): the cell sharing cancels')
D = np.load('data_fit_layers.npz'); S, N = D['S'].astype(float), D['N'].astype(float); clean = (D['nsatfail'] == 0) & (D['Ssat'].sum(axis=1) == 0) & (D['Nsat'].sum(axis=1) == 0)
say('%-12s %8s %6s %10s %10s %10s %12s' % ('layer', 'events', 'r', 'a_S', 'a_N', 'a_N/a_S', 'r/(a_N/a_S)'))
edgesL = np.array([100, 150, 200, 260, 330, 410, 500, 600, 720, 860, 1020, 1200, 1400])
resL = {}
for L in (1, 2, 3):
    m = clean & (S[:, L] > 60) & (N[:, L] > 60); xy = np.stack([S[m, L], N[m, L]], axis=1)
    R = analyse(xy, edgesL); resL[L] = R; ratio = R['aN'] / R['aS']; eratio = ratio * np.hypot(R['aS_err'] / R['aS'], R['aN_err'] / R['aN'])
    say('%-12s %8d %6.3f %5.2f+-%.2f %5.2f+-%.2f %6.3f+-%.3f %12.3f' % ('SFIL%d sum' % (L + 1), R['n'], R['r'], R['aS'], R['aS_err'], R['aN'], R['aN_err'], ratio, eratio, R['r'] / ratio))
    say('   bins: ' + ' '.join('[%.0f n%d %.0f %.0f]' % (row[0], row[2], row[3], row[4] * R['r'] ** 2) for row in R['rows']))
out.close()

# ---- linear fits sigma^2 = c + a A over 100 < A < 450 ADC (below saturation, no quadratic term): the quoted result
edgesLin = np.array([100, 130, 160, 200, 240, 290, 340, 400, 450])
def linfit(xy):
    x, N = xy[:, 0], xy[:, 1]; big = (x + N) > 200; r = np.median(N[big] / x[big]); y = N / r; s = 0.5 * (x + y); o = {}
    for end, (u, v, scale) in {'S': (x, y, 1.), 'N': (y, x, r)}.items():
        A, V, W = [], [], []
        for lo, hi in zip(edgesLin[:-1], edgesLin[1:]):
            m = (s >= lo) & (s < hi)
            if m.sum() < 300: continue
            d = u[m] - v[m]; med = np.median(d); mad = 1.4826 * np.median(np.abs(d - med)); k = np.abs(d - med) < 5 * max(mad, 1.)
            uu, vv = u[m][k], v[m][k]; var = uu.var(ddof=1) - np.cov(uu, vv)[0, 1]
            A.append(uu.mean() * scale); V.append(var * scale * scale); W.append((len(uu) - 1) / (2 * max(var * scale * scale, 1.) ** 2))
        A, V, W = map(np.array, (A, V, W)); M = np.stack([np.ones_like(A), A], axis=1) * np.sqrt(W)[:, None]
        coef, *_ = np.linalg.lstsq(M, V * np.sqrt(W), rcond=None); err = np.sqrt(np.diag(np.linalg.inv(M.T @ M))); o[end] = dict(a=coef[1], a_err=err[1], c=coef[0], A=A, V=V, n=len(A))
    return r, o
LIN = {}
out = open('photostat_ns.txt', 'a'); say('\n=== RESULT: linear fits over 100-450 ADC')
say('%-10s %8s %6s %14s %14s %9s %12s' % ('pair', 'pairs', 'N/S', 'a_S [ADC/pe]', 'a_N [ADC/pe]', 'a_N/a_S', 'pe_N/pe_S'))
for (L, c) in PAIRS:
    xy = P['p%d_%d' % (L, c)]
    if len(xy) < 3000: continue
    r, o = linfit(xy); LIN[(L, c)] = (r, o); ratio = o['N']['a'] / o['S']['a']
    say('%-10s %8d %6.3f %6.2f +- %4.2f %6.2f +- %4.2f %9.3f %12.3f' % ('SFIL%d c%d' % (L + 1, c + 1), len(xy), r, o['S']['a'], o['S']['a_err'], o['N']['a'], o['N']['a_err'], ratio, r / ratio))
say('SFIL4: not photostatistics (sigma^2 quadratic, North cells 1-2 carry a constant ~130 ADC and cell 3 only 0.35 of South): no gain from this layer')
out.close()

# ---- figure: sigma^2 vs amplitude for the beam-row pairs, South and North, with the linear fits
setup(); ROOT.gStyle.SetPadLeftMargin(0.17); ROOT.gStyle.SetPadRightMargin(0.04); ROOT.gStyle.SetPadTopMargin(0.08); ROOT.gStyle.SetPadBottomMargin(0.15); ROOT.gROOT.ForceStyle()
c = ROOT.TCanvas('c', '', 1500, 520); c.Divide(3, 1, 0.002, 0.002); keep = []
for i, L in enumerate((1, 2, 3)):
    if (L, 2) not in LIN: continue
    r, o = LIN[(L, 2)]; c.cd(i + 1); ROOT.gPad.SetTopMargin(0.1)
    ymax = max(o['S']['V'].max(), o['N']['V'].max()) * 2.2
    fr = ROOT.gPad.DrawFrame(0, 0, 520, ymax, ';fitted amplitude [ADC];#sigma^{2} of one end [ADC^{2}]'); fr.GetYaxis().SetTitleOffset(1.7); fr.GetXaxis().SetNdivisions(505); fr.GetYaxis().SetNdivisions(505); fr.GetYaxis().SetMaxDigits(4); keep.append(fr)
    for end, col, mk in (('S', ROOT.kBlack, 20), ('N', ROOT.kRed + 1, 25)):
        A, V = o[end]['A'].astype(float), o[end]['V'].astype(float); g = ROOT.TGraph(len(A), A, V); g.SetMarkerStyle(mk); g.SetMarkerSize(1.4); g.SetMarkerColor(col); g.Draw('P'); keep.append(g)
        f = ROOT.TF1('f%s%d' % (end, L), '[0]+[1]*x', 60, 500); f.SetParameters(o[end]['c'], o[end]['a']); f.SetLineColor(col); f.SetLineWidth(2); f.SetLineStyle(1 if end == 'S' else 2); f.Draw('same'); keep.append(f)
    lg = legend(0.21, 0.62, 0.7, 0.88, size=0.045); lg.SetHeader('SFIL%d cell 3, %d pairs' % (L + 1, len(P['p%d_2' % L])))
    lg.AddEntry(keep[-4], 'South: a = %.1f #pm %.1f ADC/pe' % (o['S']['a'], o['S']['a_err']), 'p'); lg.AddEntry(keep[-2], 'North: a = %.1f #pm %.1f ADC/pe' % (o['N']['a'], o['N']['a_err']), 'p'); lg.AddEntry(keep[-1], 'N/S amplitude ratio %.2f' % r, ''); lg.Draw(); keep.append(lg)
    if L == 3: keep.append(label('not photostatistics: North end anomalous', x=0.21, y=0.56, size=0.04, color=ROOT.kRed + 1))
c.cd(); keep.append(label('run 3062, photostatistics from the two ends of a cell: #sigma^{2} = var(x) #minus cov(x, y) = c + a A, a = ADC per photoelectron #times excess noise', x=0.01, y=0.95, size=0.03))
c.SaveAs('photostat_ns.png'); c.SaveAs('photostat_ns.pdf'); print('wrote photostat_ns.png')
