"""SFIL1-4 South energy resolution per paddle in run 3064 from the stability histograms (hT_p: SFIL1-4 South total of the fitted
amplitudes per 600 s bin), beam period only (before the TAC period), and its stability in time.

    python3 resolution_3064.py [stability_3064.root]  -> resolution_3064.txt, resolution_3064.png

Gaussian core of the total above 600 ADC (iterated +-2 sigma), as in sim/analysis/performance.py for run 3062 (29.9 / 26.5 / 20.4 /
18.8 % for P4 / P3 / P6 / P5 with SFIL1-4 alone). The 60 ADC bins of hT are coarse but fine for a 20-30 % width.
"""
import sys, numpy as np, ROOT
sys.path.insert(0, '/sphenix/user/tmengel/eic/bic_testbeam/common'); from plotstyle import setup, finish, legend, headroom
setup()
F = sys.argv[1] if len(sys.argv) > 1 else '/sphenix/tg/tg01/jets/tmengel/eic/kcu0/stability/stability_3064.root'
f = ROOT.TFile.Open(F); out = open('resolution_3064.txt', 'w')
def say(*a): s = ' '.join(str(x) for x in a); print(s); out.write(s + '\n')
TAC = 9.08 * 3600.
def core(y, x):
    m = x > 600; mu = (y[m] * x[m]).sum() / y[m].sum(); sig = np.sqrt((y[m] * (x[m] - mu) ** 2).sum() / y[m].sum())
    for it in range(6):
        w = m & (x > mu - 2 * sig) & (x < mu + 2 * sig); mu = (y[w] * x[w]).sum() / y[w].sum(); sig = np.sqrt((y[w] * (x[w] - mu) ** 2).sum() / y[w].sum())
    return mu, sig, y[m].sum()
keep = []; cv = ROOT.TCanvas('c', '', 1500, 650); cv.Divide(2, 1, 0.002, 0.002); keep.append(cv)
COL = {3: ROOT.kRed + 1, 4: ROOT.kBlue + 1, 5: ROOT.kGreen + 2, 6: ROOT.kMagenta + 1}
say('run 3064, beam period (< %.2f h): SFIL1-4 South total of the fitted amplitudes, Gaussian core above 600 ADC' % (TAC / 3600))
say('%-6s %10s %9s %9s %9s' % ('paddle', 'showers', 'mu [ADC]', 'sigma', 'sigma/mu'))
cv.cd(1); lg = legend(0.55, 0.66, 0.92, 0.87); keep.append(lg); first = True; hs = {}
for p in (4, 3, 6, 5):
    h = f.Get('hT_%d' % p); nx, ny = h.GetNbinsX(), h.GetNbinsY(); x = np.array([h.GetYaxis().GetBinCenter(j + 1) for j in range(ny)])
    tb = np.array([h.GetXaxis().GetBinCenter(i + 1) for i in range(nx)]); Y = np.array([[h.GetBinContent(i + 1, j + 1) for j in range(ny)] for i in range(nx)])
    y = Y[tb < TAC].sum(axis=0); mu, sig, n = core(y, x); say('P%d %10.0f %9.0f %9.0f %9.3f' % (p, n, mu, sig, sig / mu))
    hp = ROOT.TH1D('hp%d' % p, 'SFIL1-4 South total per paddle, beam period;#Sigma A_{fit} SFIL1-4 South [ADC];showers (normalised)', ny, x[0] - 30, x[-1] + 30); keep.append(hp)
    for j in range(ny): hp.SetBinContent(j + 1, y[j])
    if hp.Integral() > 0: hp.Scale(1. / hp.Integral())
    hp.SetLineColor(COL[p]); hp.SetLineWidth(2); hs[p] = hp; lg.AddEntry(hp, 'P%d: #mu %.0f, #sigma/#mu %.1f %%' % (p, mu, 100 * sig / mu), 'l')
    # per time bin
    res = []
    for i in range(nx):
        if tb[i] >= TAC or Y[i].sum() < 2000: res.append(np.nan); continue
        m_, s_, n_ = core(Y[i], x); res.append(s_ / m_)
    res = np.array(res); ok = np.isfinite(res); say('   per 600 s bin: sigma/mu mean %.3f, rms %.3f over %d bins' % (res[ok].mean(), res[ok].std(), ok.sum()))
    g = ROOT.TGraph(int(ok.sum()), np.ascontiguousarray(tb[ok] / 3600., float), np.ascontiguousarray(res[ok], float)); g.SetMarkerColor(COL[p]); g.SetLineColor(COL[p]); g.SetMarkerStyle(20); g.SetMarkerSize(0.8); keep.append(g); hs['g%d' % p] = g
mx = max(h.GetMaximum() for k, h in hs.items() if isinstance(k, int)); first = True
for p in (4, 3, 6, 5):
    if first: headroom(hs[p], 0, mx, 0.35); hs[p].Draw('hist'); first = False
    else: hs[p].Draw('hist same')
lg.Draw()
cv.cd(2); fr = ROOT.TH2D('fr', 'core #sigma/#mu of the SFIL1-4 South total per 600 s;time in run [h];#sigma/#mu', 10, 0, 9.3, 10, 0.1, 0.4); fr.SetStats(0); keep.append(fr); fr.Draw()
lg2 = legend(0.60, 0.70, 0.92, 0.87); keep.append(lg2)
for p in (4, 3, 6, 5): hs['g%d' % p].Draw('P'); lg2.AddEntry(hs['g%d' % p], 'P%d' % p, 'p')
lg2.Draw(); finish(cv); cv.SaveAs('resolution_3064.png'); out.close()
