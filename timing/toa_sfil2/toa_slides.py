"""Single-purpose figures for the TOA calibration slides (SFIL2 cell 3 on KCU0, plus the KCU0 channel survey).
Numbers that are already in the logs (toa_timing.log, toa_explore.txt, toa_baseline.txt) are parsed or copied; the code maps
are refilled from kcu0_3062.root with the same selection as toa_cdc2.py.
    python3 toa_slides.py [nmax] -> slides/*.png"""
import numpy as np, uproot, sys, os, re, ROOT
sys.path.insert(0, '/sphenix/user/tmengel/eic/bic_testbeam/common'); from plotstyle import setup, legend, colz
setup(); ROOT.gROOT.SetBatch(True); os.makedirs('slides', exist_ok=True)
NMAX = int(sys.argv[1]) if len(sys.argv) > 1 else 600000
RED, BLUE, BLK, GRN, ORG = ROOT.kRed + 1, ROOT.kBlue + 1, ROOT.kBlack, ROOT.kGreen + 2, ROOT.kOrange + 7
keep = []
def graph(x, y, mk, col, size=1.8, ls=1):
    g = ROOT.TGraph(len(x), np.ascontiguousarray(x, float), np.ascontiguousarray(y, float)); g.SetMarkerStyle(mk); g.SetMarkerSize(size); g.SetMarkerColor(col); g.SetLineColor(col); g.SetLineWidth(2); g.SetLineStyle(ls); keep.append(g); return g
def frame(x0, y0, x1, y1, axes, nd=None):
    fr = ROOT.gPad.DrawFrame(x0, y0, x1, y1, ';' + axes); keep.append(fr)
    if nd: fr.GetXaxis().SetNdivisions(nd)
    return fr

# 1. KCU0 survey: noise words on pedestal samples vs efficiency on pulses, 32 SFIL channels
rows = [re.match(r'\s+([SN])(\d)c(\d) ch(\d+)\s+noise\s+([\d.]+)\s+efficiency\s+([\d.]+)', l) for l in open('../../sync/out/toa_timing.log')]
rows = [(m.group(1), int(m.group(2)), int(m.group(3)), int(m.group(4)), float(m.group(5)), float(m.group(6))) for m in rows if m]
c = ROOT.TCanvas('c1', '', 1600, 1000); frame(-0.02, -0.05, 0.55, 1.45, 'TOA words per pedestal sample;fraction of pulses with a TOA word', 505)
lg = legend(0.2, 0.68, 0.9, 0.88, size=0.042)
for end, mk, col, lab in (('S', 20, BLUE, 'South ends (ASIC0)'), ('N', 21, RED, 'North ends (ASIC2)')):
    r = [x for x in rows if x[0] == end]; g = graph([x[4] for x in r], [x[5] for x in r], mk, col, ls=0); g.Draw('p same'); lg.AddEntry(g, lab, 'p')
s2 = [x for x in rows if x[1] == 2 and x[2] == 3]; g = graph([x[4] for x in s2], [x[5] for x in s2], 24, BLK, size=3.2); g.SetMarkerSize(3.2); g.Draw('p same'); lg.AddEntry(g, 'SFIL2 cell 3 (ch 29 / 173): this study', 'p')
lg.Draw(); keep.append(lg); c.SaveAs('slides/survey.png')

# 2. which sample carries the pulse word (toa_explore.txt)
txt = open('toa_explore.txt').read()
def wordpos(name):
    blk = txt.split(name + ':')[1].split('words at k+1')[0]
    return [float(v) for v in re.findall(r'crossing sample [+-]\d: ([\d.]+)', blk)], float(re.search(r'samples 1-3: ([\d.]+)', blk).group(1))
c = ROOT.TCanvas('c2', '', 1600, 1000); frame(-2.6, 0, 3.6, 1.0, 'TOA sample - ADC crossing sample;fraction of pulses with a word there', 506)
lg = legend(0.45, 0.66, 0.92, 0.88, size=0.042)
for name, mk, col, lab in (('SFIL2 S c3', 20, BLUE, 'South'), ('SFIL2 N c3', 21, RED, 'North')):
    v, rate = wordpos(name); g = graph(np.arange(-2, 4), v, mk, col); g.Draw('pl same'); lg.AddEntry(g, '%s (pedestal rate %.2f)' % (lab, rate), 'pl')
    l = ROOT.TLine(-2.6, rate, 3.6, rate); l.SetLineColor(col); l.SetLineStyle(2); l.SetLineWidth(2); l.Draw(); keep.append(l)
lk = ROOT.TLine(); lk.SetLineStyle(2); lk.SetLineWidth(2); keep.append(lk); lg.AddEntry(lk, 'pedestal-sample word rate (dashed)', 'l'); lg.Draw(); keep.append(lg); c.SaveAs('slides/wordpos.png')

# 3./4. code vs ADC-crossing phase, raw and unwrapped with the mode table (same selection as toa_cdc2.py)
f = uproot.open('/sphenix/user/tmengel/eic/bic_testbeam/data/kcu0_3062.root')
def load(br):
    a = f['sfil2'].arrays([br + '_adc', br + '_toa'], entry_stop=NMAX, library='np'); W = a[br + '_adc'].astype(float); T = a[br + '_toa'].astype(int)
    n = len(W); good = W[:, 0] >= 0; y = W - W[:, 1:4].mean(axis=1)[:, None]
    ipk = y[:, 4:16].argmax(axis=1) + 4; A = y[np.arange(n), ipk]; half = 0.5 * A; k = ipk.copy()
    for d in range(1, 6):
        kk = np.clip(ipk - d, 0, 19); k = np.where(y[np.arange(n), kk] >= half, kk, k)
    k = np.clip(k, 1, 18); y0 = y[np.arange(n), k - 1]; y1 = y[np.arange(n), k]; tc = 25. * (k - 1 + np.clip((half - y0) / np.maximum(y1 - y0, 1e-3), 0, 1))
    sel = good & (A > 150) & (A < 950) & (ipk >= 5) & (ipk <= 14)
    codes, ph = [], []
    for d in (0, 1):
        w = T[np.arange(n), k + d]; m = sel & (w > 0); codes.append(w[m]); ph.append(tc[m] - 25. * (k + d)[m])
    return np.concatenate(codes), np.concatenate(ph)
def mode_table(codes, phase, c0):
    cb = ((codes - c0) % 1024) // 16; xs, ys = [], []
    for b in range(64):
        m = cb == b
        if m.sum() < 200: continue
        h, e = np.histogram(phase[m], bins=np.arange(-40, 40.01, 0.5)); h = np.convolve(h, np.ones(5) / 5, 'same'); xs.append(b * 16 + 8); ys.append(e[h.argmax()] + 0.25)
    return np.array(xs), np.array(ys)
C0 = {'Sc3': 0, 'Nc3': 528}   # unwrap offsets found by toa_cdc2.py / toa_baseline.py
D = {br: load(br) for br in C0}
for tag, unwrap in (('raw', False), ('unwrapped', True)):
    c = ROOT.TCanvas('c' + tag, '', 2400, 900); c.Divide(2, 1)
    for i, (br, lab) in enumerate((('Sc3', 'South'), ('Nc3', 'North'))):
        c.cd(i + 1); colz(right=0.15); ROOT.gPad.SetLogz(); codes, ph = D[br]; cc = (codes - C0[br]) % 1024 if unwrap else codes
        h = ROOT.TH2F('h%s%d' % (tag, i), ';%s TOA code;t_{cross} - 25 s_{word} [ns]' % ('unwrapped' if unwrap else 'raw'), 128, 0, 1024, 100, -45, 15); keep.append(h)
        h.FillN(len(cc), cc.astype(float), ph.astype(float), np.ones(len(cc))); h.GetXaxis().SetNdivisions(505); h.Draw('colz')
        t = ROOT.TLatex(); t.SetNDC(); t.SetTextSize(0.05); t.DrawLatex(0.2, 0.82, lab); keep.append(t)
        if unwrap:
            xs, ys = mode_table(codes, ph, C0[br]); g = graph(xs, ys, 1, RED, ls=1); g.SetLineWidth(4); g.Draw('l same')
            lg = legend(0.2, 0.66, 0.8, 0.78, size=0.042); lg.AddEntry(g, 'mode per 16 codes (table)', 'l'); lg.Draw(); keep.append(lg)
    c.SaveAs('slides/code_%s.png' % tag)

# 5. per-trigger baseline offset: TOA - crossing shift per baseline bin and amplitude bin (toa_baseline.txt, even triggers)
bt = open('toa_baseline.txt').read().split('\n'); BC = [-40, -25, -16, -9, -3, 3, 9, 16, 25, 40]; AL = ['150-250', '250-350', '350-450', '450-600', '600-950']
c = ROOT.TCanvas('c5', '', 2400, 900); c.Divide(2, 1)
for i, (end, lab) in enumerate((('S', 'South'), ('N', 'North'))):
    j = [n for n, l in enumerate(bt) if l.startswith(end + ':')][0]; tab = np.array([[float(v) for v in bt[j + 2 + r].split()] for r in range(10)])
    c.cd(i + 1); frame(-45, -6, 45, 7, 'per-trigger baseline offset [ADC];TOA - crossing shift [ns]', 505)
    lg = legend(0.2, 0.7, 0.9, 0.86, size=0.04, ncol=3)
    for a, (mk, col) in enumerate(zip((20, 21, 22, 23, 33), (BLUE, GRN, ORG, RED, BLK))):
        g = graph(BC, tab[:, a], mk, col, size=1.6); g.Draw('pl same'); lg.AddEntry(g, 'A ' + AL[a], 'pl')
    lg.Draw(); keep.append(lg); t = ROOT.TLatex(); t.SetNDC(); t.SetTextSize(0.05); t.DrawLatex(0.2, 0.2, lab); keep.append(t)
c.SaveAs('slides/baseline.png')

# 6./7. results on odd triggers (toa_baseline.txt)
R = np.array([[float(v) for v in re.findall(r'TOA ([\d.]+)\s+ADC ([\d.]+)\s+combined ([\d.]+).*TOA-only ([\d.]+), ADC-only ([\d.]+), common ([\d.]+)', l)[0]] for l in bt if l.strip().startswith('A ') and 'combined' in l])
xm = np.array([200, 300, 400, 525, 775.])
c = ROOT.TCanvas('c6', '', 2400, 900); c.Divide(2, 1)
c.cd(1); frame(100, 0, 1000, 4.2, 'min(A_{S}, A_{N}) [ADC];core #sigma(t_{N} - t_{S}) [ns]', 505); lg = legend(0.3, 0.7, 0.9, 0.88)
for col_, mk, cl, lab in ((0, 20, RED, 'TOA (baseline-corrected)'), (1, 21, BLUE, 'ADC 50 % crossing'), (2, 22, BLK, 'TOA + ADC combined')):
    g = graph(xm, R[:, col_], mk, cl); g.Draw('pl same'); lg.AddEntry(g, lab, 'pl')
lg.Draw(); keep.append(lg)
c.cd(2); PXs = [-80, -30, 30, 80]; mT = [float(v) for v in re.search(r'TOA N - S paddle means ([\d. /-]+) ns', open('toa_baseline.txt').read()).group(1).split(' / ')]
mA = [float(v) for v in re.search(r'ADC N - S paddle means ([\d. /-]+) ns', open('toa_baseline.txt').read()).group(1).split(' / ')]
frame(-110, -1.2, 110, 3.2, 'paddle position x [mm];mean t_{N} - t_{S} [ns]', 505); lg = legend(0.2, 0.7, 0.9, 0.88)
for m_, mk, cl, lab in ((mT, 20, RED, 'TOA: v_{eff} = 226 mm/ns'), (mA, 21, BLUE, 'ADC crossing: v_{eff} = 146 mm/ns')):
    g = graph(PXs, m_, mk, cl); g.Draw('p same'); p = np.polyfit(PXs, m_, 1); fl = ROOT.TF1('f%d' % mk, 'pol1', -100, 100); fl.SetParameters(p[1], p[0]); fl.SetLineColor(cl); fl.SetLineStyle(2); fl.Draw('same'); keep.append(fl); lg.AddEntry(g, lab, 'p')
lg.Draw(); keep.append(lg); c.SaveAs('slides/result.png')
c = ROOT.TCanvas('c7', '', 1600, 1000); frame(100, 0, 1000, 3.0, 'min(A_{S}, A_{N}) [ADC];contribution to #sigma(t_{N} - t_{S}) [ns]', 505); lg = legend(0.3, 0.7, 0.9, 0.88)
for col_, mk, cl, lab in ((3, 20, RED, 'TOA only'), (4, 21, BLUE, 'ADC only'), (5, 22, BLK, 'common to both (light timing)')):
    g = graph(xm, R[:, col_], mk, cl); g.Draw('pl same'); lg.AddEntry(g, lab, 'pl')
lg.Draw(); keep.append(lg); c.SaveAs('slides/split.png')
print('wrote slides/*.png')
