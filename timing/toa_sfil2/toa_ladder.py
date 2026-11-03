"""SFIL2 cell 3 TOA: the calibration built up one step at a time, every step measured on the same showers (odd triggers, both
ends good in the final calibration) with the same measure (core sigma of N - S after removing each paddle's mean).
Steps: ADC sample only; pulse word + unwrap + nominal 24.4 ps/code; + code-to-time table; + pulse-height walk; + baseline.
    python3 toa_ladder.py -> toa_ladder.txt, slides/ladder.png"""
import numpy as np, uproot, sys, ROOT
sys.path.insert(0, '/sphenix/user/tmengel/eic/bic_testbeam/common'); from plotstyle import setup, legend
setup(); ROOT.gROOT.SetBatch(True); out = open('toa_ladder.txt', 'w')
def say(*a): print(*a, flush=True); print(*a, file=out)
z = np.load('toa_baseline.npz'); NT = len(z['pad']); pad = z['pad']; even = z['even']
f = uproot.open('/sphenix/user/tmengel/eic/bic_testbeam/data/kcu0_3062.root')
PB = {4: 2, 3: 1, 6: 8, 5: 4}; C0 = {'S': 0, 'N': 528}
def core(v, w=2.5):
    v = v[np.isfinite(v)]
    if len(v) < 30: return np.nan, np.nan
    m = np.median(v); s = 1.4826 * np.median(np.abs(v - m))
    if s == 0: return v.mean(), v.std()   # discrete 25 ns steps: no core, quote the plain std
    for _ in range(6): g = np.abs(v - m) < w * s; m, s = v[g].mean(), v[g].std()
    return m, s / 0.9
def mode_table(codes, phase):
    cb = codes // 16; xs, ys = [], []
    for b in range(64):
        m = cb == b
        if m.sum() < 200: continue
        h, e = np.histogram(phase[m], bins=np.arange(-40, 40.01, 0.5)); h = np.convolve(h, np.ones(5) / 5, 'same'); xs.append(b * 16 + 8); ys.append(e[h.argmax()] + 0.25)
    return np.interp(np.arange(1024), xs, ys)
def walk(t, tc, A, ok):
    la = np.log(np.maximum(A, 1)); bins = np.linspace(np.log(150), np.log(950), 11); xb, yb = [], []
    for lo, hi in zip(bins[:-1], bins[1:]):
        m = ok & even & (la >= lo) & (la < hi)
        if m.sum() > 300: xb.append(np.median(la[m])); yb.append(core((t - tc)[m])[0])
    return t - np.interp(la, xb, yb)
ST = {}
for end, br in (('S', 'Sc3'), ('N', 'Nc3')):
    a = f['sfil2'].arrays([br + '_adc', br + '_toa'], entry_stop=NT, library='np'); W = a[br + '_adc'].astype(float); T = a[br + '_toa'].astype(int)
    n = len(W); y = W - W[:, 1:4].mean(axis=1)[:, None]
    ipk = y[:, 4:16].argmax(axis=1) + 4; A = y[np.arange(n), ipk]; half = 0.5 * A; k = ipk.copy()
    for d in range(1, 6):
        kk = np.clip(ipk - d, 0, 19); k = np.where(y[np.arange(n), kk] >= half, kk, k)
    k = np.clip(k, 1, 18); tc = z['tc' + end]; ok = z['ok' + end]
    cand = [(T[np.arange(n), k], k), (T[np.arange(n), k + 1], k + 1)]
    def pick(fn):
        t = np.full(n, np.nan); rb = np.full(n, np.inf)
        for w, s in cand:
            tt = np.where(w > 0, fn(w, s), np.nan); rr = np.abs(tt - tc); b = np.isfinite(tt) & (rr < rb); t[b] = tt[b]; rb[b] = rr[b]
        return t
    st = {}
    st['tick'] = 25. * k                                                          # ADC sample only
    lin = pick(lambda w, s: 25. * s + ((w - C0[end]) % 1024) * 25. / 1024)        # word + unwrap, nominal 24.4 ps/code
    st['lin'] = lin - core((lin - tc)[ok & even])[0]
    codes = np.concatenate([((w - C0[end]) % 1024)[ok & even & (w > 0)] for w, s in cand]); ph = np.concatenate([(tc - 25. * s)[ok & even & (w > 0)] for w, s in cand])
    tab = mode_table(codes, ph); st['table'] = pick(lambda w, s: 25. * s + tab[(w - C0[end]) % 1024])
    st['walk'] = walk(st['table'], tc, A, ok)
    st['base'] = z['t' + end]
    for s_ in st: say('%s %-6s per-end core sigma vs ADC crossing (odd): %.2f ns' % (end, s_, core((st[s_] - tc)[ok & ~even])[1]))
    ST[end] = st
both = z['okS'] & z['okN'] & ~even; Amin = np.minimum(z['AS'], z['AN'])
def nsig(d):
    r = np.full(len(d), np.nan)
    for p in PB.values():
        m = both & (pad == p); r[m] = d[m] - core(d[m])[0]
    pooled = core(r[both])[1]; rng = []
    for lo, hi in ((150, 250), (250, 350), (350, 450), (450, 600), (600, 950)):
        m = both & (Amin >= lo) & (Amin < hi); rng.append(core(r[m])[1])
    return pooled, min(rng), max(rng), np.sqrt(np.nanmean(r[both] ** 2))
LAB = [('tick', 'ADC sample only (25 ns steps)'), ('lin', 'pulse word, unwrapped, 24.4 ps/code'), ('table', '+ code-to-time table'),
       ('walk', '+ pulse-height correction'), ('base', '+ baseline correction')]
say('\nsame %d showers (odd triggers, both ends good); N - S after paddle means removed: pooled core sigma, range over amplitude bins, rms' % both.sum())
res = []
for key, lab in LAB:
    p, lo, hi, rms = nsig(ST['N'][key] - ST['S'][key]); res.append(p); say('  %-40s %.2f ns  (%.2f-%.2f)  rms %.2f' % (lab, p, lo, hi, rms))
p, lo, hi, rms = nsig(z['tcN'] - z['tcS']); say('  %-40s %.2f ns  (%.2f-%.2f)  rms %.2f' % ('ADC 50 % crossing (reference)', p, lo, hi, rms)); adc = p
c = ROOT.TCanvas('c', '', 1600, 900); c.SetBottomMargin(0.30); c.SetLeftMargin(0.12)
h = ROOT.TH1F('h', ';;core #sigma(t_{N} - t_{S}) [ns]', len(LAB), 0, len(LAB)); h.SetMaximum(max(res) * 1.25); h.SetMinimum(0)
for i, (key, lab) in enumerate(LAB): h.SetBinContent(i + 1, res[i]); h.GetXaxis().SetBinLabel(i + 1, ['ADC sample only', 'word + unwrap', '+ code table', '+ pulse height', '+ baseline'][i])
h.GetXaxis().SetLabelSize(0.05); h.SetFillColor(ROOT.kRed - 9); h.SetLineColor(ROOT.kRed + 1); h.SetBarWidth(0.6); h.SetBarOffset(0.2); h.Draw('bar')
ln = ROOT.TLine(0, adc, len(LAB), adc); ln.SetLineColor(ROOT.kBlue + 1); ln.SetLineWidth(3); ln.SetLineStyle(2); ln.Draw()
lg = legend(0.55, 0.80, 0.92, 0.90); lg.AddEntry(h, 'TOA', 'f'); lg.AddEntry(ln, 'ADC 50 %% crossing (%.1f ns)' % adc, 'l'); lg.Draw()
c.SaveAs('slides/ladder.png'); say('wrote slides/ladder.png')
