"""Good-TOA selection for SFIL2 cell 3, ADC-free: cuts use only the TOA words themselves, never the ADC time.
  quiet    no TOA word on the channel before the pulse (samples 1 .. k-1), at both ends
  single   a word on only one of the two pulse samples k, k+1, at both ends
  pure     the word's code is in a band where pulse words outnumber noise words (noise shape from no-pulse triggers,
           noise level = the per-sample word rate before the pulse), at both ends; purity map built on even triggers
Cuts are applied cumulatively; at each stage the amplitude/baseline regression is refitted on even triggers of that stage and
the TOA and the ADC crossing are evaluated on the same odd-trigger showers.
    python3 toa_good.py -> toa_good.txt, toa_good.npz"""
import numpy as np
out = open('toa_good.txt', 'w')
def say(*a): print(*a, flush=True); print(*a, file=out)
z = np.load('toa_clean.npz'); pad = z['pad']; even = z['even']; NT = len(pad)
PB = {4: 2, 3: 1, 6: 8, 5: 4}; PX = {4: -80., 3: -30., 6: 30., 5: 80.}; PADS = (4, 3, 6, 5)
padok = np.isin(pad, list(PB.values())); has = z['hasS'] & z['hasN'] & padok
def gfit(v, lo=-25., hi=25.):
    v = v[(v > lo) & (v < hi)]; m, s, fsig = np.median(v), 2.5, 0.8
    for _ in range(200):
        g = fsig * np.exp(-0.5 * ((v - m) / s) ** 2) / (s * np.sqrt(2 * np.pi)); b = (1 - fsig) / (hi - lo); r = g / (g + b)
        fsig = r.mean(); m = (r * v).sum() / r.sum(); s = np.sqrt((r * (v - m) ** 2).sum() / r.sum())
    return m, s, fsig
# code purity per 16-code bin and end
PUR = {}
for e in 'SN':
    c0 = int(z['c0' + e]); hp = z['has' + e] & even; uc = (z['code' + e] - c0) % 1024
    hk = np.bincount(uc[hp] // 16, minlength=64).astype(float)
    fb = np.bincount(((z['bg' + e] - c0) % 1024) // 16, minlength=64).astype(float); fb /= fb.sum()
    rate = z['npre' + e][hp].sum() / np.maximum(z['k' + e][hp] - 1, 1).sum()          # noise words per sample before the pulse
    pur = np.clip(1 - rate * hp.sum() * fb / np.maximum(hk, 1), 0, 1); PUR[e] = pur
    say('%s: noise-word rate before the pulse %.3f per sample; code bins with purity > 0.5 hold %.2f of the pulse words' % (e, rate, hk[pur > 0.5].sum() / hk.sum()))
quiet = (z['npreS'] == 0) & (z['npreN'] == 0)
single = ~z['bothS'] & ~z['bothN']
pure = (PUR['S'][((z['codeS'] - int(z['c0S'])) % 1024) // 16] > 0.5) & (PUR['N'][((z['codeN'] - int(z['c0N'])) % 1024) // 16] > 0.5)
STAGES = (('all words', has), ('+ quiet before the pulse', has & quiet), ('+ one word only', has & quiet & single), ('+ pure code', has & quiet & single & pure))
X = np.stack([np.log(np.maximum(z['AS'], 1)), np.log(np.maximum(z['AN'], 1)), z['baseS'], z['baseN']], axis=1)
dT0 = z['tN'] - z['tS']; dA = z['tcN'] - z['tcS']; n0 = (has & ~even).sum(); rows = []
say('\nodd triggers, one paddle, clean pulse and a TOA word at both ends: %d showers' % n0)
for lab, m in STAGES:
    # amplitude/baseline regression refitted on the even triggers of this stage
    pm = {p: np.median(dT0[m & even & (pad == PB[p])]) for p in PADS}; res = dT0 - np.select([pad == PB[p] for p in PADS], [pm[p] for p in PADS], np.nan)
    fm = m & even & (np.abs(res) < 10); A_ = np.column_stack([np.ones(fm.sum()), X[fm]]); beta = np.linalg.lstsq(A_, res[fm], rcond=None)[0]
    dT = dT0 - X @ beta[1:]; mo = m & ~even
    r = {}
    for key, d in (('TOA', dT), ('ADC', dA)):
        fits = [gfit(d[mo & (pad == PB[p])]) for p in PADS]; mu = np.array([f[0] for f in fits]); sg = np.array([f[1] for f in fits]); fr = np.array([f[2] for f in fits])
        c1 = np.polyfit([PX[p] for p in PADS], mu, 1)[0]; r[key] = (sg.min(), sg.max(), fr.mean(), 2 / c1, (sg / c1).min(), (sg / c1).max(), (mu[3] - mu[0]) / sg.mean())
    dd = dT0[mo] - np.median(dT0[mo]); diag = (np.abs(dd) < 5).mean(); corr = np.corrcoef(z['tN'][mo], z['tS'][mo])[0, 1]
    say('%-26s n %6d (%.2f)  corr %.3f  within 5 ns %.2f | TOA sigma %.2f-%.2f ns core %.2f  %3.0f mm/ns  sigma_x %3.0f-%3.0f mm  sep %.2f | ADC same showers: sigma %.2f-%.2f ns  sigma_x %3.0f-%3.0f mm  sep %.2f' %
        (lab, mo.sum(), mo.sum() / n0, corr, diag, r['TOA'][0], r['TOA'][1], r['TOA'][2], r['TOA'][3], r['TOA'][4], r['TOA'][5], r['TOA'][6], r['ADC'][0], r['ADC'][1], r['ADC'][4], r['ADC'][5], r['ADC'][6]))
    rows.append((mo.sum() / n0, corr, diag) + r['TOA'] + r['ADC'])
np.savez('toa_good.npz', dTgood=dT, mgood=m, rows=np.array(rows), labels=np.array([s[0] for s in STAGES]), purS=PUR['S'], purN=PUR['N'], quiet=quiet, single=single, pure=pure)
say('wrote toa_good.npz')
