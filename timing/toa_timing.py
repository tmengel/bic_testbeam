"""TOA / TOT timing study, run 3062: can the H2GCROC3 TDC give ps timing between the North and South
ends of the SFIL fibres (both on KCU0)?

The 10-bit TOA is the TDC word (datasheet 1.2.2): fine = bits 2:0 (FTDC, 25 ps LSB), mid = bits 7:3
(CTDC, 0.2 ns LSB), coarse = bits 9:8 (160 MHz Gray counter, 6.25 ns). In these data the TOA word is
stored in the sample AFTER the one in which the ADC crosses the pulse edge, and the mapping of the
fields to time is not the naive one (coarse is neither plain nor Gray, the mid slope depends on the
coarse value), so each channel is calibrated empirically here against the ADC 20%-crossing time:
t_TOA = 25 x sample + intercept(coarse) + slope(coarse) x (mid + fine/8). That reference is itself
non-linear at the ns level (linear interpolation of a 2-sample rising edge), so the TOA resolution is
measured TOA-vs-TOA between two calibrated channels of the same waveform (KCU1 single SiPMs).

Inputs:  output/kcutime_0001.root  (ts_dump kcutime:  KCU0, the 36 kcuadc channels, 400k waveforms)
         output/kcu1time_0000.root (ts_dump kcu1time: KCU1 ch 144-179, SFIL5 South single SiPMs)
Outputs: out/toa_timing.png, out/toa_timing.log

    python3 toa_timing.py
"""
import numpy as np, uproot, ROOT
from scipy.stats import spearmanr
import sys; sys.path.insert(0, '/sphenix/user/tmengel/eic/bic_testbeam/common')
from plotstyle import setup, finish, legend, headroom, colz

LOG = open('out/toa_timing.log', 'w')
def say(*a):
    s = ' '.join(str(x) for x in a); print(s); LOG.write(s + '\n')
def mad(x):
    x = x[np.isfinite(x)]
    return 1.4826 * np.median(np.abs(x - np.median(x))) if len(x) else np.nan

kAmpCh = [2, 21, 63, 43, 6, 25, 59, 37, 10, 29, 52, 0, 14, 33, 48, 4, 146, 165, 207, 187, 150, 169, 203, 186, 154, 173, 196, 144, 158, 177, 192, 148, 84, 86, 81, 97]
def lab0(k):
    if k < 16: return f'S{k%4+1}c{k//4+1} ch{kAmpCh[k]}'
    if k < 32: return f'N{k%4+1}c{(k-16)//4+1} ch{kAmpCh[k]}'
    return f'P{[3,4,5,6][k-32]} ch{kAmpCh[k]}'

def load(fin, tree, usecrc, nmax):
    ch = next(uproot.iterate(fin + ':' + tree, ['wf', 'toa', 'tot', 'ncrc'], library='np', step_size=nmax))
    ok = (ch['ncrc'] == 0) if usecrc else np.ones(len(ch['ncrc']), bool)
    wf = ch['wf'][ok].astype(np.float32); toa = ch['toa'][ok].astype(np.int32); tot = ch['tot'][ok].astype(np.int32)
    ped = wf[:, :, 1:4].mean(2); v = wf - ped[:, :, None]; A = v[:, :, 4:13].max(2)
    return dict(toa=toa, tot=tot, v=v, A=A, ped=ped)

def per_channel(D, k, frac=0.2):
    """ADC crossing time at frac x amplitude (linear interpolation), and the first TOA word in the crossing sample or the two after it"""
    v, toa, A = D['v'], D['toa'], D['A']; n = len(A)
    seg = v[:, k, 4:13]; lvl = frac * A[:, k]; above = seg >= lvl[:, None]; s = above.argmax(1) + 4; has = above.any(1) & (s >= 5)
    vs = v[np.arange(n), k, s]; vp = v[np.arange(n), k, s - 1]; t = 25.0 * ((s - 1) + (lvl - vp) / np.maximum(vs - vp, 1e-3))
    w = np.stack([toa[np.arange(n), k, np.clip(s + d, 0, 19)] for d in (0, 1, 2)], 1); nz = w != 0; have = nz.any(1); d = nz.argmax(1)
    return dict(ok=has & have & (A[:, k] > 60), t=t, s=s, st=s + d, code=w[np.arange(n), d], A=A[:, k], toa=toa[:, k, :], tot=D['tot'][:, k, :])

def calibrate(ch, amin=250, amax=900):
    """intercept and slope vs mid, per coarse value, from the ADC crossing time"""
    m = ch['ok'] & (ch['A'] > amin) & (ch['A'] < amax); base = (ch['t'] - 25.0 * ch['st'])[m]; c = ch['code'][m]
    coarse = c >> 8; mid = (c >> 3) & 31
    cal = {}
    for q in range(4):
        mq = coarse == q
        if mq.sum() >= 150: cal[q] = np.polyfit(mid[mq], base[mq], 1)
    return cal

def model_time(ch, cal):
    c = ch['code']; coarse = c >> 8; mid = (c >> 3) & 31; fine = c & 7; tt = np.full(len(c), np.nan)
    for q, p in cal.items():
        mq = coarse == q; tt[mq] = 25.0 * ch['st'][mq] + p[1] + p[0] * (mid[mq] + fine[mq] / 8.0)
    return tt

# ======================================================================================= KCU0
K0 = load('../../output/kcutime_0001.root', 'kcutime', True, 100000)
say(f'KCU0: {len(K0["A"])} waveforms (ncrc == 0)')
say('\n1. TOA thresholds per channel: noise TOA words per pedestal sample (samples 0-5 of waveforms with A < 30),')
say('   and the fraction of hits (A > 100) with a TOA word in samples 6-8')
noise0, eff0, nhit0 = np.zeros(36), np.zeros(36), np.zeros(36, int)
for k in range(36):
    nohit = K0['A'][:, k] < 30; hit = K0['A'][:, k] > 100
    noise0[k] = (K0['toa'][nohit, k, 0:6] != 0).mean() if nohit.sum() else np.nan
    eff0[k] = (K0['toa'][hit, k, 6:9] != 0).any(1).mean() if hit.sum() else np.nan
    nhit0[k] = hit.sum()
cls = lambda k: 'threshold in the noise' if noise0[k] > 0.05 else ('no TOA on pulses' if eff0[k] < 0.3 else 'OK')
for k in range(36):
    say(f'   {lab0(k):14s}  noise {noise0[k]:6.3f}  efficiency {eff0[k]:5.2f}  (hits {nhit0[k]:6d})  {cls(k)}')
n_noise = sum(1 for k in range(32) if cls(k) == 'threshold in the noise'); n_dead = sum(1 for k in range(32) if cls(k) == 'no TOA on pulses'); n_ok = 32 - n_noise - n_dead
say(f'   SFIL channels: {n_noise} with the threshold in the noise, {n_dead} that never fire on pulses, {n_ok} usable ({[lab0(k) for k in range(32) if cls(k) == "OK"]})')

# ======================================================================================= KCU1
K1 = load('../../output/kcu1time_0000.root', 'kcu1time', False, 100000)
say(f'\nKCU1: {len(K1["A"])} waveforms (ASIC2 lines are CRC-clean; the dead ASIC0/1 board is not used)')
good1 = []
for k in range(36):
    if ((K1['ped'][:, k] > 1000) | (K1['ped'][:, k] < 5)).mean() > 0.5: continue
    nohit = K1['A'][:, k] < 30; hit = K1['A'][:, k] > 100
    if hit.sum() < 300: continue
    noise = (K1['toa'][nohit, k, 0:6] != 0).mean(); eff = (K1['toa'][hit, k, 5:9] != 0).any(1).mean()
    if noise < 0.02 and eff > 0.9: good1.append(k)
say(f'   channels with a sane TOA threshold (noise < 2%, efficiency > 90%): {len(good1)}: {[144 + k for k in good1]}')

# ---- 2. where the TOA word is, and which bits carry time (KCU1 pooled, hits at A 150-400)
say('\n2. KCU1: position of the TOA word and bit-by-bit correlation with the ADC 20% crossing time (hits with A 150-400)')
T, C, DD = [], [], []
for k in good1:
    ch = per_channel(K1, k); m = ch['ok'] & (ch['A'] > 150) & (ch['A'] < 400)
    T.append(ch['t'][m] - 25.0 * ch['st'][m]); C.append(ch['code'][m]); DD.append((ch['st'] - ch['s'])[m])
t = np.concatenate(T); c = np.concatenate(C); dd = np.concatenate(DD)
say('   TOA sample - ADC crossing sample: ' + ', '.join(f'{d:+d}: {100*(dd == d).mean():.1f}%' for d in (0, 1, 2)))
bitrho = [spearmanr(t, (c >> b) & 1)[0] for b in range(9, -1, -1)]
bitdt = [t[((c >> b) & 1) == 1].mean() - t[((c >> b) & 1) == 0].mean() for b in range(9, -1, -1)]
say('   bit:                 ' + ' '.join(f'{b:6d}' for b in range(9, -1, -1)))
say('   rank correlation:    ' + ' '.join(f'{r:+6.2f}' for r in bitrho))
say('   <t | bit=1> - <t | bit=0> [ns]: ' + ' '.join(f'{d:+6.2f}' for d in bitdt))
say('   (a working 0.2 ns CTDC would put ~3 ns on bit 7 and ~1.6 ns on bit 6; the coarse bits carry +5 / -2 ns)')

# ---- 3. per-channel calibration against the ADC and its verification
say('\n3. per-channel calibration t_TOA = 25 x sample + intercept(coarse) + slope(coarse) x (mid + fine/8), hits A 250-900')
cal1 = {}; ch1 = {}
for k in good1:
    ch = per_channel(K1, k); cal = calibrate(ch); cal1[k] = cal; ch1[k] = ch
for k in good1[:3] + [good1[-1]]:
    cal = cal1[k]; ch = ch1[k]; tt = model_time(ch, cal); m = ch['ok'] & np.isfinite(tt) & (ch['A'] > 250) & (ch['A'] < 900)
    say(f'   KCU1 ch{144 + k}: ' + '; '.join(f'coarse {q}: {p[1]:6.2f} {p[0]:+.3f} x mid' for q, p in sorted(cal.items())) + f'  | ADC - model MAD {mad(ch["t"][m] - tt[m]):.2f} ns (n={m.sum()})')
say('   The slope vs mid changes with the coarse value and the intercepts are not multiples of 6.25 ns: the fields carry')
say('   sub-ns information, but the mapping is not the datasheet one and the ADC reference itself is non-linear at the ns level.')
c196 = per_channel(K0, 26); cal196 = calibrate(c196); tt196 = model_time(c196, cal196)
for lo, hi in ((250, 500), (500, 900), (900, 2000)):
    m = c196['ok'] & np.isfinite(tt196) & (c196['A'] >= lo) & (c196['A'] < hi)
    say(f'   KCU0 N3c3 ch196 (the one usable SFIL channel), A {lo}-{hi}: ADC - model MAD {mad(c196["t"][m] - tt196[m]):.2f} ns (n={m.sum()}) - limited by the ADC reference')

# ---- 4. TOA resolution: TOA - TOA between two calibrated KCU1 channels of the same waveform
say('\n4. TOA resolution from pairs of calibrated KCU1 channels in the same waveform (A 250-900 on both)')
PAIRS = [(k, k2) for i, k in enumerate(good1) for k2 in good1[i + 1:] if (144 + k2) - (144 + k) in (4, 8)][:6]
pair_res = []
for a, b in PAIRS:
    ta = model_time(ch1[a], cal1[a]); tb = model_time(ch1[b], cal1[b])
    m = ch1[a]['ok'] & ch1[b]['ok'] & np.isfinite(ta) & np.isfinite(tb) & (ch1[a]['A'] > 250) & (ch1[a]['A'] < 900) & (ch1[b]['A'] > 250) & (ch1[b]['A'] < 900)
    if m.sum() < 200: continue
    dt = ta[m] - tb[m]; g = 1 / np.sqrt(ch1[a]['A'][m]) - 1 / np.sqrt(ch1[b]['A'][m])
    core = np.abs(dt - np.median(dt)) < 3 * max(mad(dt), 0.1); p = np.polyfit(g[core], dt[core], 1); dtc = dt - p[0] * g
    dta = ch1[a]['t'][m] - ch1[b]['t'][m]
    pair_res.append((a, b, m.sum(), mad(dt), mad(dtc), mad(dta), dtc - np.median(dtc), dta - np.median(dta)))
    say(f'   ch{144 + a} - ch{144 + b}: n={m.sum():5d}  TOA MAD {mad(dt):.2f} ns, timewalk-corrected {mad(dtc):.2f} ns -> {mad(dtc)/np.sqrt(2)*1e3:.0f} ps per channel  |  ADC 20% crossing difference MAD {mad(dta):.2f} ns -> {mad(dta)/np.sqrt(2)*1e3:.0f} ps per channel')

# ---- 5. North - South on KCU0: what exists
say('\n5. KCU0 North - South with TOA on both ends')
c52 = per_channel(K0, 10); both = c196['ok'] & c52['ok']
say(f'   the usable channel N3c3 ch196 has {c196["ok"].sum()} TOA hits; its South partner S3c3 ch52 has {c52["ok"].sum()} (threshold above the pulses); pairs {both.sum()}')
cal52 = calibrate(c52, 100, 2000); t52 = model_time(c52, cal52); mm = both & np.isfinite(t52) & np.isfinite(tt196)
if mm.sum() > 50:
    dt = t52[mm] - tt196[mm]
    say(f'   S3c3 - N3c3 TOA - TOA: n={mm.sum()}, MAD {mad(dt):.2f} ns (S3c3 fires at the very top of its pulses only); ADC difference MAD {mad(c52["t"][mm] - c196["t"][mm]):.2f} ns')
say('   -> no North - South pair with working TOA on both ends exists in run 3062; the ADC crossing gives 1.6 - 2.0 ns per end (method.md 4.3).')

# ---- 6. TOT
say('\n6. TOT (fires only on saturated pulses, Tc set): fraction of hits with TOT, KCU0 beam row')
for k in (8, 9, 10, 11, 24, 25, 26, 27):
    hit = K0['A'][:, k] > 100; say(f'   {lab0(k):14s} hits {hit.sum():6d}  with TOT {100 * (K0["tot"][hit, k, 4:13] > 0).any(1).mean():5.1f}%  with ADC at 1023 {100 * (K0["A"][hit, k] + K0["ped"][hit, k] >= 1023).mean():5.1f}%')
say('   TOT is a usable amplitude for the saturated pulses (10-bit code, decompress 8(v-512)+4 for v >= 512), hence a timewalk')
say('   handle where the ADC is clipped; it adds nothing to timing by itself.')

# ======================================================================================= figure
setup()
cv = ROOT.TCanvas('c', '', 1800, 1100); cv.Divide(3, 2); keep = []
# (a) KCU0 thresholds
cv.cd(1)
h1 = ROOT.TH2D('h1', 'KCU0 TOA thresholds: noise rate vs efficiency, 32 SFIL channels;noise TOA words per pedestal sample;fraction of hits with a TOA word at the pulse', 10, -0.02, 0.55, 10, -0.05, 1.05); keep.append(h1); h1.Draw()
gS = ROOT.TGraph(); gN = ROOT.TGraph(); g196 = ROOT.TGraph()
for k in range(32):
    if not np.isfinite(noise0[k]) or not np.isfinite(eff0[k]): continue
    (g196 if k == 26 else (gS if k < 16 else gN)).AddPoint(noise0[k], eff0[k])
for g, col, st in ((gS, ROOT.kBlue + 1, 20), (gN, ROOT.kRed + 1, 21), (g196, ROOT.kGreen + 2, 29)):
    g.SetMarkerColor(col); g.SetMarkerStyle(st); g.SetMarkerSize(1.4 if st != 29 else 2.4); g.Draw('P same'); keep.append(g)
lg = legend(0.45, 0.2, 0.93, 0.42); lg.AddEntry(gS, 'South (ASIC0)', 'p'); lg.AddEntry(gN, 'North (ASIC2)', 'p'); lg.AddEntry(g196, 'N3c3 ch196: the one usable channel', 'p'); lg.Draw(); keep.append(lg)
# (b) per-bit correlation KCU1
cv.cd(2)
hb = ROOT.TH1D('hb', 'KCU1: time carried by each TOA bit, <t | bit = 1> - <t | bit = 0>;TOA bit;#Delta t [ns]', 10, -0.5, 9.5); keep.append(hb)
for i, b in enumerate(range(9, -1, -1)): hb.SetBinContent(b + 1, bitdt[i])
hb.SetFillColor(ROOT.kAzure - 9); hb.SetLineColor(ROOT.kAzure + 2); hb.SetMinimum(-3); hb.SetMaximum(7.5); hb.Draw('hist')
lb = ROOT.TLatex(); lb.SetNDC(); lb.SetTextFont(42); lb.SetTextSize(0.038); lb.DrawLatex(0.22, 0.8, 'bits 9:8 (160 MHz counter) carry time'); lb.DrawLatex(0.22, 0.74, 'bits 7:0 (CTDC, FTDC): < 0.2 ns each'); keep.append(lb)
# (c) ch196 ADC vs TOA model
cv.cd(3); colz()
m = c196['ok'] & np.isfinite(tt196) & (c196['A'] > 250)
h3 = ROOT.TH2D('h3', 'KCU0 N3c3 ch196: ADC 20% crossing vs calibrated TOA;TOA time [ns];ADC 20% crossing [ns]', 100, 90, 240, 100, 90, 240); keep.append(h3)
for x, y in zip(tt196[m][:200000], c196['t'][m][:200000]): h3.Fill(x, y)
h3.Draw('colz')
# (d) KCU1 ch: residual vs mid per coarse value
cv.cd(4)
k = max(good1, key=lambda kk: ch1[kk]['ok'].sum()); ch = ch1[k]; mm = ch['ok'] & (ch['A'] > 250) & (ch['A'] < 900)
base = (ch['t'] - 25.0 * ch['st'])[mm]; code = ch['code'][mm]; coarse = code >> 8; mid = (code >> 3) & 31
gs = []; lg4 = legend(0.2, 0.62, 0.6, 0.86, size=0.036); ymin, ymax = -50.0, -25.0
for q in range(4):
    g = ROOT.TGraph(); keep.append(g); gs.append(g)
    for x in range(32):
        sel = (coarse == q) & (mid == x)
        if sel.sum() > 10: y = float(np.median(base[sel])); g.AddPoint(x, y); ymin = min(ymin, y); ymax = max(ymax, y)
    g.SetMarkerStyle(20 + q); g.SetMarkerColor([ROOT.kBlack, ROOT.kRed + 1, ROOT.kBlue + 1, ROOT.kGreen + 2][q]); g.SetLineColor(g.GetMarkerColor())
    lg4.AddEntry(g, f'coarse bits = {q}', 'p')
h4 = ROOT.TH2D('h4', f'KCU1 ch{144 + k}: ADC crossing - 25 x TOA sample, by TOA fields;mid bits 7:3;median [ns]', 10, -1, 32, 10, ymin - 2, ymax + (ymax - ymin) * 0.8 + 2); keep.append(h4); h4.Draw()
for g in gs: g.Draw('PL same')
lg4.Draw(); keep.append(lg4)
# (e) KCU1 pair TOA-TOA vs ADC-ADC
cv.cd(5)
if pair_res:
    a, b, n, m1, m2, m3, dtc, dta = pair_res[0]
    h5a = ROOT.TH1D('h5a', f'KCU1 ch{144 + a} - ch{144 + b}, same waveform;t_{{a}} - t_{{b}} [ns];pairs', 100, -5, 5); h5b = ROOT.TH1D('h5b', '', 100, -5, 5); keep += [h5a, h5b]
    for d in dtc: h5a.Fill(d)
    for d in dta: h5b.Fill(d)
    h5a.SetLineColor(ROOT.kRed + 1); h5a.SetLineWidth(2); h5b.SetLineColor(ROOT.kBlue + 1); h5b.SetLineWidth(2)
    headroom(h5a, 0, max(h5a.GetMaximum(), h5b.GetMaximum())); h5a.Draw('hist'); h5b.Draw('hist same')
    lg5 = legend(0.2, 0.62, 0.8, 0.86, size=0.036); lg5.AddEntry(h5a, f'TOA, calibrated: MAD {m2:.2f} ns ({m2/np.sqrt(2)*1e3:.0f} ps/ch)', 'l'); lg5.AddEntry(h5b, f'ADC 20% crossing: MAD {m3:.2f} ns ({m3/np.sqrt(2)*1e3:.0f} ps/ch)', 'l'); lg5.Draw(); keep.append(lg5)
# (f) KCU0 ADC-based N-S per layer, cell 3 (what is available)
cv.cd(6)
hs = []; lg6 = legend(0.5, 0.62, 0.93, 0.86, size=0.036); ymax = 0
for i, k in enumerate((8, 9, 10, 11)):
    ca = per_channel(K0, k); cb = per_channel(K0, k + 16)
    m = (ca['A'] > 60) & (cb['A'] > 60) & (ca['s'] >= 5) & (cb['s'] >= 5)
    dta = ca['t'][m] - cb['t'][m]; dta = dta - np.median(dta)
    h = ROOT.TH1D(f'h6_{k}', 'KCU0 South - North from the ADC 20% crossing, cell 3;t_{S} - t_{N} [ns];pairs', 100, -12, 12); keep.append(h); hs.append(h)
    for d in dta: h.Fill(d)
    h.SetLineColor([ROOT.kBlack, ROOT.kRed + 1, ROOT.kBlue + 1, ROOT.kGreen + 2][i]); h.SetLineWidth(2); ymax = max(ymax, h.GetMaximum())
    lg6.AddEntry(h, f'SFIL{k % 4 + 1}  MAD {mad(dta):.2f} ns', 'l')
headroom(hs[0], 0, ymax); hs[0].Draw('hist')
for h in hs[1:]: h.Draw('hist same')
lg6.Draw(); keep.append(lg6)
finish(cv)
cv.SaveAs('out/toa_timing.png')
say('\nwrote out/toa_timing.png')
LOG.close()
