"""Trigger paddles P3-P6 as seen by the CAEN (channels 12-15, cell-corrected) and by KCU0 ASIC1 (channels 84/86/81/97).
Pulse-height spectrum per paddle for single-paddle triggers, most-probable value (Landau x Gauss fit around the peak),
clipping at the CAEN floor, charge/peak (pulse shape), timing; KCU0: negative pulses on a 77 ADC pedestal that clip at
ADC 0, TOA coverage, no TOT.  If the four paddles are identical, the MPVs are their relative gains.
    python3 paddle_gains.py  -> paddle_gains.txt, paddle_gains.png, paddle_waveforms.png
"""
import numpy as np, uproot, sys, ROOT
sys.path.insert(0, '/sphenix/user/tmengel/eic/bic_testbeam/common'); from plotstyle import setup, legend, label
setup(); out = open('paddle_gains.txt', 'w')
def say(*a): print(*a, flush=True); print(*a, file=out)
segs = [np.load('/sphenix/tg/tg01/jets/tmengel/eic/caen/calib/caen_calib_seg%d.npz' % s) for s in (0, 1, 2)]
cat = lambda k: np.concatenate([s[k] for s in segs])
peak, qn, q140, tpk, base = cat('peak'), cat('qn'), cat('q140'), cat('tpk'), cat('base')
names = ['P3', 'P4', 'P5', 'P6']; P = peak[:, 12:16]; Q = q140[:, 12:16]; T = tpk[:, 12:16]
fired = P > 150; nf = fired.sum(axis=1)
say('CAEN paddles, %d triggers; paddles above 150 counts per trigger: %s' % (len(P), np.bincount(nf, minlength=5)))
# CAEN pedestal of the paddle channels ~3730-3780, so the floor (ADC 0) is a depth of ~3750: pulses with peak > 3600 are clipped
FLOOR = 3600.
ROOT.gROOT.SetBatch(True)
c = ROOT.TCanvas('c', '', 2400, 1600); c.Divide(3, 2); keep = []
cols = {'P3': ROOT.kBlue + 1, 'P4': ROOT.kGreen + 2, 'P5': ROOT.kOrange + 7, 'P6': ROOT.kRed + 1}
c.cd(1); fr = ROOT.gPad.DrawFrame(0, 0, 4000, 0.06, ';CAEN peak depth [counts];fraction of single-paddle triggers per 50 counts'); keep.append(fr); lg = legend(0.45, 0.6, 0.88, 0.88)
mpv = {}
fl = ROOT.TF1('fl', 'landau', 0, 4000)
for i, nm in enumerate(names):
    m = fired[:, i] & (nf == 1); v = P[m, i]
    h = ROOT.TH1F('h' + nm, '', 80, 0, 4000); [h.Fill(x) for x in v]; h.Scale(1. / h.GetEntries()); h.SetLineColor(cols[nm]); h.SetLineWidth(2); h.Draw('hist same'); keep.append(h)
    # MPV: fit a Landau around the mode of the unclipped part
    hh = ROOT.TH1F('hh' + nm, '', 160, 0, 4000); [hh.Fill(x) for x in v[v < FLOOR]]; j = hh.GetMaximumBin(); x0 = hh.GetBinCenter(j)
    fl.SetParameters(hh.GetMaximum(), x0, 0.15 * x0); hh.Fit(fl, 'QN', '', max(200, 0.6 * x0), min(FLOOR, 1.6 * x0)); mpv[nm] = fl.GetParameter(1)
    clip = (v >= FLOOR).mean(); qp = np.median(Q[m, i] / np.maximum(P[m, i], 1))
    say('%s: single-paddle triggers %7d; mode %4.0f, Landau MPV %4.0f (width %3.0f); median %4.0f; clipped at the floor %.3f; charge/peak median %.1f ns; peak sample median %.0f' % (nm, m.sum(), x0, mpv[nm], fl.GetParameter(2), np.median(v), clip, qp, np.median(T[m, i])))
    lg.AddEntry(h, '%s: MPV %.0f, clipped %.0f %%' % (nm, mpv[nm], 100 * clip), 'l')
lg.Draw(); keep.append(lg); keep.append(label('single-paddle triggers, DRS4 cell-corrected', x=0.15, y=0.93, size=0.035))
ref = mpv['P3']; say('relative gains from the MPV (P3 = 1): ' + '  '.join('%s %.2f' % (nm, mpv[nm] / ref) for nm in names))
# two-paddle triggers: the same pulse in both? (adjacent paddles overlap?) and the pulse height in the second paddle
say('\ntwo-paddle triggers by pair (P3 P4 P5 P6 order is North to South: P4 P3 P6 P5):')
for i in range(4):
    for j in range(i + 1, 4):
        m = fired[:, i] & fired[:, j] & (nf == 2)
        if m.sum() > 100: say('  %s+%s: %6d  median peaks %4.0f / %4.0f' % (names[i], names[j], m.sum(), np.median(P[m, i]), np.median(P[m, j])))
# charge spectrum (q140) as the second estimator, less sensitive to clipping
c.cd(2); fr = ROOT.gPad.DrawFrame(0, 0, 200000, 0.06, ';CAEN charge, samples 20-159 [counts x ns];fraction per bin'); keep.append(fr); lg2 = legend(0.45, 0.65, 0.88, 0.88)
for i, nm in enumerate(names):
    m = fired[:, i] & (nf == 1); v = Q[m, i]; h = ROOT.TH1F('hq' + nm, '', 80, 0, 200000); [h.Fill(x) for x in v]; h.Scale(1. / h.GetEntries()); h.SetLineColor(cols[nm]); h.SetLineWidth(2); h.Draw('hist same'); keep.append(h); lg2.AddEntry(h, '%s: median %.0f' % (nm, np.median(v)), 'l')
lg2.Draw(); keep.append(lg2)
# peak vs charge (clipping bends it)
c.cd(3); fr = ROOT.gPad.DrawFrame(0, 0, 4000, 200000, ';CAEN peak depth [counts];CAEN charge [counts x ns]'); keep.append(fr)
m = fired[:, 3] & (nf == 1); h2 = ROOT.TH2F('h2', '', 80, 0, 4000, 80, 0, 200000); [h2.Fill(x, y) for x, y in zip(P[m, 3][::10], Q[m, 3][::10])]; h2.Draw('col same'); keep.append(h2); keep.append(label('P6: charge against peak; the vertical edge at the floor is clipping', x=0.15, y=0.93, size=0.032))
# peak time per paddle
c.cd(4); fr = ROOT.gPad.DrawFrame(60, 0, 120, 0.3, ';CAEN peak sample [ns];fraction'); keep.append(fr); lg4 = legend(0.6, 0.65, 0.88, 0.88)
for i, nm in enumerate(names):
    m = fired[:, i] & (nf == 1); h = ROOT.TH1F('ht' + nm, '', 60, 60, 120); [h.Fill(x) for x in T[m, i]]; h.Scale(1. / h.GetEntries()); h.SetLineColor(cols[nm]); h.SetLineWidth(2); h.Draw('hist same'); keep.append(h); lg4.AddEntry(h, nm, 'l')
lg4.Draw(); keep.append(lg4)
# ---- KCU0 paddles: 200k triggers
t = uproot.open('/sphenix/user/tmengel/eic/bic_testbeam/data/kcu0_3062.root')['paddles']
a = t.arrays([k for k in t.keys() if k.endswith('_adc') or k.endswith('_toa') or k.endswith('_tot')], entry_stop=200000, library='np')
say('\nKCU0 paddles (ASIC1, 200k triggers): negative pulses on a ~77 ADC pedestal')
c.cd(5); fr = ROOT.gPad.DrawFrame(0, 0, 120, 0.5, ';KCU0 depth below pedestal [ADC] (clipped at the pedestal, ~77);fraction'); keep.append(fr); lg5 = legend(0.4, 0.6, 0.88, 0.88)
for nm in names:
    W = a[nm + '_adc'].astype(float); g = W[:, 0] >= 0; W = W[g]; ped = W[:, 1:4].mean(axis=1); d = ped - W.min(axis=1); clip = (W.min(axis=1) <= 0).mean()
    toa = (a[nm + '_toa'][g] > 0).any(axis=1).mean(); tot = (a[nm + '_tot'][g] > 0).any(axis=1).mean()
    say('  %s: clipped at ADC 0 in %.3f of triggers, TOA word in %.3f, TOT word in %.3f; depth median %.0f' % (nm, clip, toa, tot, np.median(d)))
    h = ROOT.TH1F('hk' + nm, '', 60, 0, 120); [h.Fill(x) for x in d]; h.Scale(1. / h.GetEntries()); h.SetLineColor(cols[nm]); h.SetLineWidth(2); h.Draw('hist same'); keep.append(h); lg5.AddEntry(h, '%s: clipped %.0f %%' % (nm, 100 * clip), 'l')
lg5.Draw(); keep.append(lg5)
c.cd(6); fr = ROOT.gPad.DrawFrame(0, -10, 19, 110, ';sample (25 ns);ADC'); keep.append(fr)
W = a['P5_adc'][:2000].astype(float); shown = 0
for i in range(len(W)):
    if W[i, 0] < 0: continue
    gr = ROOT.TGraph(20, np.arange(20, dtype=float), W[i]); gr.SetLineColor(ROOT.kGray + 1 + (shown % 3)); gr.Draw('l same'); keep.append(gr); shown += 1
    if shown >= 40: break
keep.append(label('KCU0 P5: 40 raw waveforms, negative pulses clipped at ADC 0', x=0.15, y=0.93, size=0.032))
c.SaveAs('paddle_gains.png'); c.SaveAs('paddle_gains.pdf'); say('wrote paddle_gains.png')
# ---- CAEN waveform gallery per paddle (seg1, first 20k events, cell-corrected)
f = uproot.open('/sphenix/tg/tg01/jets/tmengel/eic/caen/caen_3062_seg1.root'); N = 20000
hd = f['header'].arrays(['g_index_cell'], entry_stop=N, library='np')['g_index_cell']; pdw = f['paddles'].arrays(['P3_adc', 'P4_adc', 'P5_adc', 'P6_adc'], entry_stop=N, library='np')
ped = np.load('/sphenix/user/tmengel/eic/bic_testbeam/calibration/caen_calib/drs4_pedestals.npz')['ped']
cw = ROOT.TCanvas('cw', '', 2400, 700); cw.Divide(4, 1); keep2 = []
for i, nm in enumerate(names):
    cw.cd(i + 1); fr = ROOT.gPad.DrawFrame(0, -200, 256, 4000, '%s on the CAEN;sample [ns];baseline - ADC [counts]' % nm); keep2.append(fr)
    Wp = pdw[nm + '_adc'].astype(float); cells = (hd[:, 1, None] + np.arange(256)[None, :]) % 1024; corr = Wp - ped[12 + i][cells]; b = np.median(corr[:, :20], axis=1); D = b[:, None] - corr
    pk = D[:, 20:230].max(axis=1); idx = np.where(pk > 150)[0][:30]
    for k, j in enumerate(idx):
        gr = ROOT.TGraph(256, np.arange(256, dtype=float), D[j]); gr.SetLineColor(ROOT.kGray + 1 + (k % 3)); gr.Draw('l same'); keep2.append(gr)
cw.SaveAs('paddle_waveforms.png'); say('wrote paddle_waveforms.png')
