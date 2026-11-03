"""South-side SFIL5 / SFIL4 relative gain: KCU1 ADC (SFIL5 South, arrays B+C+D, raw peaks) per KCU0 ADC (SFIL4 South, template-fit
layer sum) for the same light, using the model only for the SFIL5/SFIL4 light ratio of adjacent layers, which depends weakly on energy.
Join: fit cache trig -> KCU0 header ts0 -> events_3062 trig_ts (same clock, same values).
    python3 sfil5_sfil4_gain.py [simdir]  -> sfil5_sfil4_gain.txt
"""
import sys, glob, numpy as np, uproot
SIMDIR = sys.argv[1] if len(sys.argv) > 1 and sys.argv[1] != 'rec' else '/sphenix/tg/tg01/jets/tmengel/eic/sim'
SUF = '_rec' if 'rec' in sys.argv[1:] else ''
out = open('sfil5_sfil4_gain%s.txt' % SUF, 'w')
def say(*a): print(*a); print(*a, file=out)
D = np.load('data_fit_layers.npz'); S = D['S'].astype(float); trig = D['trig']; pad = D['pad']; sat = D['satS'] > 0
ts0 = uproot.open('/sphenix/user/tmengel/eic/bic_testbeam/data/kcu0_3062.root')['header']['ts0'].array(library='np')
ev = uproot.open('/sphenix/user/tmengel/eic/bic_testbeam/data/events_3062.root')['events'].arrays(['trig_ts', 'has_kcu1', 'sfil5_s', 'pad'], library='np')
ets = ev['trig_ts']; order = np.argsort(ets); ets_s = ets[order]
tfit = ts0[trig]; pos = np.searchsorted(ets_s, tfit); pos = np.clip(pos, 0, len(ets_s) - 1); hit = ets_s[pos] == tfit; row = order[pos]
has5 = hit & (ev['has_kcu1'][row] == 1)
s5 = np.clip(ev['sfil5_s'][row].astype(float), 0, None); A5 = s5[:, 18:35].sum(axis=1) + s5[:, 36:53].sum(axis=1) + s5[:, 54:71].sum(axis=1)
PBIT = {3: 1, 4: 2, 5: 4, 6: 8}
say('fit-tree triggers matched to the event table: %.3f; with KCU1 present: %d' % (hit.mean(), has5.sum()))
REC = 'rec' in sys.argv[1:]
if REC:   # recovered SFIL5 sums (kcu1_sat_recovery.py: TOT / neighbour predictor for saturated single SiPMs), keyed by event row
    Rz = np.load('sfil5_recovered.npz'); arr = np.full(len(ets), np.nan); arr[Rz['row']] = Rz['A5']; A5 = arr[row]; has5 = has5 & np.isfinite(A5)
shower = (S.sum(axis=1) > 600) & has5
# model: South light SFIL4 (cells 1-3) and SFIL5 (cells 1-3), per energy, beam under P6 (x = 30) and under P5 (x = 80) where available
sim = {}
for f in sorted(glob.glob(SIMDIR + '/bic_sim_e*_s*.root')):
    g = uproot.open(f)['sim'].arrays(['ebeam', 'sfil_ls'], library='np'); E = int(round(g['ebeam'].mean())); x = 80 if '_x80' in f else 30
    if 2500 <= E <= 6500: sim[(x, E)] = g['sfil_ls'][:, :5, :3].sum(axis=2)
E6 = sorted({E for x, E in sim if x == 30})
say('\nmodel South light ratio <L5>/<L4> (cells 1-3, beam under P6): ' + '  '.join('%d GeV %.3f' % (E // 1000, sim[(30, E)][:, 4].mean() / sim[(30, E)][:, 3].mean()) for E in E6))
say('model, beam under P5 (x = 80): ' + '  '.join('%d GeV %.3f' % (E // 1000, sim[(80, E)][:, 4].mean() / sim[(80, E)][:, 3].mean()) for E in E6 if (80, E) in sim))
say('\nData <A5 (KCU1 ADC)> / <A4 (KCU0 ADC, fitted)> and the relative gain g5/g4 = ratio / model ratio, per paddle')
say('%-8s %8s %8s %8s %8s | ' % ('sample', 'n', '<A5>', '<A4>', 'A5/A4') + ' '.join('%10s' % ('g5/g4@%dGeV' % (E // 1000)) for E in E6))
rows = []
for name, m in [('all', shower)] + [('P%d' % p, shower & (pad == p)) for p in (4, 3, 6, 5)] + [('all,nosat4', shower & ~sat[:, 3])] + [('P%d,nosat4' % p, shower & (pad == p) & ~sat[:, 3]) for p in (4, 3, 6, 5)]:
    if m.sum() < 100: continue
    a5, a4 = A5[m].mean(), S[m, 3].mean(); r = a5 / a4
    x = 80 if 'P5' in name else 30
    gs = [r / (sim[(x, E)][:, 4].mean() / sim[(x, E)][:, 3].mean()) if (x, E) in sim else np.nan for E in E6]
    say('%-10s %8d %8.0f %8.1f %8.3f | ' % (name, m.sum(), a5, a4, r) + ' '.join('%10.2f' % g for g in gs)); rows.append((name, m.sum(), r, gs))
say('\nNote: the "nosat4" rows drop showers whose SFIL4 South beam-row channel saturates (a selection against large SFIL4 pulses; the model\nratio is not re-cut, so those rows are a cross-check of the fit recovery, not a calibration). SFIL5 South single SiPMs DO saturate (ceiling ~943 ADC above pedestal; an earlier check used a 1000 ADC threshold on pedestal-subtracted peaks and missed it).')
out.close()

# ---- figure: per-event A5/A4 distributions per paddle (left) and g5/g4 versus the energy hypothesis (right)
import ROOT
sys.path.insert(0, '/sphenix/user/tmengel/eic/bic_testbeam/common'); from plotstyle import setup, legend, label
setup(); ROOT.gStyle.SetPadLeftMargin(0.15); ROOT.gStyle.SetPadRightMargin(0.04); ROOT.gStyle.SetPadTopMargin(0.1); ROOT.gStyle.SetPadBottomMargin(0.15); ROOT.gStyle.SetTitleOffset(1.1, 'x'); ROOT.gROOT.ForceStyle()
c = ROOT.TCanvas('c', '', 1500, 650); c.Divide(2, 1, 0.002, 0.002); keep = []
COL = {4: ROOT.kAzure + 2, 3: ROOT.kGreen + 2, 6: ROOT.kOrange + 7, 5: ROOT.kRed + 1}
c.cd(1); hs = {}; hmax = 0
for p in (4, 3, 6, 5):
    m = shower & (pad == p) & (S[:, 3] > 150); r = A5[m] / S[m, 3]
    h = ROOT.TH1D('h%d' % p, '', 60, 0, 40); [h.Fill(v) for v in r[:200000]]; h.Scale(1. / h.Integral()); h.SetLineColor(COL[p]); h.SetLineWidth(3); hs[p] = h; hmax = max(hmax, h.GetMaximum()); keep.append(h)
fr = ROOT.gPad.DrawFrame(0, 0, 40, hmax * 1.9, ';A_{5} / A_{4} per shower  (KCU1 ADC / KCU0 ADC, South);fraction of showers'); fr.GetYaxis().SetTitleOffset(1.4); fr.GetXaxis().SetNdivisions(505); fr.GetYaxis().SetNdivisions(505); keep.append(fr)
lg = legend(0.5, 0.6, 0.93, 0.88, size=0.04); lg.SetHeader('run 3062, SFIL4 South > 150 ADC')
for p in (4, 3, 6, 5):
    hs[p].Draw('hist same'); m = shower & (pad == p); lg.AddEntry(hs[p], 'P%d: <A5>/<A4> = %.2f' % (p, A5[m].mean() / S[m, 3].mean()), 'l')
lg.Draw(); keep.append(lg)
c.cd(2); fr2 = ROOT.gPad.DrawFrame(2.5, 8, 6.5, 22, ';beam energy hypothesis [GeV];g_{5} / g_{4}  [KCU1 ADC per KCU0 ADC]'); fr2.GetYaxis().SetTitleOffset(1.3); fr2.GetXaxis().SetNdivisions(505); fr2.GetYaxis().SetNdivisions(505); keep.append(fr2)
lg2 = legend(0.5, 0.62, 0.93, 0.88, size=0.04); lg2.SetHeader('data ratio / model light ratio')
for p in (4, 3, 6, 5):
    m = shower & (pad == p); r = A5[m].mean() / S[m, 3].mean(); x = 80 if p == 5 else 30
    Es = [E for E in E6 if (x, E) in sim]; gs = np.array([r / (sim[(x, E)][:, 4].mean() / sim[(x, E)][:, 3].mean()) for E in Es]); g = ROOT.TGraph(len(Es), np.array(Es, float) / 1e3, gs); g.SetLineColor(COL[p]); g.SetMarkerColor(COL[p]); g.SetLineWidth(3); g.SetMarkerStyle(20); g.SetMarkerSize(1.6); g.Draw('LP'); keep.append(g); lg2.AddEntry(g, 'P%d showers' % p, 'lp')
lg2.Draw(); keep.append(lg2)
c.cd(); keep.append(label('SFIL5 / SFIL4 South relative gain: the model supplies only the light ratio of two adjacent layers (0.77-0.87 from 3 to 6 GeV)', x=0.01, y=0.95, size=0.032))
c.SaveAs('sfil5_sfil4_gain%s.png' % SUF); c.SaveAs('sfil5_sfil4_gain%s.pdf' % SUF); print('wrote sfil5_sfil4_gain.png')
