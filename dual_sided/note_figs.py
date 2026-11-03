"""Extra figures for the dual-sided note (mirroring the 2017 EmCal note's Figs. 3, 4, 5 and 8), SFIL2 cell 3, run 3062.
    python3 note_figs.py -> fig_waveform.pdf, fig_clipped.pdf, fig_selection.pdf, fig_cfd.pdf, note_figs.txt"""
import numpy as np, uproot, sys, ROOT
sys.path.insert(0, '/sphenix/user/tmengel/eic/bic_testbeam/common'); from plotstyle import setup, legend, headroom, finish
setup(); out = open('note_figs.txt', 'w'); keep = []
def say(*a): print(*a, flush=True); print(*a, file=out)
NMAX = 400000
f = uproot.open('/sphenix/user/tmengel/eic/bic_testbeam/data/kcu0_3062.root'); ts0 = f['header']['ts0'].array(library='np', entry_stop=NMAX)
ev = uproot.open('/sphenix/user/tmengel/eic/bic_testbeam/data/events_3062.root')['events'].arrays(['trig_ts', 'pad'], library='np')
order = np.argsort(ev['trig_ts']); es = ev['trig_ts'][order]; pos = np.clip(np.searchsorted(es, ts0), 0, len(es) - 1); hit = es[pos] == ts0
pad = np.where(hit, ev['pad'][order[pos]], 0)
W = {e: f['sfil2'][br + '_adc'].array(library='np', entry_stop=NMAX).astype(float) for e, br in (('S', 'Sc3'), ('N', 'Nc3'))}
good = (W['S'][:, 0] >= 0) & (W['N'][:, 0] >= 0); one = np.isin(pad, [1, 2, 4, 8])
Y = {e: W[e] - W[e][:, 1:4].mean(axis=1)[:, None] for e in W}; RAW = {e: W[e].max(axis=1) for e in W}
A = {e: Y[e][np.arange(len(pad)), Y[e][:, 4:16].argmax(axis=1) + 4] for e in W}
tns = 25. * np.arange(20)
def wf_graph(y, col, mk):
    g = ROOT.TGraph(20, tns, np.ascontiguousarray(y)); g.SetLineColor(col); g.SetMarkerColor(col); g.SetMarkerStyle(mk); g.SetMarkerSize(1.5); g.SetLineWidth(3); keep.append(g); return g
# Fig. waveform: one unclipped P6 shower with both ends 400-700 ADC
tail = lambda e: np.all(np.diff(Y[e][:, 10:], axis=1) < 15, axis=1)   # tail falls, no second pulse
cand = np.where(good & (pad == 8) & (A['S'] > 400) & (A['S'] < 700) & (A['N'] > 400) & (A['N'] < 700) & (RAW['S'] < 1000) & (RAW['N'] < 1000) & tail('S') & tail('N'))[0]; i = cand[5]
c = ROOT.TCanvas('c1', '', 1000, 750)
fr = ROOT.gPad.DrawFrame(0, -60, 500, 1, 'SFIL2 cell 3 waveforms, one shower;time [ns];amplitude [ADC, pedestal subtracted]'); headroom(fr, -60, max(Y['S'][i].max(), Y['N'][i].max()), 0.30); keep.append(fr)
lg = legend(0.55, 0.74, 0.92, 0.86)
for e, col, mk, nm in (('S', ROOT.kBlue + 1, 20, 'South end'), ('N', ROOT.kMagenta + 1, 21, 'North end')): g = wf_graph(Y[e][i], col, mk); g.Draw('pl same'); lg.AddEntry(g, nm, 'pl')
lg.Draw(); finish(c); c.SaveAs('fig_waveform.png'); say('waveform example: trigger %d, paddle P6, A_S %.0f, A_N %.0f' % (i, A['S'][i], A['N'][i]))
# Fig. clipped: both ends reach the 10-bit ceiling
cand = np.where(good & one & (RAW['S'] >= 1023) & (RAW['N'] >= 1023) & (W['S'][:, 1:4].max(axis=1) < 300) & (W['N'][:, 1:4].max(axis=1) < 300) & tail('S') & tail('N'))[0]; i = cand[3]
c = ROOT.TCanvas('c2', '', 1000, 750)
fr = ROOT.gPad.DrawFrame(0, -100, 500, 1, 'SFIL2 cell 3 waveforms, clipped;time [ns];amplitude [ADC, raw]'); headroom(fr, -100, 1023, 0.25); keep.append(fr)
lg = legend(0.55, 0.78, 0.92, 0.88)
for e, col, mk, nm in (('S', ROOT.kBlue + 1, 20, 'South end'), ('N', ROOT.kMagenta + 1, 21, 'North end')): g = wf_graph(W[e][i], col, mk); g.Draw('pl same'); lg.AddEntry(g, nm, 'pl')
ln = ROOT.TLine(0, 1023, 500, 1023); ln.SetLineStyle(2); ln.Draw(); keep.append(ln)
lg.Draw(); finish(c); c.SaveAs('fig_clipped.png')
# Fig. selection: raw maximum of each end, with the clipping cut
c = ROOT.TCanvas('c3', '', 1000, 750); ROOT.gPad.SetLogy()
hs = []
for e, col, st in (('S', ROOT.kBlue + 1, 1), ('N', ROOT.kMagenta + 1, 1)):
    h = ROOT.TH1D('hr' + e, '', 105, 0, 1050); v = np.ascontiguousarray(RAW[e][good & one]); h.FillN(len(v), v, ROOT.nullptr); h.SetLineColor(col); h.SetLineWidth(3); h.SetLineStyle(st); hs.append(h); keep.append(h)
fr = ROOT.gPad.DrawFrame(0, 1, 1050, 10, 'raw maximum per waveform, single-paddle triggers;largest raw sample [ADC];waveforms'); headroom(fr, 1, max(h.GetMaximum() for h in hs), 0.35, log=True); keep.append(fr)
lg = legend(0.20, 0.74, 0.60, 0.86)
for h, nm in zip(hs, ('South end', 'North end')): h.Draw('hist same'); lg.AddEntry(h, nm, 'l')
ln = ROOT.TLine(1000, 1, 1000, max(h.GetMaximum() for h in hs)); ln.SetLineColor(ROOT.kRed + 1); ln.SetLineStyle(2); ln.SetLineWidth(3); ln.Draw(); keep.append(ln); lg.AddEntry(ln, 'cut: below 1000', 'l')
lg.Draw(); finish(c); c.SaveAs('fig_selection.png')
m = good & one & (A['S'] > 100) & (A['N'] > 100)
say('single-paddle triggers with both ends above 100 ADC: %d; of these clipped at either end: %.3f' % (m.sum(), ((RAW['S'] >= 1000) | (RAW['N'] >= 1000))[m].mean()))
# Fig. CFD: waveform, CFD function, zoom on the zero crossing (South end of the waveform example)
cand = np.where(good & (pad == 8) & (A['S'] > 400) & (A['S'] < 700) & (RAW['S'] < 1000) & tail('S'))[0]; i = cand[5]; y = Y['S'][i]
cf = y[:-1] - 0.8 * y[1:]; k0 = np.argmax(y)
j = next(j for j in range(2, 18) if cf[j] < 0 and cf[j + 1] >= 0 and j >= k0 - 3); tz = 25. * (j + (-cf[j]) / (cf[j + 1] - cf[j]))
say('CFD example: zero crossing between samples %d and %d at %.2f ns; peak sample %d' % (j, j + 1, tz, k0))
c = ROOT.TCanvas('c4', '', 1600, 1200)
P = [ROOT.TPad('p1', '', 0, 0.5, 0.5, 1), ROOT.TPad('p2', '', 0.5, 0.5, 1, 1), ROOT.TPad('p3', '', 0.2, 0, 0.8, 0.5)]
for p in P: c.cd(); p.Draw(); keep.append(p)
P[0].cd(); fr = ROOT.gPad.DrawFrame(0, -60, 500, y.max() * 1.25, 'South waveform;time [ns];amplitude [ADC]'); keep.append(fr); wf_graph(y, ROOT.kBlue + 1, 20).Draw('pl same')
P[1].cd(); fr = ROOT.gPad.DrawFrame(0, cf.min() * 1.3, 500, cf.max() * 1.3, 'CFD signal c[n] = w[n] - 0.8 w[n+1];time [ns];CFD signal [ADC]'); keep.append(fr)
g = ROOT.TGraph(19, tns[:19], np.ascontiguousarray(cf)); g.SetMarkerStyle(20); g.SetMarkerSize(1.5); g.SetLineWidth(3); g.SetLineColor(ROOT.kBlue + 1); g.SetMarkerColor(ROOT.kBlue + 1); g.Draw('pl same'); keep.append(g)
z0 = ROOT.TLine(0, 0, 500, 0); z0.SetLineStyle(3); z0.Draw(); keep.append(z0)
P[2].cd(); lo, hi = 25. * j - 12, 25. * (j + 1) + 12
fr = ROOT.gPad.DrawFrame(lo, min(cf[j], 0) * 1.6 - 5, hi, max(cf[j + 1], 0) * 1.6 + 5, 'CFD signal, zoom on the zero crossing;time [ns];CFD signal [ADC]'); keep.append(fr)
g.Draw('p same'); fl = ROOT.TLine(25. * j, cf[j], 25. * (j + 1), cf[j + 1]); fl.SetLineColor(ROOT.kRed + 1); fl.SetLineWidth(3); fl.Draw(); keep.append(fl)
z1 = ROOT.TLine(lo, 0, hi, 0); z1.SetLineStyle(3); z1.Draw(); keep.append(z1)
mk = ROOT.TMarker(tz, 0, 29); mk.SetMarkerColor(ROOT.kRed + 1); mk.SetMarkerSize(0); keep.append(mk)
tl = ROOT.TLine(tz, fr.GetMinimum(), tz, fr.GetMaximum()); tl.SetLineColor(ROOT.kRed + 1); tl.SetLineStyle(2); tl.SetLineWidth(2); tl.Draw(); keep.append(tl)
lg = legend(0.56, 0.20, 0.95, 0.32, size=0.045); lg.AddEntry(fl, 'linear interpolation', 'l'); lg.AddEntry(tl, 'zero crossing %.1f ns' % tz, 'l'); lg.Draw(); keep.append(lg)
from plotstyle import label
for p in P:
    p.cd(); t = [o for o in p.GetListOfPrimitives() if isinstance(o, ROOT.TH1)][0].GetTitle(); label(t, size=0.05, pad=p)
c.SaveAs('fig_cfd.png')
