"""Imaging-layer performance in run 3062: cluster multiplicity per matched shower, beam profile per chip (rows), matched
fraction per paddle, hit rate per chip, and the ToT of matched hits (MIP) against all hits.
    python3 astro_imaging.py -> astro_imaging.txt/.png
"""
import numpy as np, uproot, sys, ROOT
sys.path.insert(0, '/sphenix/user/tmengel/eic/bic_testbeam/common'); from plotstyle import setup, legend, label
setup(); ROOT.gROOT.SetBatch(True); out = open('astro_imaging.txt', 'w')
def say(*a): print(*a, flush=True); print(*a, file=out)
ev = uproot.open('/sphenix/user/tmengel/eic/bic_testbeam/data/events_3062.root')['events'].arrays(['has_kcu0', 'kcu0_frag', 'pad', 'sfil_s', 'n_astro', 'astro_x', 'astro_y', 'astro_chip', 'astro_nhits', 'astro_purity', 'astro_log10LR'], library='np')
S = np.clip(ev['sfil_s'], 0, None).astype(float).sum(axis=(1, 2)); base = (ev['has_kcu0'] == 1) & (ev['kcu0_frag'] == 0) & np.isin(ev['pad'], [1, 2, 4, 8]); shower = base & (S > 600)
PB = {4: 2, 3: 1, 6: 8, 5: 4}
say('showers %d; with an AstroPix match (n_astro > 0) %d (%.3f); purity >= 0.3: %d; >= 0.7: %d' % (shower.sum(), (shower & (ev['n_astro'] > 0)).sum(), (shower & (ev['n_astro'] > 0)).mean() / shower.mean(), (shower & (ev['astro_purity'] >= 0.3)).sum(), (shower & (ev['astro_purity'] >= 0.7)).sum()))
say('per paddle: showers, matched fraction, matched with purity >= 0.3, chips of those, cluster hits (median / 90 %%), beam rows (median, 16-84 %%)')
for p in (4, 3, 6, 5):
    m = shower & (ev['pad'] == PB[p]); mm = m & (ev['astro_purity'] >= 0.3)
    chips = dict(zip(*np.unique(ev['astro_chip'][mm], return_counts=True))) if mm.sum() else {}
    nh = ev['astro_nhits'][mm]; yy = ev['astro_y'][mm]; yy = yy[np.isfinite(yy)]
    say('  P%d: %7d  matched %.3f  purity>=0.3 %6d  chips %s  nhits %s / %s  y [mm] %s' % (p, m.sum(), (m & (ev['n_astro'] > 0)).sum() / max(1, m.sum()), mm.sum(), {int(k): int(v) for k, v in chips.items()}, np.median(nh) if len(nh) else None, np.percentile(nh, 90) if len(nh) else None, (round(float(np.median(yy)), 1), round(float(np.percentile(yy, 16)), 1), round(float(np.percentile(yy, 84)), 1)) if len(yy) else None))
# per chip: matched clusters, row profile width (beam height on the chip), cluster size
mm = shower & (ev['astro_purity'] >= 0.3)
say('\nper chip (purity >= 0.3): clusters, median hits per cluster, rows: median and rms (0.5 mm pixels) -> beam height sigma [mm]')
for c in range(9):
    s_ = mm & (ev['astro_chip'] == c)
    if s_.sum() < 50: continue
    yy = ev['astro_y'][s_]; yy = yy[np.isfinite(yy)]; rows = yy / 0.5
    say('  chip %d: %6d clusters, hits/cluster median %.1f (90 %% %.0f), row median %.1f, rms %.2f rows = %.2f mm' % (c, s_.sum(), np.median(ev['astro_nhits'][s_]), np.percentile(ev['astro_nhits'][s_], 90), np.median(rows), rows.std(), 0.5 * rows.std()))
# hit records per chip over the run (all records, from the hits tree), rate per chip
t = uproot.open('/sphenix/tg/tg01/jets/tmengel/eic/astropix/astropix_3062.root')['hits']; a = t.arrays(['id', 'iscol', 'valid', 'layer'], library='np')
v = a['valid'] & (a['layer'] == 1) & (a['id'] <= 8); ids = a['id'][v]; cnt = np.bincount(ids, minlength=9)
say('\nhit records per chip over the run (%d valid, 1486 s): %s  -> row records per second %s' % (v.sum(), cnt.tolist(), (cnt / 2 / 1486.).round(0).tolist()))
c = ROOT.TCanvas('c', '', 2000, 700); c.Divide(3, 1); keep = []
c.cd(1); fr = ROOT.gPad.DrawFrame(0, 0, 10, 1.0, ';hits per matched cluster;fraction'); keep.append(fr); h = ROOT.TH1F('hn', '', 10, 0, 10); [h.Fill(x) for x in ev['astro_nhits'][mm]]; h.Scale(1. / h.GetEntries()); h.SetLineWidth(2); h.Draw('hist same'); keep.append(h); fr.SetMaximum(h.GetMaximum() * 1.3)
c.cd(2); fr2 = ROOT.gPad.DrawFrame(0, 0, 35, 0.3, ';row (0.5 mm pixels);fraction of matched clusters'); keep.append(fr2); lg = legend(0.6, 0.55, 0.88, 0.88)
for ch, colr in zip((3, 5, 6, 8), (ROOT.kBlue + 1, ROOT.kGreen + 2, ROOT.kOrange + 7, ROOT.kRed + 1)):
    s_ = mm & (ev['astro_chip'] == ch); yy = ev['astro_y'][s_]; yy = yy[np.isfinite(yy)]; hh = ROOT.TH1F('hr%d' % ch, '', 35, 0, 35); [hh.Fill(x) for x in yy / 0.5]; hh.Scale(1. / max(1, hh.GetEntries())); hh.SetLineColor(colr); hh.SetLineWidth(2); hh.Draw('hist same'); keep.append(hh); lg.AddEntry(hh, 'chip %d' % ch, 'l')
fr2.SetMaximum(max(k.GetMaximum() for k in keep[-4:]) * 1.3); lg.Draw(); keep.append(lg)
c.cd(3); fr3 = ROOT.gPad.DrawFrame(-0.5, 0, 8.5, cnt.max() * 1.2, ';chip;hit records in the run'); keep.append(fr3); hc = ROOT.TH1F('hc', '', 9, -0.5, 8.5); [hc.SetBinContent(i + 1, cnt[i]) for i in range(9)]; hc.SetLineWidth(2); hc.Draw('hist same'); keep.append(hc)
c.SaveAs('astro_imaging.png'); c.SaveAs('astro_imaging.pdf'); say('wrote astro_imaging.png')
