"""Figures and LaTeX tables for the decoding / synchronization note from output/sync_runs/run<RUN>/{sync,decode}.npz.

    python3 run_sync_figs.py OUTDIR
writes OUTDIR/sync_fig_{kcu1,caen,paddle,astro}.png and OUTDIR/tab_{decode,sync,astro}.tex
"""
import sys, os, numpy as np, ROOT
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import sys; sys.path.insert(0, '/sphenix/user/tmengel/eic/bic_testbeam/common')
from plotstyle import setup, finish, legend, headroom
setup()
OUT = sys.argv[1]; os.makedirs(OUT, exist_ok=True)
S = '/sphenix/user/tmengel/eic/output/sync_runs/run%d/%s.npz'
RUNS = [3062, 3061, 3063, 3064, 3065, 3066, 3067, 3068, 3069, 3088, 3090, 3091]
syn, dec = {}, {}
for r in RUNS:
    if os.path.exists(S % (r, 'sync')): syn[r] = dict(np.load(S % (r, 'sync'), allow_pickle=True))
    if os.path.exists(S % (r, 'decode')): dec[r] = np.load(S % (r, 'decode'), allow_pickle=True)['R'].item()
print('runs with sync:', sorted(syn), ' with decode:', sorted(dec))
COL = {3062: ROOT.kBlack, 3063: ROOT.kBlue + 1, 3064: ROOT.kRed + 1, 3069: ROOT.kGreen + 2, 3088: ROOT.kMagenta + 1, 3066: ROOT.kOrange + 7}
MRK = {3062: 20, 3063: 21, 3064: 22, 3069: 33, 3088: 23, 3066: 34}
keep = []
def H(name, title, edges, vals):
    h = ROOT.TH1D(name, title, len(edges) - 1, np.asarray(edges, float)); keep.append(h)
    for i, v in enumerate(vals): h.SetBinContent(i + 1, float(v))
    return h

# ---------------------------------------------------------------- KCU1
cv = ROOT.TCanvas('c1', '', 1500, 650); cv.Divide(2, 1, 0.002, 0.002); keep.append(cv)
cv.cd(1); ROOT.gPad.SetLogy()
show = [r for r in (3062, 3063, 3088) if r in syn and 'off1_spike' in syn[r]]
leg = legend(0.55, 0.70, 0.93, 0.88); hs = []
for r in show:
    h = H('sp%d' % r, ';ts(KCU0) - ts(KCU1) - coarse offset [ms];pairs / 10 #mus', np.linspace(-20, 20, 4001), syn[r]['off1_spike'])
    h.Rebin(10); h.Scale(1 / h.Integral()); h.GetYaxis().SetTitle('fraction of pairs / 100 #mus'); h.SetLineColor(COL[r]); hs.append(h); leg.AddEntry(h, 'run %d' % r, 'l')
if hs:
    headroom(hs[0], 1e-6, max(h.GetMaximum() for h in hs), 0.35, log=True); hs[0].Draw('hist')
    for h in hs[1:]: h.Draw('hist same')
    leg.Draw()
cv.cd(2)
fr = ROOT.TH1D('fr1', ';ts(KCU0) - ts(KCU1) - offset [6.25 ns ticks];fraction of matched pairs', 17, -8.5, 8.5); keep.append(fr)
fr.SetMaximum(1.35); fr.SetMinimum(0); fr.Draw('axis')
leg2 = legend(0.22, 0.70, 0.60, 0.88)
for k, r in enumerate(show):
    v = syn[r]['off1_res'].astype(float); v /= v.sum()
    g = ROOT.TGraph(17, np.arange(-8, 9, dtype=float) + 0.12 * (k - 1), v); g.SetTitle(''); g.SetMarkerStyle(MRK[r]); g.SetMarkerColor(COL[r]); g.SetMarkerSize(1.5); g.Draw('P same'); keep.append(g)
    leg2.AddEntry(g, 'run %d' % r, 'p')
leg2.Draw(); finish(cv); cv.SaveAs(OUT + '/sync_fig_kcu1.png')

# ---------------------------------------------------------------- CAEN
cv = ROOT.TCanvas('c2', '', 2100, 650); cv.Divide(3, 1, 0.002, 0.002); keep.append(cv)
cv.cd(1)
long_runs = sorted([r for r in syn if 'caen_slice_t' in syn[r]], key=lambda r: -syn[r]['caen_slice_t'].max())[:3]
fr = ROOT.TH1F('fr2', ';time in run [h];slice offset - line [#mus]', 10, 0, max([syn[r]['caen_slice_t'].max() / 3600 for r in long_runs] + [0.5]) * 1.05); keep.append(fr)
ym = max([np.abs(syn[r]['caen_slice_res_us']).max() for r in long_runs] + [5]); fr.SetMinimum(-1.3 * ym); fr.SetMaximum(2.0 * ym); fr.Draw('axis')
leg = legend(0.22, 0.72, 0.9, 0.88)
for r in long_runs:
    g = ROOT.TGraph(len(syn[r]['caen_slice_t']), syn[r]['caen_slice_t'] / 3600, syn[r]['caen_slice_res_us'].astype(float)); g.SetMarkerStyle(MRK.get(r, 20))
    g.SetMarkerColor(COL.get(r, 1)); g.SetMarkerSize(0.7); g.SetTitle(''); g.Draw('P same'); keep.append(g); leg.AddEntry(g, 'run %d, drift %+.2f ppm' % (r, syn[r]['caen_drift_ppm']), 'p')
leg.Draw()
ref = 3063 if 3063 in syn else sorted(syn)[0]
cv.cd(2); ROOT.gPad.SetLogy()
h1 = H('cr', ';KCU0 - predicted CAEN time [#mus];CAEN triggers / 50 ns', np.linspace(-5, 5, 201), syn[ref]['caen_res_hist'])
h2 = H('cc', '', np.linspace(-5, 5, 201), syn[ref]['caen_ctl_hist'])
h1.SetLineColor(ROOT.kBlue + 1); h2.SetLineColor(ROOT.kRed + 1); h2.SetLineStyle(2)
headroom(h1, 0.5, max(h1.GetMaximum(), h2.GetMaximum()), 0.35, log=True); h1.Draw('hist'); h2.Draw('hist same')
leg = legend(0.22, 0.74, 0.92, 0.88); leg.AddEntry(h1, 'run %d, nearest KCU0 trigger' % ref, 'l'); leg.AddEntry(h2, 'CAEN shifted +10 ms (control)', 'l'); leg.Draw()
cv.cd(3)
h3 = H('cz', ';KCU0 - predicted CAEN time [ns];fraction / ns', np.linspace(-50, 50, 101), syn[ref]['caen_zoom_hist']); h3.Scale(1 / h3.Integral())
h3.SetLineColor(ROOT.kBlue + 1); headroom(h3, 0, h3.GetMaximum(), 0.3); h3.GetXaxis().SetNdivisions(505); h3.Draw('hist')
leg = legend(0.22, 0.78, 0.92, 0.88); leg.AddEntry(h3, 'run %d, MAD %.1f ns' % (ref, syn[ref]['caen_mad_ns']), 'l'); leg.Draw()
finish(cv); cv.SaveAs(OUT + '/sync_fig_caen.png')

# ---------------------------------------------------------------- paddle patterns
pr = [r for r in RUNS if r in syn and 'pad_same' in syn[r]]
cv = ROOT.TCanvas('c3', '', 1100, 650); keep.append(cv)
fr = ROOT.TH1F('fr3', ';;KCU0 and CAEN paddle pattern equal', len(pr), 0, len(pr)); keep.append(fr)
for i, r in enumerate(pr): fr.GetXaxis().SetBinLabel(i + 1, str(r))
fr.SetMinimum(0); fr.SetMaximum(1.45); fr.GetXaxis().SetLabelSize(0.05); fr.Draw('axis')
g1 = ROOT.TGraph(len(pr), np.arange(len(pr)) + 0.5, np.array([syn[r]['pad_same'] for r in pr], float))
g2 = ROOT.TGraph(len(pr), np.arange(len(pr)) + 0.5, np.array([syn[r]['pad_decoy'] for r in pr], float))
g1.SetMarkerStyle(20); g1.SetMarkerSize(1.6); g1.SetMarkerColor(ROOT.kBlue + 1); g2.SetMarkerStyle(25); g2.SetMarkerSize(1.6); g2.SetMarkerColor(ROOT.kRed + 1)
g1.SetTitle(''); g2.SetTitle(''); g1.Draw('P same'); g2.Draw('P same'); keep += [g1, g2]
leg = legend(0.18, 0.76, 0.9, 0.88, ncol=2); leg.AddEntry(g1, 'same trigger', 'p'); leg.AddEntry(g2, 'neighbouring trigger (decoy)', 'p'); leg.Draw()
finish(cv); cv.SaveAs(OUT + '/sync_fig_paddle.png')

# ---------------------------------------------------------------- AstroPix
ar = [r for r in (3062, 3063, 3064, 3069, 3066) if r in syn and 'astro_z' in syn[r]]
cv = ROOT.TCanvas('c4', '', 1500, 650); cv.Divide(2, 1, 0.002, 0.002); keep.append(cv)
cv.cd(1)
fr = ROOT.TH1F('fr4', ';phase - peak phase;excess over flat [#sigma]', 10, -12.5, 12.5); keep.append(fr)
zm = max([syn[r]['astro_z'].max() for r in ar] + [5]); zl = min([syn[r]['astro_z'].min() for r in ar] + [-1])
fr.SetMinimum(zl - 0.05 * zm); fr.SetMaximum(1.5 * zm); fr.Draw('axis')
leg = legend(0.18, 0.70, 0.93, 0.88, ncol=2)
for k, r in enumerate(ar):
    z = syn[r]['astro_z']; ph = int(syn[r]['astro_phase']); x = np.arange(-12, 13)
    g = ROOT.TGraph(len(x), x + 0.12 * (k - 2.0), z[(ph + x) % 256].astype(float)); g.SetMarkerStyle(MRK.get(r, 20)); g.SetMarkerColor(COL.get(r, 1)); g.SetLineColor(COL.get(r, 1))
    g.SetMarkerSize(1.2); g.SetTitle(''); g.Draw('PL same'); keep.append(g); leg.AddEntry(g, '%d: phase %d' % (r, ph), 'p')
leg.Draw()
cv.cd(2)
lr = [r for r in ar if 'astro_lat_m' in syn[r]]
fr = ROOT.TH1F('fr5', ';K - AstroPix stamp - O [ms];(matched - control) / matched total', 10, -250, 250); keep.append(fr); fr.GetXaxis().SetNdivisions(505)
mx = 0; gs = []
for r in lr:
    e = syn[r]['astro_lat_edges']; ex = (syn[r]['astro_lat_m'] - syn[r]['astro_lat_x']).astype(float); ex /= max(ex.sum(), 1)
    h = H('lat%d' % r, '', e, ex); h.SetLineColor(COL.get(r, 1)); gs.append((r, h)); mx = max(mx, h.GetMaximum())
fr.SetMinimum(-0.05 * mx); fr.SetMaximum(1.45 * mx); fr.Draw('axis')
leg = legend(0.18, 0.74, 0.93, 0.88, ncol=2)
for r, h in gs: h.Draw('hist same'); leg.AddEntry(h, 'run %d' % r, 'l')
leg.Draw(); finish(cv); cv.SaveAs(OUT + '/sync_fig_astro.png')

# ---------------------------------------------------------------- LaTeX tables
def f(x, fmt='%d'): return fmt % x if x is not None else '--'
with open(OUT + '/tab_decode.tex', 'w') as t:
    for r in RUNS:
        if r not in dec: continue
        d = dec[r]; k0 = d.get('kcu0', {}); k1 = d.get('kcu1', {}); c = d.get('caen', {}); a = d.get('astropix', {})
        cell = lambda k: ('%s & %.1f & %.2f' % ('{:,}'.format(k['n']).replace(',', '\\,'), 100 * k['complete'], 100 * k['crc_bad'] / max(k['lines'], 1))) if k.get('n') else '-- & -- & --'
        t.write('%d & %s & %s & %s & %s \\\\\n' % (r, cell(k0), cell(k1), '{:,}'.format(c.get('n', 0)).replace(',', '\\,'), '{:,}'.format(a.get('hits', 0)).replace(',', '\\,')))
with open(OUT + '/tab_sync.tex', 'w') as t:
    for r in RUNS:
        if r not in syn: continue
        s = syn[r]
        off = ('%d' % s['off1']) if 'off1' in s else '--'
        sh = ('%.1f/%.1f' % (100 * s['off1_shares'][0], 100 * s['off1_shares'][1])) if 'off1_shares' in s else '--'
        cst = len(set(s['off1_slices'].tolist())) == 1 if 'off1_slices' in s else None
        drift = ('%+.2f' % s['caen_drift_ppm']) if 'caen_drift_ppm' in s else '--'
        mt = ('%.1f' % s['caen_frac_matched'][0]) if 'caen_frac_matched' in s else '--'
        ct = ('%.2f' % s['caen_frac_shifted'][0]) if 'caen_frac_shifted' in s else '--'
        mad = ('%.1f' % s['caen_mad_ns']) if 'caen_mad_ns' in s else '--'
        pd = ('%.3f/%.3f' % (s['pad_same'], s['pad_decoy'])) if 'pad_same' in s else '--'
        t.write('%d & %s%s & %s & %s & %s & %s & %s & %s \\\\\n' % (r, off, '' if cst in (None, True) else '$^*$', sh, drift, mt, ct, mad, pd))
with open(OUT + '/tab_astro.tex', 'w') as t:
    for r in RUNS:
        if r not in syn or 'astro_zmax' not in syn[r]: continue
        s = syn[r]; lat = '--'
        if 'astro_lat_m' in s:
            e = s['astro_lat_edges']; ex = s['astro_lat_m'] - s['astro_lat_x']; sig = ex > 3 * np.sqrt(s['astro_lat_x'] + 1)
            lat = ', '.join('%+d' % (0.5 * (e[i] + e[i + 1])) for i in np.flatnonzero(sig & (ex > 0.05 * ex.max())))
        t.write('%d & %s & %.1f & %+.3f & %.1f & %d & %.1f & %s \\\\\n' % (r, '{:,}'.format(int(s['astro_ncl'])).replace(',', '\\,'), s['astro_wide_z'], s['astro_wide_dt_s'],
                s['astro_zmax'], s['astro_phase'], s['astro_zctl'], lat))
print('wrote figures and tables to', OUT)
