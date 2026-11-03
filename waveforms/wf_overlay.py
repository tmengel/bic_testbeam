"""Individual KCU0 waveforms (ADC vs sample) overlaid per channel, from the ts_dump kcutime file.
    python3 wf_overlay.py [out/kcutime_0001.root]
Writes out/kcuhists/wf_overlay_pulses.png (waveforms with a pulse) and wf_overlay_random.png (any waveforms).
"""
import sys, uproot, numpy as np, ROOT
import sys; sys.path.insert(0, '/sphenix/user/tmengel/eic/bic_testbeam/common')
from plotstyle import setup, label
FIN = sys.argv[1] if len(sys.argv) > 1 else 'out/kcutime_0001.root'
K = [2, 21, 63, 43, 6, 25, 59, 37, 10, 29, 52, 0, 14, 33, 48, 4, 146, 165, 207, 187, 150, 169, 203, 186, 154, 173, 196, 144, 158, 177, 192, 148, 84, 86, 81, 97]
NAME = {10: 'S1 cell 3', 29: 'S2 cell 3', 52: 'S3 cell 3', 0: 'S4 cell 3', 154: 'N1 cell 3', 173: 'N2 cell 3', 196: 'N3 cell 3', 144: 'N4 cell 3',
        84: 'P3', 86: 'P4', 81: 'P5', 97: 'P6', 63: 'S3 cell 1', 207: 'N3 cell 1', 59: 'S3 cell 2', 203: 'N3 cell 2'}
SHOW = [10, 29, 52, 0, 154, 173, 196, 144, 63, 207, 59, 203, 84, 86, 81, 97]
t = uproot.open(FIN)['kcutime']
wf = t['wf'].array(library='np', entry_stop=100000)
rng = np.random.default_rng(3062)
setup()
ROOT.gStyle.SetPadTopMargin(0.08); ROOT.gStyle.SetPadLeftMargin(0.16); ROOT.gStyle.SetPadBottomMargin(0.15); ROOT.gStyle.SetPadRightMargin(0.04)
ROOT.gStyle.SetTitleSize(0.06, 'xy'); ROOT.gStyle.SetLabelSize(0.05, 'xy'); ROOT.gStyle.SetTitleOffset(1.1, 'x'); ROOT.gStyle.SetTitleOffset(1.3, 'y'); ROOT.gStyle.SetNdivisions(505, 'xy')
ROOT.gROOT.ForceStyle()
keep = []
for tag, sel in (('pulses', True), ('random', False)):
    c = ROOT.TCanvas('c' + tag, '', 1600, 1560); band = ROOT.TPad('b' + tag, '', 0, 1 - 40 / 1560, 1, 1); band.Draw(); grid = ROOT.TPad('g' + tag, '', 0, 0, 1, 1 - 40 / 1560); grid.Draw(); grid.Divide(4, 4, 0.002, 0.002); keep += [band, grid]
    for i, ch in enumerate(SHOW):
        k = K.index(ch); w = wf[:, k, :]
        ped = np.median(w[:, :4], axis=1)
        if sel: idx = np.flatnonzero((ped - w.min(axis=1) > 40) if ch in (84, 86, 81, 97) else (w.max(axis=1) - ped > 100))
        else: idx = np.arange(len(w))
        idx = rng.choice(idx, min(30, len(idx)), replace=False)
        grid.cd(i + 1); fr = ROOT.gPad.DrawFrame(-0.5, -20, 19.5, 1060, ';sample (25 ns);ADC'); keep.append(fr)
        for j, n in enumerate(idx):
            g = ROOT.TGraph(20, np.arange(20, dtype=float), w[n].astype(float)); g.SetLineColor(ROOT.TColor.GetColorPalette(int(255 * j / 30))); g.SetLineWidth(1); g.SetMarkerStyle(20); g.SetMarkerSize(0.4); g.SetMarkerColor(g.GetLineColor()); g.Draw('LP'); keep.append(g)
        label(f'ch {ch}  {NAME[ch]}', x=0.2, y=0.94, size=0.065)
        label(f'{len(idx)} of {len(w)} waveforms' + ('  (with a pulse)' if sel else ''), x=0.2, y=0.87, size=0.045)
    band.cd(); tx = ROOT.TLatex(); tx.SetNDC(); tx.SetTextFont(42); tx.SetTextSize(0.4); tx.SetTextAlign(12)
    tx.DrawLatex(0.01, 0.5, 'KCU0 run 3062 seg 0001: individual 20-sample waveforms (stock decoder + ROC_pool), ' + ('30 random waveforms with a pulse (> 100 ADC above the pedestal; paddles > 40 below)' if sel else '30 random waveforms, no selection')); keep.append(tx)
    c.SaveAs(f'out/kcuhists/wf_overlay_{tag}.png')
print('wrote out/kcuhists/wf_overlay_{pulses,random}.png')
