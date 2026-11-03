"""The LED-macro figures (PlotEventBICSingle.C) for the run-3062 KCU0 data, from kcu_hists output.

4 x 4 grids per channel group (SFIL South, SFIL North, paddles), one figure per quantity:
  *_ADC.png        ADC vs sample                 *_TOA.png   TOA code vs sample
  *_TOT.png        decompressed TOT vs sample    *_ADCcorr.png  waveform time-aligned by the first TOA word
and per ASIC the channel x value maps (ASIC_<q>.png). Panels are labelled "ch <n> <detector>".

    python3 kcu_hists_plot.py [out/kcuhists/kcu0_hists.root]
"""
import sys, ROOT
import sys; sys.path.insert(0, '/sphenix/user/tmengel/eic/bic_testbeam/common')
from plotstyle import setup, finish, label

FIN = sys.argv[1] if len(sys.argv) > 1 else 'out/kcuhists/kcu0_hists.root'
OUT = 'out/kcuhists/'
f = ROOT.TFile.Open(FIN)
nwf = int(f.Get('hNwf').GetBinContent(1))

# channel groups, in the kcuadc order: k = 4*(cell-1) + (layer-1)
S = [2, 21, 63, 43, 6, 25, 59, 37, 10, 29, 52, 0, 14, 33, 48, 4]
N = [146, 165, 207, 187, 150, 169, 203, 186, 154, 173, 196, 144, 158, 177, 192, 148]
P = [84, 86, 81, 97]
def det(ch):
    if ch in S: k = S.index(ch); return f'S{k%4+1} cell {k//4+1}'
    if ch in N: k = N.index(ch); return f'N{k%4+1} cell {k//4+1}'
    if ch in P: return 'P' + str([3, 4, 5, 6][P.index(ch)])
    return ''
GROUPS = [('SFIL_South', S, 4, 4), ('SFIL_North', N, 4, 4), ('paddles', P, 2, 2)]
QUANT = [('ADC', 'hADC_%d', 'sample', 'ADC', None, True),
         ('TOA', 'hTOA_%d', 'sample', 'TOA code', None, True),
         ('TOT', 'hTOT_%d', 'sample', 'TOT (decompressed)', None, True),
         ('ADCcorr', 'hADCcorr_%d', 'time from the first TOA word [ns]', 'ADC', (-50, 450), True),
         ('ADCpulse', 'hADCpulse_%d', 'time from the pulse TOA word [ns]', 'ADC', (-50, 450), True)]

setup()
ROOT.gStyle.SetPadTopMargin(0.08); ROOT.gStyle.SetPadLeftMargin(0.16); ROOT.gStyle.SetPadBottomMargin(0.15); ROOT.gStyle.SetPadRightMargin(0.04)
ROOT.gStyle.SetTitleSize(0.06, 'xy'); ROOT.gStyle.SetLabelSize(0.05, 'xy'); ROOT.gStyle.SetTitleOffset(1.1, 'x'); ROOT.gStyle.SetTitleOffset(1.3, 'y'); ROOT.gStyle.SetNdivisions(505, 'xy')
ROOT.gROOT.ForceStyle()
keep = []
for gname, chans, nx, ny in GROUPS:
    for qname, hname, xt, yt, xr, logz in QUANT:
        c = ROOT.TCanvas('c' + gname + qname, '', 400 * nx, 380 * ny + 40)
        band = ROOT.TPad('band', '', 0, 1 - 40.0 / (380 * ny + 40), 1, 1); band.Draw(); keep.append(band)
        grid = ROOT.TPad('grid', '', 0, 0, 1, 1 - 40.0 / (380 * ny + 40)); grid.Draw(); keep.append(grid)
        grid.Divide(nx, ny, 0.002, 0.002)
        for i, ch in enumerate(chans):
            grid.cd(i + 1); ROOT.gPad.SetLogz(logz)
            h = f.Get(hname % ch); keep.append(h)
            h.SetTitle(''); h.GetXaxis().SetTitle(xt); h.GetYaxis().SetTitle(yt)
            if xr: h.GetXaxis().SetRangeUser(*xr)
            if qname == 'TOT': h.GetYaxis().SetRangeUser(0, 4096)
            h.Draw('col')
            label(f'ch {ch}  {det(ch)}', x=0.2, y=0.94, size=0.065)
        band.cd(); t = ROOT.TLatex(); t.SetNDC(); t.SetTextFont(42); t.SetTextSize(0.4 if nx == 4 else 0.26); t.SetTextAlign(12)
        note = {'ADCcorr': ';  aligned to the FIRST TOA word: t = 25 s - (25 s_w + TOA_ns), forward convention', 'ADCpulse': ';  aligned to the pulse TOA word (sample after the largest ADC step): t = 25 s - (25 s_w + TOA_ns) + 200',
                'TOT': ';  decompressed 12-bit TOT, fires on saturated pulses', 'TOA': ';  raw 10-bit code', 'ADC': ''}[qname]
        t.DrawLatex(0.01, 0.5, f'KCU0 run 3062 seg 0001, {nwf} complete waveforms, {gname.replace("_", " ")}: {yt} vs {xt}{note}'); keep.append(t)
        c.SaveAs(OUT + f'{gname}_{qname}.png')

# per ASIC: channel x value
for q in range(4):
    c = ROOT.TCanvas(f'cA{q}', '', 1500, 450); c.Divide(3, 1, 0.004, 0.004)
    for i, (hname, yt) in enumerate((('hADCASIC_%d', 'ADC'), ('hTOAASIC_%d', 'TOA code'), ('hTOTASIC_%d', 'TOT (decompressed)'))):
        c.cd(i + 1); ROOT.gPad.SetLogz(True); ROOT.gPad.SetRightMargin(0.13)
        h = f.Get(hname % q); keep.append(h); h.SetTitle(''); h.GetXaxis().SetTitle(f'ASIC {q} channel (KCU0 ch {72*q}-{72*q+71})'); h.GetYaxis().SetTitle(yt)
        if 'TOT' in hname: h.GetYaxis().SetRangeUser(0, 4096)
        h.Draw('colz')
        label({0: 'ASIC0: SFIL1-4 South', 1: 'ASIC1: paddles P3-P6 (ch 84, 86, 81, 97)', 2: 'ASIC2: SFIL1-4 North', 3: 'ASIC3: Baby BCAL North (silent)'}[q], x=0.17, y=0.94, size=0.045)
    c.SaveAs(OUT + f'ASIC_{q}.png')
print('wrote', OUT, 'figures for', nwf, 'waveforms')
