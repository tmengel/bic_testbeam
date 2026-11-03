"""KCU0 channel map drawn along the test-beam setup -> docs/kcu0_channel_map.png (ROOT primitives).
Facts: docs/channel_map.md, maps/map_run3062.json (8/13 hookup table, verified in data).
    python3 channel_map_fig.py
"""
import ROOT
ROOT.gROOT.SetBatch(True); ROOT.gErrorIgnoreLevel = ROOT.kWarning
W, H = 1500, 800
c = ROOT.TCanvas('map', '', W, H); c.Range(0, 0, W, H)   # pixel-like coordinates; Y() flips so y runs downwards as in the layout
keep = []
def Y(y): return H - y
def col(hexs):
    return ROOT.TColor.GetColor(hexs)
def rect(x, y, w, h, fill='#ffffff', stroke='#333333', lw=1):
    b = ROOT.TBox(x, Y(y + h), x + w, Y(y)); b.SetFillColor(col(fill)); b.SetLineColor(col(stroke)); b.SetLineWidth(lw); b.Draw('l'); keep.append(b)
    bl = ROOT.TBox(x, Y(y + h), x + w, Y(y)); bl.SetFillStyle(0); bl.SetLineColor(col(stroke)); bl.SetLineWidth(lw); bl.Draw('l'); keep.append(bl)
def text(x, y, t, size=12, align=22, bold=False, color='#111111'):
    l = ROOT.TLatex(x, Y(y), t); l.SetTextFont(62 if bold else 42); l.SetTextSize(size * 1.3 / H); l.SetTextAlign(align); l.SetTextColor(col(color)); l.Draw(); keep.append(l)
def line(x1, y1, x2, y2, color='#333333', lw=1):
    l = ROOT.TLine(x1, Y(y1), x2, Y(y2)); l.SetLineColor(col(color)); l.SetLineWidth(lw); l.Draw(); keep.append(l)

rect(0, 0, W, H, '#ffffff', '#ffffff')
text(W / 2, 26, 'Run 3062: KCU0 channels along the BIC test-beam setup (beam from the left; North = top, South = bottom of each column)', 16, bold=True)
text(W / 2, 50, 'KCU0 ASIC0 = SFIL1-4 South (ch 0-71), ASIC1 = trigger paddles (ch 72-143), ASIC2 = SFIL1-4 North (ch 144-215), ASIC3 = Baby BCAL North (ch 216-287, silent).', 10.5, color='#444444')
text(W / 2, 66, 'Mapping from the 8/13 shift table (maps/channel_mapping.txt), verified in the data; full table in docs/channel_map.md.', 10.5, color='#444444')
a = ROOT.TArrow(40, Y(350), 225, Y(350), 0.02, '|>'); a.SetLineColor(col('#c0392b')); a.SetFillColor(col('#c0392b')); a.SetLineWidth(4); a.Draw(); keep.append(a)
text(130, 330, 'beam e^{+} (3-6 GeV)', 13, bold=True, color='#c0392b')
# paddles
x = 260; text(x + 55, 112, 'Trigger paddles', 14, bold=True); text(x + 55, 131, 'R580 PMTs, NIM OR', 9.5, color='#444444')
for i, (p, ch, caen, chips) in enumerate((('P4', '86 (J4)', 'CAEN 13', 'chip 0'), ('P3', '84 (J18)', 'CAEN 12', 'chips 2-3'), ('P6', '97 (J6)', 'CAEN 15', 'chips 5-6'), ('P5', '81 (J17)', 'CAEN 14', 'chip 8'))):
    y = 160 + i * 95; rect(x, y, 110, 80, '#fdebd0', '#b9770e', 2)
    text(x + 55, y + 20, p, 15, bold=True); text(x + 55, y + 42, 'A1 ch ' + ch, 12); text(x + 55, y + 58, caen, 10.5, color='#444444'); text(x + 55, y + 72, 'AstroPix ' + chips, 9, color='#666666')
text(x + 55, 556, 'neg. pulses, rail at 0', 9.5, color='#666666'); text(x + 55, 570, 'other A1 ch unused', 9.5, color='#666666')
# AstroPix
x = 425; rect(x, 150, 70, 390, '#eaf2f8', '#2e86c1', 2); text(x + 35, 112, 'AstroPix', 14, bold=True); text(x + 35, 131, '9-chip v3', 9.5, color='#444444')
for i in range(9):
    y = 160 + i * 41; rect(x + 10, y, 50, 34, '#ffffff', '#2e86c1', 1); text(x + 35, y + 21, f'chip {i}', 11)
text(x + 35, 556, 'chip 0 = North', 9.5, color='#666666')
# SFIL1-4
sfil = {1: ((2, 146), (6, 150), (10, 154), (14, 158)), 2: ((21, 165), (25, 169), (29, 173), (33, 177)), 3: ((63, 207), (59, 203), (52, 196), (48, 192)), 4: ((43, 187), (37, 186), (0, 144), (4, 148))}
x0 = 560
for L in range(1, 5):
    x = x0 + (L - 1) * 170; text(x + 70, 112, f'SFIL{L}', 14, bold=True); text(x + 70, 131, 'Pb/SciFi', 9.5, color='#444444')
    rect(x, 150, 140, 390, '#f4f6f7', '#566573', 2)
    for cell in range(1, 5):
        y = 160 + (cell - 1) * 95; beam = cell == 3
        rect(x + 8, y, 124, 84, '#f9e79f' if beam else '#ffffff', '#b7950b' if beam else '#7f8c8d', 2 if beam else 1)
        text(x + 70, y + 16, f'cell {cell}' + (' (beam)' if beam else ''), 11, bold=beam)
        s, n = sfil[L][cell - 1]
        text(x + 70, y + 38, f'N: A2 ch {n}', 12, color='#1a5276'); text(x + 70, y + 57, f'S: A0 ch {s}', 12, color='#7b241c')
        note = 'SiPM1+2 sum (v2)' if cell in (1, 2) and L != 1 else ('SiPM3' if cell == 3 else ('SiPM4' if cell == 4 else 'SiPM1 / SiPM2 (v1)'))
        text(x + 70, y + 74, note, 9, color='#666666')
    if L == 1: text(x + 70, 570, 'North column noisy', 9.5, color='#666666')
    if L == 4: text(x + 70, 570, 'N/S gain offset', 9.5, color='#666666')
text(x0 + 325, 556, 'every cell: North end -> KCU0 ASIC2, South end -> KCU0 ASIC0; beam row = cell 3 (SiPM3)', 9.5, color='#666666')
# SFIL5 and Baby BCAL
x = 1250; rect(x, 150, 110, 390, '#fef9e7', '#9a7d0a', 2); text(x + 55, 112, 'SFIL5', 14, bold=True); text(x + 55, 131, 'single SiPMs, KCU1', 9.5, color='#444444')
for i, (name, chs, st) in enumerate((('array C', '180-196', 'weak'), ('array D', '198-214', 'live'), ('array B', '162-178', 'live'), ('array A', '144-160', 'mostly dead'))):
    y = 160 + i * 95; rect(x + 8, y, 94, 84, '#ffffff', '#9a7d0a', 1); text(x + 55, y + 20, name, 12, bold=True); text(x + 55, y + 42, 'KCU1 A2', 10.5); text(x + 55, y + 57, 'ch ' + chs, 10.5); text(x + 55, y + 73, st, 9.5, color='#666666')
text(x + 55, 556, 'S end only (N dead)', 9.5, color='#666666'); text(x + 55, 570, 'order: A<B<D<C', 9.5, color='#666666')
x = 1385; rect(x, 150, 100, 390, '#f2f3f4', '#626567', 2); text(x + 50, 112, 'Baby BCAL', 14, bold=True); text(x + 50, 131, 'downstream', 9.5, color='#444444')
text(x + 50, 250, 'North:', 12, bold=True); text(x + 50, 270, 'KCU0 A3', 11); text(x + 50, 287, 'ch 216-287', 11); text(x + 50, 303, 'silent', 9.5, color='#666666'); text(x + 50, 317, '252-287 railed', 9.5, color='#666666')
text(x + 50, 400, 'South:', 12, bold=True); text(x + 50, 420, 'CAEN ch 0-11', 11)
line(260, 598, 1485, 598, '#555555', 1); text(870, 588, 'z along the beam (schematic: spacing and the AstroPix position are not to scale)', 10, color='#555555')
y = 630; text(60, y, 'Reading the diagram', 14, align=12, bold=True)
for i, t in enumerate(['Every SFIL cell is read at both fibre ends: North end -> KCU0 ASIC2 (ch 144-215), South end -> KCU0 ASIC0 (ch 0-71); 16 + 16 live channels.',
                       'The beam sits in cell 3 (SiPM3) of every layer. On v2 boards cells 1 and 2 read the same SiPM1+SiPM2 sum (treat them as one row); North SFIL1 is a v1 board.',
                       'The paddle OR triggers KCU0, KCU1 and the CAEN (no common busy). Paddles sit along the AstroPix module, North to South: P4 (chip 0), P3 (chips 2-3), P6 (chips 5-6), P5 (chip 8).',
                       'ASIC1 ch 84/86/81/97 = P3/P4/P5/P6 (8/13 hookup; the pre-8/13 numbers 82/93 are quiet). ASIC3 = Baby BCAL North BBS01-40 on J1-J40: no signal in this run (use the CAEN South end).',
                       'Pedestals 79-102 ADC. SiPM pulses are positive, peak at sample 7 and saturate at 1023; paddle pulses are negative and rail at 0 around samples 8-16.']):
    text(60, y + 26 + i * 22, t, 10.5, align=12, color='#222222')
c.SaveAs('../docs/kcu0_channel_map.png')
print('wrote docs/kcu0_channel_map.png')
