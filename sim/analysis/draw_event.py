"""Event display of one simulated shower in the test-beam geometry (PyROOT): side view (beam z horizontal, vertical y) and top view
(z, x along the fibres). Charged steps as polylines (e- blue, e+ red, other charged grey), photons as thin green lines, over the
outlines of the paddles, the AstroPix module, SFIL1-5 with their four cells and the Baby BCAL with its aluminium plates.
    python3 draw_event.py <bic_sim file with a 'steps' tree> [event index] [out.png]
"""
import sys, numpy as np, uproot, ROOT
sys.path.insert(0, '/sphenix/user/tmengel/eic/bic_testbeam/common'); from plotstyle import setup
fn = sys.argv[1]; iev = int(sys.argv[2]) if len(sys.argv) > 2 else 0; out = sys.argv[3] if len(sys.argv) > 3 else 'event_display.png'
f = uproot.open(fn); st = f['steps'].arrays(library='np'); ev = f['sim'].arrays(['ebeam', 'sfil_evis', 'bb_evis', 'e_escape'], library='np')
x, y, z, e, pdg, trk = (st[k][iev].astype(float) if k != 'pdg' and k != 'trk' else st[k][iev] for k in ('x', 'y', 'z', 'edep', 'pdg', 'trk')); E = ev['ebeam'][iev]
pitch, tmod, hh = 20.74 + 2 * 0.61 + 43.1, 20.74 + 2 * 0.61, 38.1; bb0 = 343.2; bbfront, bbdepth, bbback = 8., 10 * 18 * 1.22, 31.75; bbh = 2 * 25.25
# paddles as the MODEL places them (DetectorConstruction.cc): over AstroPix chips 0, 2-3, 5-6, 8 (P4, P3, P6, P5) at x = (chip - 4) x 20 mm,
# z = paddleZpos + p x (paddleZ + 2 mm); the data-driven paddle map (channel_map.md, PADX) supports these x positions, the official
# DD4hep file's paddle X values do not (they put P3 and P6 within 1 cm of each other, which the two-paddle rates exclude)
PX = [(0 - 4) * 20., (2.5 - 4) * 20., (5.5 - 4) * 20., (8 - 4) * 20.]; PZ = [-180. + 120. + 15., -180., -180. + 120. + 15., -180. + 15.]; PW, PH, PT = 23.8, 137., 14.   # P4, P3, P6, P5: two z planes (official), x from the paddle map
setup(); ROOT.gStyle.SetPadLeftMargin(0.07); ROOT.gStyle.SetPadRightMargin(0.02); ROOT.gStyle.SetPadTopMargin(0.12); ROOT.gStyle.SetPadBottomMargin(0.16); ROOT.gROOT.ForceStyle()
c = ROOT.TCanvas('c', '', 1800, 1000); c.Divide(1, 2, 0.001, 0.001); keep = []
def box(x1, y1, x2, y2, fill, line=ROOT.kGray + 2, alpha=1.0):
    b = ROOT.TBox(x1, y1, x2, y2); b.SetFillColorAlpha(fill, alpha); b.SetLineColor(line); b.SetLineWidth(1); b.Draw('l'); keep.append(b)
for ipad, (vert, lab, half, is_side) in enumerate(((y, 'y, vertical [mm]', hh, True), (x, 'x along fibres, N (-) to S (+) [mm]', 255., False))):
    c.cd(ipad + 1); fr = ROOT.gPad.DrawFrame(-230, -half * 1.2, 640, half * 1.2, ';z along the beam [mm]  (SFIL1 front face at 0);' + lab); fr.GetXaxis().SetTitleSize(0.07); fr.GetYaxis().SetTitleSize(0.07); fr.GetXaxis().SetLabelSize(0.06); fr.GetYaxis().SetLabelSize(0.06); fr.GetYaxis().SetTitleOffset(0.42); fr.GetYaxis().CenterTitle(); fr.GetXaxis().SetTitleOffset(1.0); fr.GetYaxis().SetNdivisions(505); keep.append(fr)
    for L in range(5):
        box(L * pitch, -half, L * pitch + tmod, half, ROOT.kGray)
        if is_side:
            for cc in range(1, 4): ln = ROOT.TLine(L * pitch, hh - cc * 19.05, L * pitch + tmod, hh - cc * 19.05); ln.SetLineColor(ROOT.kGray + 2); ln.SetLineStyle(3); ln.Draw(); keep.append(ln)
    box(-22.2, -17.5 if is_side else -90, -21.0, 17.5 if is_side else 90, ROOT.kOrange - 3, ROOT.kOrange + 3)
    box(bb0, -bbh if is_side else -290, bb0 + bbfront, bbh if is_side else 290, ROOT.kGray + 1); box(bb0 + bbfront, -bbh if is_side else -290, bb0 + bbfront + bbdepth, bbh if is_side else 290, ROOT.kGray); box(bb0 + bbfront + bbdepth, -bbh if is_side else -290, bb0 + bbfront + bbdepth + bbback, bbh if is_side else 290, ROOT.kGray + 1)
    for p in range(4): box(PZ[p], -PH / 2 if is_side else PX[p] - PW / 2, PZ[p] + PT, PH / 2 if is_side else PX[p] + PW / 2, ROOT.kYellow - 9, ROOT.kOrange + 3)
    order = np.lexsort((np.arange(len(trk)), trk)); tz, tv, tp, tt = z[order], vert[order], pdg[order], trk[order]
    inside = (tz > -230) & (tz < 640) & (np.abs(tv) < half * 1.2)
    for sel, col, width, alpha in ((tp == 22, ROOT.kGreen + 2, 1, 0.18), (~np.isin(tp, (22, 11, -11)), ROOT.kGray + 3, 1, 0.6), (tp == 11, ROOT.kAzure + 2, 1, 0.7), (tp == -11, ROOT.kRed + 1, 2, 0.9)):
        idx = np.where(sel)[0]
        if len(idx) == 0: continue
        for seg0 in np.split(idx, np.where(np.diff(tt[idx]) != 0)[0] + 1):
          for seg in np.split(seg0, np.where(np.diff(inside[seg0].astype(int)) != 0)[0] + 1):
            if len(seg) < 2 or not inside[seg[0]]: continue
            pl = ROOT.TPolyLine(len(seg), tz[seg], tv[seg]); pl.SetLineColorAlpha(col, alpha); pl.SetLineWidth(width); pl.Draw('l'); keep.append(pl)
    if is_side:
        t = ROOT.TLatex(); t.SetTextFont(42); t.SetTextSize(0.055); t.SetTextAlign(21)
        for z0, name in [(-120, 'paddles'), (-45, 'AstroPix'), (bb0 + 130, 'Baby BCAL')] + [(L * pitch + 11, 'SFIL%d' % (L + 1)) for L in range(5)]: t.DrawLatex(z0, hh * 1.03, name)
        keep.append(t)
    ROOT.gPad.RedrawAxis()
c.cd(1); sf = ev['sfil_evis'][iev][:5, :4].sum(axis=1); bb = ev['bb_evis'][iev].sum()
tt = ROOT.TLatex(); tt.SetNDC(); tt.SetTextFont(42); tt.SetTextSize(0.055); tt.DrawLatex(0.07, 0.93, 'bic_sim, %.0f GeV e^{+} (event %d): visible energy SFIL1-5 = %s MeV, Baby BCAL %.0f MeV, escaping %.0f MeV' % (E / 1000., iev, ' / '.join('%.0f' % v for v in sf), bb, ev['e_escape'][iev])); keep.append(tt)
c.cd(2); t2 = ROOT.TLatex(); t2.SetNDC(); t2.SetTextFont(42); t2.SetTextSize(0.055); t2.DrawLatex(0.07, 0.93, 'top view;  e^{-} blue, e^{+} red, photons green, other charged grey; paddles yellow, AstroPix orange, Pb/ScFi grey, Al dark grey'); keep.append(t2)
c.SaveAs(out); c.SaveAs(out.replace('.png', '.pdf')); print('wrote', out)
