#!/bin/bash
# Rebuild figs/ for the November 3 deck from the analysis sources: trim the white margins, add an 8 px WHITE border
# (never ImageMagick's default grey), crop the multi-panel figures to the panels shown. Run from docs/slides_nov03.
set -e; B=../..; F=figs; mkdir -p $F
tb() { convert "$1" -fuzz 2% -trim +repage -bordercolor white -border 8 "$2"; }      # trim + white border
plain="calibration/stability/gain_stability_kcu0.png sim/analysis/photostat_ns.png calibration/cells/cell_equalisation.png sim/analysis/bbcal_gains_yfine.png sim/analysis/sfil5_sfil4_gain_rec.png sim/analysis/stack_shares_sfil_pe_v4.png sim/analysis/event_display_4GeV.png sim/analysis/sim_summary_v4.png sim/analysis/performance_3064.png astropix/astro_imaging.png waveforms/kcu0_sfil3_pedsub.png docs/beam_rate.png timing/toa_clean/fig_xres.png astropix/astro_tot_calib.png docs/kcu0_channel_map.png dual_sided/fig_asym_pos.png dual_sided/fig_dt_pos.png timing/toa_clean/fig_t0.png"
for f in $plain; do tb $B/$f $F/$(basename $f); done
# crops: first trim the source, then cut, then trim + border
T=$(mktemp -d)
convert $B/calibration/attenuation/astro_x_attenuation.png -fuzz 2% -trim +repage $T/a.png; W=$(identify -format %w $T/a.png); H=$(identify -format %h $T/a.png); convert $T/a.png -crop ${W}x$((H/2))+0+0 +repage $T/a2.png; tb $T/a2.png $F/attenuation_top.png
convert $B/calibration/cosmics/cosmic_mip_gm_3088.png -fuzz 2% -trim +repage $T/g.png; W=$(identify -format %w $T/g.png); H=$(identify -format %h $T/g.png); convert $T/g.png -crop $((W*75/100))x$((H/4))+$((W*25/100))+$((H/2)) +repage $T/g2.png; tb $T/g2.png $F/cosmic_gm_3panels.png
convert $B/calibration/cosmics/cosmic_mip_sim_3088tgt.png -fuzz 2% -trim +repage $T/c.png; W=$(identify -format %w $T/c.png); H=$(identify -format %h $T/c.png); convert $T/c.png -crop $((W*60/100))x$((H/4))+$((W*40/100))+$((H/2)) +repage $T/c2.png; tb $T/c2.png $F/cosmic_sfil3_3panels.png
convert $B/sim/analysis/performance.png -fuzz 2% -trim +repage $T/p.png; W=$(identify -format %w $T/p.png); H=$(identify -format %h $T/p.png); convert $T/p.png -crop ${W}x$((H/2))+0+0 +repage $T/p2.png; tb $T/p2.png $F/performance_top.png
rm -rf $T
# the template-fit figure is made here from the fit file
source /sphenix/user/tmengel/eic/setup.sh >/dev/null 2>&1
python3 - <<'PY'
import ROOT, sys
sys.path.insert(0, '/sphenix/user/tmengel/eic/bic_testbeam/common'); from plotstyle import setup, legend, finish
setup(); ROOT.gROOT.SetBatch(True); ROOT.gStyle.SetPadRightMargin(0.04)
f = ROOT.TFile.Open('/sphenix/tg/tg01/jets/tmengel/eic/kcu0/fit_full/kcu0_sfil_fit_full.root'); keep = []
c = ROOT.TCanvas('c', '', 1800, 760); c.Divide(2, 1, 0.002, 0.002)
c.cd(1); h = f.Get('hTemplate_18'); h.SetTitle('SFIL3 South beam cell: pulse template;t - t_{cross} [ns];normalised amplitude'); h.SetLineWidth(3); h.SetLineColor(ROOT.kBlue + 1); h.SetStats(0); h.GetXaxis().SetRangeUser(-60, 320); h.Draw('hist'); keep.append(h)
lg = legend(0.45, 0.72, 0.92, 0.86); lg.AddEntry(h, 'mean of clean unsaturated pulses', 'l'); lg.Draw(); keep.append(lg)
c.cd(2); ROOT.gPad.SetRightMargin(0.14); ROOT.gPad.SetLogz(); h2 = f.Get('hAfitVsAmp_18'); h2.SetTitle('SFIL3 South beam cell: fitted amplitude vs raw peak;raw peak - pedestal [ADC];template-fit amplitude [ADC]'); h2.SetStats(0); h2.Draw('colz'); keep.append(h2)
l = ROOT.TLine(930, 0, 930, h2.GetYaxis().GetXmax()); l.SetLineStyle(2); l.SetLineWidth(2); l.Draw(); keep.append(l)
t = ROOT.TLatex(); t.SetNDC(); t.SetTextSize(0.04); t.DrawLatex(0.18, 0.82, '10-bit ceiling: the fit'); t.DrawLatex(0.18, 0.77, 'extrapolates through it'); keep.append(t)
finish(c); c.SaveAs('figs/tplfit_sfil3_raw.png')
PY
tb $F/tplfit_sfil3_raw.png $F/tplfit_sfil3.png; rm $F/tplfit_sfil3_raw.png
for f in $F/*.png; do printf "%-28s %s edge %s\n" $(basename $f) $(identify -format '%wx%h' $f) "$(convert $f -format '%[pixel:p{0,20}]' info:)"; done
