# What `sim/analysis/stack_shares_calib6_v1.pdf` shows, and how it was made (2026-10-23)

Script: `sim/analysis/stack_shares_calib6.py [simdir] [tag] [anchor_GeV]`, run as
`stack_shares_calib6.py /sphenix/tg/tg01/jets/tmengel/eic/sim/v1_purePb _v1 4`. Numbers: `stack_shares_calib6_v1.txt`.
Inputs: `data_fit_layers.npz` (from `extract_fit_layers.py`), `data_sfil5_bbcal.npz`, the v1 simulation files.

## 1. The figure

Horizontal axis: the six elements of the stack in beam order, SFIL1, SFIL2, SFIL3, SFIL4 (KCU0), SFIL5 (KCU1), Baby BCAL
(CAEN). Vertical axis: the share of the visible energy, i.e. each element's visible (fibre) energy divided by the sum of
the six. Energy that escapes the stack, energy in lead, AstroPix and paddles are not part of the denominator, in the
simulation or in the data.

- Coloured lines: Geant4 positrons at 3, 4, 5, 6 GeV (5000 each) in the **v1 pure-lead geometry** (`sim/v1_purePb`),
  beam under paddle P6 (x = +30 mm, y = -6.4 mm). Per element the mean over events of the South-end light,
  `sfil_ls` summed over the layer's cells (5 cells of 12.8 mm in the model) and `bb_ls` summed over the Baby BCAL, divided
  by the sum of the six means. South-end light = Birks-quenched fibre energy times exp(-d_S / 3000 mm), with the fibre
  length assumed 300 mm; because every layer has the same fibre length, the South-light shares equal the visible-energy
  shares in the model to 3 decimals. Each point of a line is one simulation; the lines join them.
- Black points with the black line: run 3062, all four paddles together, the "data" described in section 2.
- Inner (thick) bar: "measurement" systematics plus the anchor-energy range, section 4.
- Outer (thin) bar: the same plus the gain-uniformity term, section 4.
- Grey band: the lowest to highest value among the four single-paddle samples P3, P4, P5, P6 (section 3). It is a
  physics spread (the beam energy differs between paddles), not an error.

## 2. The data points

**SFIL1-4 (KCU0).** From the full-run template fits, `kcu0/fit_full/kcu0_sfil_fit_full.root` (`kcu0_sfil_fit.C`):
2,907,386 triggers kept with: complete 20-sample waveforms on all SFIL links, no CRC or Hamming error, no backward or
overflowing timestamps, exactly one paddle with a dip > 40 ADC (pile-up and no-paddle triggers removed). Per channel the
amplitude is

- the template-fit amplitude when the fit was accepted (pedestal from samples 1-3, raw peak > 60 ADC, peak sample 7-9,
  chi2/ndf < 8 with 5 ADC noise and a 5 % shape error, saturated samples masked, so saturated pulses are extrapolated
  from the template; valid to ~1800 ADC);
- otherwise the raw pedestal-subtracted peak (max of samples 4-12 minus the pedestal) minus that channel's noise-peak
  bias, the mean raw peak in triggers without any pulse (between -1.5 and +2 ADC on the South channels).

The four SFIL1 North channels are masked (noisy board), so the North end cannot be used for SFIL1 and the points use
the **South end only**. Per trigger and layer the four South cells are summed; the shower sample is the 1,971,542
triggers with a South total above 600 ADC (the rest are empty or edge triggers; the 25th percentile of the total is
290 ADC, the median 1840). The all-paddle means are 56.8 / 634.4 / 862.1 / 732.9 ADC for SFIL1-4; by paddle, from the
fit tree's own paddle number (P3 445,888, P4 447,470, P5 579,989, P6 498,195 showers).

**SFIL5 and Baby BCAL.** From `bic_testbeam/data/events_3062.root`, the triple-matched events (KCU0 complete, KCU1 and
CAEN matched): exactly one paddle bit set (P3 = 1, P4 = 2, P5 = 4, P6 = 8; see `docs/pitfalls.md` for the bitmask
mistake this replaced), South peak sum of SFIL1-4 above 600 ADC: 136,584 showers (P3 30,766, P4 30,617, P5 40,894,
P6 34,307). SFIL5 = sum of the raw pedestal-subtracted peaks of the KCU1 South arrays B, C, D (channels 18-34, 36-52,
54-70, clipped at 0), mean 8137 KCU1 ADC. Baby BCAL = sum of the CAEN South charges of channels 0-11 (clipped at 0),
mean 37,443 CAEN units. These are raw peaks and charges, not template fits.

**Conversion to shares.** Three gains, all defined on the P6 showers against the simulation with the beam under P6
at the anchor energy (4 GeV here):

- one common gain for SFIL1-4, g24 = (P6 mean of SFIL2 + SFIL3 + SFIL4) / (simulated South light of SFIL2 + SFIL3 +
  SFIL4) = 13.25 ADC per MeV. SFIL1 is excluded from the definition because its implied gain is 0.4 of the others
  (section 5), and the same number is applied to all four layers: **equal gains assumed**;
- g5 = P6 mean of SFIL5 / simulated SFIL5 South light = 212.5 KCU1 ADC per MeV;
- gbb = P6 mean of the Baby BCAL / simulated Baby BCAL South light = 703.2 CAEN units per MeV.

The all-paddle means divided by these gains give 4.3 / 47.9 / 65.1 / 55.3 / 38.3 / 53.2 MeV (sum 264 MeV; the model
at 4 GeV has 277 MeV, the difference is the SFIL1 deficit and the all-paddle-versus-P6 difference). The plotted shares
are these divided by their sum: 0.016 / 0.181 / 0.247 / 0.209 / 0.145 / 0.202.

By construction the SFIL5 and Baby BCAL points reproduce the model at the anchor energy up to the all-paddle/P6
difference and the SFIL1 deficit: they are a **normalisation, not a measurement**, until KCU1 and CAEN units are tied
to KCU0 ADC by a calibration run. The shape of the SFIL1-4 points is a measurement, with the equal-gain assumption.

## 3. The paddle band

Each paddle's shares use that paddle's own means for all six elements (SFIL1-4 from the fit tree, SFIL5 and Baby BCAL
from the triple-matched file) with the gains of section 2 held fixed. Results: P3 0.018 / 0.185 / 0.261 / 0.205 /
0.141 / 0.190; P4 0.020 / 0.193 / 0.270 / 0.194 / 0.141 / 0.182; P5 0.014 / 0.173 / 0.230 / 0.220 / 0.147 / 0.216;
P6 0.015 / 0.179 / 0.239 / 0.212 / 0.148 / 0.206. The band is the min-max envelope. The ordering P4 < P3 < P6 < P5 in
shower depth is the same in all three readouts and is the energy gradient across the paddles.

## 4. The error bars

All are absolute on the share, evaluated by recomputing the whole chain with one input changed and taking half the
range of the variations (or the full shift when there is one), then added in quadrature:

| source | what was varied | SFIL1 | SFIL2 | SFIL3 | SFIL4 | SFIL5 | Baby BCAL |
| --- | --- | --- | --- | --- | --- | --- | --- |
| shower threshold | South total > 300 and > 1200 instead of 600 ADC | 0.002 | 0.004 | 0.007 | 0.005 | 0.006 | 0.009 |
| saturation recovery | fitted amplitude of saturated channels x 0.9 and x 1.1 (they are 11 / 47 / 38 % of the SFIL2-4 sums) | <0.001 | 0.002 | 0.006 | 0.003 | 0.003 | 0.004 |
| noise-peak bias | +-1 ADC per channel taken from the raw peak (3.3-4 per layer) | 0.001 | <0.001 | <0.001 | <0.001 | 0.001 | 0.001 |
| statistics | means + 1 sigma | 0 | 0 | 0 | 0 | <0.001 | <0.001 |
| anchor energy | the three gains defined at 3, 5, 6 GeV instead of 4 | +-0.001 | -0.013/+0.008 | -0.017/+0.011 | -0.014/+0.009 | -0.009/+0.011 | -0.021/+0.034 |
| **inner bar** | quadrature of the above | 0.002 | -0.013/+0.009 | -0.020/+0.015 | -0.016/+0.011 | -0.011/+0.013 | -0.023/+0.036 |
| gain uniformity | each SFIL layer's gain x exp(+-0.175), one layer at a time, quadrature over layers | 0.003 | 0.028 | 0.035 | 0.031 | 0.009 | 0.013 |
| **outer bar** | inner (+) gain uniformity | | | | | | |

The gain-uniformity number is the only estimate of the equal-gain assumption available without a calibration run:
the North and South ends of the same layers give layer sums in the ratio 0.858 / 0.972 / 0.640 for SFIL2-4 (all
paddles; attenuation is symmetric for the all-paddle sample), and the rms of ln(N/S) over the three layers, 0.175, is
taken as the uncertainty of any one layer's gain relative to the others. The SFIL5 and Baby BCAL entries of that row
are not their own gain uncertainty (their gains are anchored) but the effect of the SFIL1-4 variations on the
denominator. The anchor-energy row is the full change of the SFIL5 and Baby BCAL points when they are anchored to a
different model energy; since the data energy is unknown within 3-6 GeV this is the honest size of those two points'
meaning.

## 5. Caveats and assumptions

1. **Equal gains within SFIL1-4.** The data show it is not exact: the two ends of SFIL4 differ by 36 %, of SFIL2 by
   14 %. Whether that is SiPM/electronics gain or light collection is unknown; route-3 photostatistics (N-S) or the
   LED runs 3179-3190 can tell, the cosmic runs 3088-3172 can replace the assumption outright.
2. **SFIL1.** Its implied gain is 0.4 of the SFIL2-4 gain at every anchor energy (4.5 vs 10-13 ADC/MeV at 5 GeV). The
   plotted SFIL1 point therefore sits at 0.016 against 0.03-0.05 in the model. Either SFIL1 South has a real deficit
   (gain, light yield, or a dead cell) or the model has too much energy in the first layer (upstream material missing in
   the model would move energy the other way). Not resolved.
3. **SFIL5 and Baby BCAL are anchored** to the model (section 2). Their points cannot disagree with the model at the
   anchor; only their paddle-to-paddle variation and the gain-free double ratios carry information.
4. **The model is the v1 pure-lead geometry** (radiator 11.35 g/cm3, 1.9 X0 per layer). The corrected v2 geometry
   (lead/glue 8.56 g/cm3, 1.43 X0 per layer, Baby BCAL per the FTBF paper) is being simulated; it moves the shower
   maximum one layer later, so the model lines will shift and the anchored gains will change. The data ADC shares do
   not change; the SFIL1-4 points move only through the common gain (a renormalisation) and the anchored SFIL5/BBCAL
   points move with the model.
5. **Lateral containment.** Data layer sums cover four cells of unknown height; the model sums five cells of an
   assumed 12.8 mm. The simulated beam at y = -6.4 mm sits on a cell boundary, which is why cell-level comparisons
   (the 2026-10-22 gains) failed and layer sums are used. If the real cells are larger, the data contain more of the
   lateral tail than the model's five cells (SFIL height 64 mm in the model, 16 % of the energy escapes laterally).
6. **South end only, attenuation.** The model attenuates with 3000 mm over fibres of 300 mm (assumed), a 1-2 % effect
   between the P5 and P6 positions; the data have the same fibre length for all four layers so the South-only choice
   cancels in the shares to that level. The paddle band includes whatever real attenuation does.
7. **Beam energy.** Unknown per paddle (3-6 GeV, rising from P4 to P5); the all-paddle sample is a mixture. The model
   lines are mono-energetic.
8. **Two samples.** SFIL1-4 come from 1.97 M KCU0 showers, SFIL5 and the Baby BCAL from the 137 k triple-matched
   showers (a 7 % subsample selected by KCU1 and CAEN availability). The per-element means are assumed to be the same
   in both; the paddle mix differs slightly (P5 is 30 % of the triple sample and 29 % of the full one).
9. **Raw peaks for SFIL5 and the Baby BCAL.** No saturation recovery, no pulse fit, no noise-peak bias correction on the
   51 KCU1 and 12 CAEN channels; clipping at 0 biases small sums up. The anchor absorbs the mean of this but not its
   paddle dependence.
10. **Simulation statistics.** 5000 events per energy; the model shares have statistical errors of about 0.002.
11. **What the model does not contain.** (v1/v2 as used here) the AstroPix module was placed after SFIL1 and the Baby BCAL
    50 mm behind SFIL5, cells were 12.8 mm and fibres 300 mm; v3 (2026-10-23 evening) takes these from the epic BeamTest
    branch and puts the AstroPix module in front of SFIL1. v1/v2 also carried a fibre-to-cell bookkeeping bug (pitfalls.md)
    that affects cell-level numbers only.

## 6. v2 version (`stack_shares_calib6_v2.pdf`, 2026-10-23)

Same procedure against the corrected geometry: gains at the 4 GeV anchor 15.73 ADC/MeV (SFIL1-4 common), 173.0 KCU1
ADC/MeV, 315.5 CAEN units/MeV; shares 0.012 / 0.130 / 0.176 / 0.150 / 0.151 / 0.382 against the v2 model at 4 GeV
0.023 / 0.097 / 0.166 / 0.180 / 0.152 / 0.383. The model lines moved (shower maximum one layer later, Baby BCAL share
0.35-0.43); the SFIL1-4 data points keep their ADC shape, so SFIL2 is now above and SFIL4 below every model line under
equal gains. The photostatistics measurement of method 4.13 (SFIL2 gain 1.49 x SFIL3) removes most of the SFIL2
excess; see `stack_shares_sfil_pe_v2.pdf`.
