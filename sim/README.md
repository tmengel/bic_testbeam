# bic_sim: Geant4 model of the Hall D BIC test-beam setup (August 2026, runs 3060-3064)

Geant4 10.7 (sPHENIX cvmfs build), physics list FTFP_BERT_EMZ, ROOT output. Standalone: `source
../../setup.sh; mkdir build; cd build; cmake ..; make; ./bic_sim run.mac` (batch) or `./bic_sim` (OpenGL).
Energy scan on condor: `condor/scan.job` (plain `condor_submit`), outputs in
`/sphenix/tg/tg01/jets/tmengel/eic/sim/`, summary `analysis/sim_summary.py`.

## Geometry (beam along +z, x along the fibres from North (-x) to South (+x), y vertical, cell 1 at the top)

| Element | Model | Source / status |
| --- | --- | --- |
| Pb/SciFi matrix | radiator = ePIC `SciFiPb_PbGlue` (lead sheets + glue, 96.2 % Pb by mass, 8.558 g/cm3; X0 of the matrix 1.45 cm as in GlueX) with 1 mm fibres (0.96 mm polystyrene core 1.049 g/cm3, 0.04 mm PMMA cladding), fibre pitch 1.34 mm in y, rows 1.22 mm apart in z, rows staggered by half a pitch, 17 rows per layer = 20.74 mm = 1.43 X0, plus a 0.61 mm Pb/glue edge skin (`SciFiPb_PbGlue_Edge`, 4.231 g/cm3) on each face | ePIC `compact/ecal/bic/bic.xml` and `compact/materials.xml` (github.com/eic/epic, main, 2026-10); TIC deck "17 rows of fiber between corrugated lead". **Before 2026-10-23 the radiator was pure lead (11.35 g/cm3), which made the layers 1.9 X0; those outputs are in `sim/v1_purePb/`.** |
| SFIL1-5 | five such layers, 4 readout cells stacked in y, SFIL1 front face at z = 0; 17 mm slots between layers with 2 mm carbon-fibre frame plates on both sides (the ePIC imaging-layer frame) | layer count and cells from the hookup; slot and frames from ePIC `bic.xml`; **fibre length 300 mm and cell height 12.8 mm are ASSUMED** |
| AstroPix | one 9-chip module, chip pitch 20 mm, 35 x 35 pixels of 0.5 mm (17.5 mm active), stack Si 500 um (sensitive) + Si 150 um + Cu 100 um + Kapton 200 um + epoxy 100 um + 0.5 mm carbon stave, chip 0 at the North end | ePIC stave stack; **placed in the slot after SFIL1 (`astroSlot 1`): ASSUMED** |
| Baby BCAL | GlueX-type module as measured at Fermilab (arXiv:2504.03079): 22.2 cm of Pb/ScFi (184 lead sheets; modelled as 10 depth segments of 18 fibre rows = 21.96 cm, 15.1 X0), fibres 58 cm long, trapezoid 8.4-11.8 cm wide modelled as a 10.1 cm box of 4 cells, 0.8 cm aluminium plate upstream and 3.175 cm downstream, 40 light guides per end read as 10 depth rows x 4 columns (the paper sums each column into 4 channels of 1, 2, 3, 4 arrays); 50 mm behind SFIL5 | the paper; **distance ASSUMED**; the same paper gives SF 9.5 %, Moliere radius 3.63 cm, X0 1.45 cm, double-clad Kuraray SCSF-78MJ fibres, Birks 0.132 mm/MeV |
| Paddles | four 20 x 60 x 10 mm plastic scintillators at z = -200 mm placed under the AstroPix chips they cover (P4 chip 0, P3 chips 2-3, P6 chips 5-6, P5 chip 8) | positions from the run-3062 analysis; **sizes ASSUMED** |

With the corrected radiator the 4 GeV shower peaks in SFIL4 (layer visible energies 8 / 34 / 58 / 61 / 53 MeV), the SFIL fibre sampling fraction is 0.112, the stack holds 78 % and 13 % escapes; the full scan with these materials is `condor/scan.job` (both beam positions).

Every number is a `/bic/geom/` command (see `macros/geometry_default.mac` and `include/BicGeometry.hh`);
set them before `/run/initialize`. Replace the ASSUMED values with the setup metrology (TIC slide 37,
J. Richards' simulation) when available. The real SFIL cells are probably larger than 12.8 mm: with the
default height (4 x 12.8 = 51 mm) a 4 GeV shower leaks about 16% of its energy sideways and out of the back.

Local copies of the public sources (ePIC `bic.xml`/`materials.xml`, the Baby BCAL and AstroPix FTBF papers, the PDR2 recap and
the 2026 EICUG talk) are in `../../context/`, indexed with URLs in `../../context/SOURCES.md`.

## Readout model

- Fibre energy is recorded as deposited (`efib`) and Birks-quenched (`evis`, kB = 0.126 mm/MeV).
- North / South light: `evis x exp(-d/lambda)` with d the distance to the fibre end, lambda = 3000 mm
  (`attenuationLength`). Position: `sfil_x` is the evis-weighted x of the deposits in the cell.
- Timing: `sfil_tn`, `sfil_ts` are the earliest arrival times at the two ends, deposit time +
  distance / v_eff (170 mm/ns); -1 when the cell is empty. No SiPM, shaping or TDC model: this is the
  physics limit for the North-South method of the analysis plan.
- AstroPix: energy per pixel above `pixelThreshold` (1 keV default; a MIP leaves ~190 keV in 500 um),
  with the first deposit time.
- `e_escape` is the kinetic energy leaving the world, `e_other` deposits outside the listed volumes.

## Output tree `sim` (one entry per beam particle; MeV, ns, mm)

`ebeam x0 y0 pdg`; SFIL per layer (index 0-4) and cell (0-3, top to bottom): `sfil_efib sfil_evis
sfil_epb sfil_ln sfil_ls sfil_x sfil_tn sfil_ts` (arrays `[8][8]`, unused slots 0), `sfil_e[8]`
(fibre + lead per layer), `sfil_etot`; Baby BCAL `bb_efib bb_evis bb_epb bb_ln bb_ls` (`[16][8]`),
`bb_e[16]`, `bb_etot`; `paddle_e[4]` (P3 P4 P5 P6); `astro_n astro_chip astro_row astro_col astro_e
astro_t`; `e_escape e_other e_dep`. The `info` TNamed holds the geometry used.

First numbers (100 positrons of 4 GeV on the beam row, default geometry): SFIL fibre sampling fraction
0.083 (fibre / fibre+lead), 5 SFILs hold 65% of the energy, the Baby BCAL 17%, 16% escapes; visible
fibre energy 274 MeV with sigma/mu = 0.094; North-South arrival-time difference rms 0.10 ns for a 5 mm
spot; 8 AstroPix pixels per event above threshold after 1.4 X0.

## Energy scan (condor, 5000 e+ per point) and comparison with run 3062

`analysis/sim_summary.png` / `sim_summary.txt`: sigma/mu of the visible fibre energy = 6.0%/sqrt(E) (+) 2.1%,
72.5 MeV visible per GeV, SFIL fibre sampling fraction 0.089, 80% of the energy in SFILs + Baby BCAL and
16% escaping (lateral, from the assumed 51 mm SFIL height). The fourth panel shows the shares of the whole stack
(SFIL1-5 and the Baby BCAL) for the 3-6 GeV simulation and for run 3062 (`analysis/data_shares.npz`:
triple-matched showering triggers of `bic_testbeam/data/events_3062.root`; SFIL1-4 = sqrt(S x N) per cell summed
over cells 1, 3, 4; SFIL5 = KCU1 arrays B+C+D peak sum; BBCAL = CAEN charge of ch 0-11). The KCU1 and CAEN scales
are fixed so that the all-paddle SFIL5 and BBCAL means relative to SFIL1-4 equal the 4 GeV simulation, so those two
points are a normalisation and only their paddle-to-paddle variation (P4 North end, P5 South end shown) is a
measurement; the relative gain calibration of the three readouts replaces this. The SFIL1-4 part of the data: the data peak one layer later
than the 3-6 GeV simulation and carry more energy in SFIL4, as expected if the real cells are larger (less
lateral leakage) or there is more material upstream; the data are uncalibrated ADC with 20-36% of all beam-row waveforms saturated (53-66% of shower waveforms)
on the beam row, so the comparison is qualitative until step 3 of the analysis plan.

## Paddle P5 comparison and the South-end gain calibration (2026-10-22)

**Correction 2026-10-23: the P5/P6 selection of this section was wrong.** `pad` in `events_3062.root` is a bitmask
(P3 = 1, P4 = 2, P5 = 4, P6 = 8); `pad == 5` / `pad == 6` selected two-paddle events (1084 / 958 of them). With the
single-paddle masks (33468 P5, 27856 P6 showers; `south_gain_calibration.py [simdir]`, old output kept as
`south_gain_calibration_WRONGPAD_2026-10-22.txt`) the gain-free double ratios are 1.135 (SFIL5) and 1.201 (Baby
BCAL), above even the scaled May map (1.081, 1.103): P5 showers are deeper, P5 is at higher energy than P6, and
the SFIL5 / SFIL1-4 and Baby BCAL / SFIL1-4 South ratios rise monotonically from P4 (North) through P3 and P6 to P5
(3.37, 3.53, 3.73, 3.91 for SFIL5). The "P5 = P6 within 0.3 GeV" conclusion below is withdrawn; `data_shares.npz`
per-paddle entries and `data_p5_south.npz` were regenerated with the correct masks and `stack_shares.py`,
`stack_shares_p5.py` rerun.


`analysis/stack_shares_p5.png` (`stack_shares_p5.py`): P5-tagged showers (South end of the module) with South-end
amplitudes only, against the simulation with the beam under P5 (x = +80 mm, `condor/scan.job` with the beam
position argument) using the South-end light; 5 and 6 GeV curves. `analysis/south_gain_calibration.{py,txt,png}`:
per-element gains data/sim for the P5 and P6 samples under four energy hypotheses, and the gain-free double ratios
[element / SFIL3 cell 3]_P5 / [...]_P6. Result: the double ratios of SFIL5 and the Baby BCAL (1.006, 1.039 in data)
sit on the equal-energy prediction (1.010, 0.988), not on the scaled May map (1.081, 1.103), and the P5/P6 gain
ratios are 0.98 +- 0.04 under equal energies against 0.84 under the map: P5 and P6 see nearly the same energy
(within ~0.3 GeV), so the August paddles do not follow the scaled May calibration. Preliminary gains at 5 GeV,
ADC per MeV of simulated South-end light: SFIL1-4 beam row 20 / 33 / 23.5 / 22.8, SFIL5 South 176 (KCU1 ADC,
sum of arrays B+C+D), Baby BCAL South 504 (CAEN charge units); 10-15% uncertain from the absolute energy, and
tied to the assumed 12.8 mm cell size (the simulation puts more light in the cells next to the beam row than the
data show, pointing to larger real cells).

## Improved data points and preliminary per-layer gains (2026-10-23)

`analysis/extract_fit_layers.py` -> `data_fit_layers.npz`: per trigger of the full-run template fits
(`kcu0/fit_full/kcu0_sfil_fit_full.root`) the South-end amplitude of each SFIL layer summed over its four cells, the
fitted amplitude (saturation recovered) where a fit was accepted and the raw pedestal-subtracted peak minus the channel's
noise-peak bias (-1.5 to +2 ADC, measured in pulse-free triggers) elsewhere; showers = South total > 600 ADC
(1.97 M of 2.91 M fitted triggers). `analysis/stack_shares_calib.py [simdir] [tag]` -> `stack_shares_sfil_calib<tag>.pdf`,
`sfil_gains<tag>.pdf`, `sfil_gain_calibration<tag>.txt`. `analysis/stack_shares_sfil.py` is the same share figure
with the old points (3-6 GeV simulation, SFIL1-4 only, ticks on the layer labels), used for the presentation PDF.

Against the pure-Pb (v1) simulation, South light summed over the cells of a layer:

| SFIL1-4 shares | SFIL1 | SFIL2 | SFIL3 | SFIL4 |
| --- | --- | --- | --- | --- |
| data, previous (peak amplitude sqrt(SxN), beam-row cells) | 0.063 | 0.312 | 0.403 | 0.222 |
| data, fitted amplitude, South end, all cells, all paddles | 0.025 | 0.277 | 0.377 | 0.321 |
| data P4 / P6 / P5 | 0.030 / 0.023 / 0.022 | 0.286 / 0.278 / 0.272 | 0.398 / 0.371 / 0.361 | 0.286 / 0.328 / 0.345 |
| sim v1 3 / 4 / 5 / 6 GeV | 0.070 / 0.061 / 0.054 / 0.049 | 0.267 / 0.251 / 0.239 / 0.230 | 0.362 / 0.366 / 0.367 / 0.368 | 0.301 / 0.323 / 0.340 / 0.354 |

The old SFIL4 deficit was the measurement, not the detector: the beam-row cell alone with saturated peaks and the
North/South geometric mean lost a third of SFIL4. With the fitted, all-cell South sums SFIL2-4 fall between the 3 and
5 GeV curves; SFIL1 stays at less than half of the simulated share.

Per-layer gain = <data South amplitude> / <simulated South light> (ADC per MeV, P6 showers against the beam under P6;
P5 against x = 80 mm gives the same within 5-10 %):

| beam energy hypothesis | SFIL1 | SFIL2 | SFIL3 | SFIL4 | SFIL2-4 mean (rms) |
| --- | --- | --- | --- | --- | --- |
| 3 GeV | 5.3 | 16.9 | 16.6 | 17.7 | 17.1 (3 %) |
| 4 GeV | 4.8 | 14.1 | 12.9 | 13.0 | 13.3 (4 %) |
| 5 GeV | 4.5 | 12.3 | 10.7 | 10.2 | 11.1 (8 %) |
| 6 GeV | 4.3 | 11.1 | 9.2 | 8.5 | 9.6 (11 %) |

These replace the beam-row-cell gains quoted above (20 / 33 / 23.5 / 22.8): those compared one data cell with one
simulated cell while the simulated beam at y = -6.4 mm straddles two 12.8 mm cells, so they were about twice too high.
SFIL2-4 need the same gain within 3-4 % if the energy is 3-4 GeV; SFIL1 needs 0.4 of it under every hypothesis
(either a real gain or light-yield deficit of SFIL1 South, or the model has too much energy in the first layer).
Applying per-layer gains anchored to the simulation reproduces the simulated shares by construction, so the data
points are drawn with equal gains and the gains are reported separately. Equal-gain energy estimate from the SFIL2-4
profile (chi2 of the renormalised shares against the interpolated model, 2 % per share): all paddles 3.2 GeV
(3.0-3.9), P3 3.0, P4 3.0 (at the scan edge), P6 3.4 (3.0-4.0), P5 4.0 (3.3-5.0) with the v1 geometry; the v2 geometry,
which moves the shower maximum later, will push these lower. The P4 -> P6 -> P5 ordering is the May-map direction, but the
P5-P6 step is within the gain-free double-ratio result (P5 = P6 within ~0.3 GeV).

## Six-element shares with gains and systematic uncertainties (2026-10-23)

`analysis/stack_shares_calib6.py [simdir] [tag] [anchor_GeV]` -> `stack_shares_calib6<tag>.pdf`, `.txt`
(`data_sfil5_bbcal.npz`: SFIL5 South arrays B+C+D and Baby BCAL South means of the 136.6k triple-matched single-paddle
showers, per paddle). SFIL1-4 use one common gain (equal gains assumed; 13.25 ADC/MeV at a 4 GeV anchor from the P6
layer sums 2-4), SFIL5 and the Baby BCAL their own gains anchored to the simulation with P6 showers against the beam
under P6 (212 KCU1 ADC/MeV, 703 CAEN units/MeV at 4 GeV): those two points sit on the model at the anchor by
construction and carry the 3-6 GeV anchor range as their bar. Systematics on the SFIL1-4 shares (absolute, v1 model,
all paddles): shower threshold 300/1200 instead of 600 ADC <= 0.007, saturation recovery (+-10 % on the fitted
amplitude of saturated channels, 11 / 47 / 38 % of the SFIL2-4 sums) <= 0.006, noise-peak bias +-1 ADC per raw-peak
channel <= 0.001, statistics < 0.001; the equal-gain assumption, estimated from the North/South spread of the same
layers (N/S = 0.86 / 0.97 / 0.64, rms of ln = 0.175 taken as the per-layer gain uncertainty), 0.028 / 0.035 / 0.031
on SFIL2-4 and dominant; the paddle-to-paddle spread (physics, grey band) up to 0.02. Shares all paddles at the 4 GeV
anchor: 0.016 / 0.181 / 0.247 / 0.209 / 0.145 / 0.202 against the 4 GeV model 0.040 / 0.164 / 0.239 / 0.211 / 0.145 /
0.201.

## v2 scan results and the photostatistics cross-check (2026-10-23)

v2 scan (`sim/bic_sim_e{3000,4000,5000,6000}_s2*.root`, `e{5000,6000}_x80`): `sim_summary.py` -> sigma/mu
6.2 %/sqrt(E) (+) 2.5 %, 86.6 MeV visible per GeV, SFIL sampling fraction 0.112, SFIL + Baby BCAL 77-78 %, escape
13-14 %, shower maximum in SFIL4. The v1 figures are kept as `sim_summary_v1.*`, `south_gain_calibration_v1.*`; every
share/gain script takes the simulation directory and a tag (`_v1`, `_v2`). `stack_shares_calib.py` on v2: implied layer
gains at 4 GeV 7.2 / 20.5 / 15.9 / 13.0 ADC/MeV (SFIL2-4 rms/mean 19 %); the equal-gain profile fit pins to the 3 GeV
edge with chi2 14.5, i.e. equal gains and v2 are incompatible. `stack_shares_calib6.py` on v2 (anchor 4 GeV): common
SFIL1-4 gain 15.7, SFIL5 173 KCU1 ADC/MeV, Baby BCAL 316 CAEN units/MeV.

`photostat_ns.py` (method.md 4.13): from the North/South fluctuations of run 3062, ADC per photoelectron of the South
beam-row channels SFIL2 14.09 +- 0.16 and SFIL3 9.47 +- 0.06, ratio 1.49; the v2 model implies 1.29-1.41, v1 1.09-1.20.
Both ends of SFIL2 and SFIL3 see equal photoelectrons; SFIL4 North is anomalous (channels 186/187 constant ~130 ADC).
`stack_shares_pe.py` -> `stack_shares_sfil_pe_v{1,2}.pdf`: shares with the SFIL2/SFIL3 ratio applied (SFIL1 and SFIL4
assumed at the SFIL3 gain): 0.027 / 0.205 / 0.415 / 0.353; v2 at 4 GeV 0.049 / 0.208 / 0.356 / 0.387. The profile
energy fit is degenerate with the unmeasured SFIL4 gain (v2 -> 3 GeV edge, v1 -> 6 GeV edge) and is not an energy
measurement.

SFIL5 South and Baby BCAL: no data-only photoelectron calibration possible from 3062 (intra-array pairs uncorrelated, no dark-count peak; method.md 4.13). `stack_shares_calib6.py <simdir> _pe_v2 4 pe` -> `stack_shares_calib6_pe_v2.pdf`: SFIL1-4 with the photostatistics ratio, SFIL5 / Baby BCAL anchored.

## Geometry v3 from the epic BeamTest branch (2026-10-23)

The collaboration's DD4hep description of the August stack (`eic/epic`, branch `BeamTest`,
`compact/ecal/baby_barrel_interlayers_Summer2026.xml`, J. Richards; table in `../../context/SOURCES.md`) replaces the
remaining assumptions:

| parameter | v2 (assumed) | v3 (official) |
| --- | --- | --- |
| SFIL readout cells | 5 x 12.8 mm (64 mm) | 4 x 19.05 mm (76.2 mm slab) |
| fibre length | 300 mm | 510 mm (505-519 per layer) |
| SFIL centre pitch | 37.0 mm (17 mm slot) | 65.1 mm (SFIL4-5 61.9 mm); air gap 43.1 mm |
| SFIL5 back to Baby BCAL front | 50 mm | 61 mm |
| fibre spacing, cladding | 1.34 mm, 0.04 mm | 1.35 mm, 0.06 mm |
| paddles | 20 x 60 x 10 mm at -200 mm | 23.8 x 137 x 14 mm slabs, P5/P3 ~180 mm and P6/P4 ~60 mm upstream of SFIL1 |
| beam row | y = -6.4 mm (a cell boundary) | y = -9.5 mm (centre of cell 3) |
| AstroPix module | slot after SFIL1 | **kept**: one 9-chip AstroPix v3 module was in the August stack (its data are in `astropix_3062.root`); the official DD4hep file has the imaging detector commented out, which is a gap in that description, not a reason to drop it |

Not changed: the Pb/glue radiator, 17 rows per layer, the Baby BCAL (the official file has 6 x 17 + 4 x 20 rows = 182,
the model 10 x 18 = 180), the readout model. Carbon frames of the ePIC imaging slot are still placed at the slot edges;
in the test stand the 43 mm gaps are air. v2 outputs are archived in `sim/v2_slots17mm/`; the v3 scan writes to `sim/`.
Official paddle x positions relative to SFIL3 (model frame, North negative): P4 -43, P3 +46, P6 +55, P5 +144 mm; the
scan keeps x = 30 and 80 mm for continuity (attenuation over 510 mm fibres makes the difference below 1 %).

**Cell-assignment bug fixed with v3 (2026-10-23):** the fibre copy number overflowed its cell digit for rows 10-16, booking
41 % of each layer's fibres one cell too low in every v1 and v2 output. Cell-level numbers from those files (off-row light,
per-cell gains, the "cells larger than 12.8 mm" inference, cell-level double ratios) are void; layer sums are unaffected.
The beam now sits at y = -5 mm, 5 mm inside cell 3 on the cell-2 side, as the official slab offset (`relativeY = -0.5 cm`)
and the data pattern (cell 3 dominant, some light in cell 2, none in cell 4) indicate. The AstroPix module is in front of
SFIL1 (`astroSlot 0`): the matched data clusters are 1-2 pixels, a single track, where a module behind 1.4 X0 would give 5-6.

**v3 scan results (2026-10-23 night, `sim_summary.py`, `*_v3` figures):** SF 0.102, 76 MeV visible per GeV, sigma/mu 4.6 % (3 GeV)
to 3.8 % (6 GeV), SFIL + Baby BCAL 74-76 %, escape 15-16 %, shower maximum SFIL4, AstroPix 0.7 pixels above 25 keV per event
(in front of SFIL1; the data clusters average 1.36 pixels). Six-element shares at 4 GeV 0.024 / 0.100 / 0.169 / 0.183 /
0.151 / 0.372 (v2: 0.023 / 0.097 / 0.166 / 0.180 / 0.152 / 0.383). Implied layer gains at 4 GeV 7.8 / 22.8 / 18.0 / 14.8
ADC/MeV; SFIL2/SFIL3 1.27 (4 GeV) to 1.40 (6 GeV) against the photostatistics 1.49. SFIL3 cell shares 0.03 / 0.24 / 0.65 /
0.08 against the data 0.28 (cells 1+2) / 0.72 / 0.005. Gain-free double-ratio predictions for the scaled May map: SFIL5
1.062, Baby BCAL 1.095 (data 1.135, 1.201).

**Unsaturated sample (2026-10-23 night):** `stack_shares_pe.py <simdir> <tag> nosat|nosat_raw` keeps the 30 % of showers without any
saturated channel and applies the equivalent per-cell threshold (930 ADC / implied gain) to the model (80 % kept: the data
saturate much more often than the model predicts for the mean gains). Raw peaks and fits agree within 0.004 there. Data
0.037 / 0.200 / 0.391 / 0.373; model with the same cut 0.050 / 0.197 / 0.353 / 0.400 (4 GeV), and 0.054 / 0.205 / 0.354 / 0.387
when cell 4, empty in the data, is dropped from the model as well. Residual: SFIL3 +0.035, SFIL1 -0.02.

**Saturation handled by quantile matching, saturation-edge gains, SFIL5/SFIL4 (2026-10-24):** `stack_shares_pe.py <simdir> <tag> nosat`
now drops from the model, per layer, the same fraction of its highest-light events as the data lose to the 10-bit ceiling (equivalent
selection without a gain); `sat_gain.py` compares the gain from the mean amplitude with the gain from the saturation fraction (ratio
1.3, energy-independent: the data amplitude distribution is broader than the model's); `sfil5_sfil4_gain.py` joins the SFIL4 fits to
the KCU1 SFIL5 South sums through the KCU0 timestamp and gives g5/g4 = 13.4-15.7 KCU1 ADC per KCU0 ADC (3-6 GeV). Relative layer
gains from the amplitude estimators: SFIL2 : SFIL3 : SFIL4 = 1.35 : 1 : 0.85 (+-0.07); photostatistics SFIL2/SFIL3 = 1.49.
`draw_event.py` renders one simulated shower (steps tree, `/bic/stepsEvents N`) over the geometry: `event_display_4GeV.png`.

**SFIL1-5 shares with SFIL5 on the KCU0 scale (2026-10-24, `stack_shares_sfil5.py [simdir] [tag] [anchor]` -> `stack_shares_sfil5_pe_v3.pdf`):**
480 k showers whose trigger also has KCU1; SFIL5 = KCU1 arrays B+C+D divided by g5/g4 = 14.5 (the data ratio over the model's
L5/L4 at 4 GeV, so <SFIL5>/<SFIL4> equals the model ratio at the anchor by construction; the bar spans the 3-6 GeV range, +-8 %).
Shares 0.021 / 0.160 / 0.324 / 0.275 / 0.220 against the model at 4 GeV 0.042 / 0.169 / 0.274 / 0.286 / 0.229: SFIL2, SFIL4 and
SFIL5 within 0.01-0.02, SFIL3 +0.05, SFIL1 -0.02. Per paddle P4 -> P5 the SFIL4 and SFIL5 shares rise (0.252/0.202 -> 0.291/0.233).

## Geometry v4: the August survey (2026-10-29)

`context/MeasureMap2_Complete.pdf` (3 pages: x-z positions, AstroPix/paddle relative positions, dimensions) is the as-built survey
of the setup. It is the same survey the epic BeamTest file encodes (`SFILn_X/Z/Length/Thickness`, `BBCal_offsetX`, the paddle
`relativeX/Z`), so v3 was already close to it in z; v3 had averaged it (one fibre length, equal pitch, everything centred). v4
(`/bic/geom/layout 1`, the default; `layout 0` = v3) takes every number from the survey (`BicGeometry.hh`, `sv*` constants).

Survey conventions: x along the fibres, measured from the beam-pipe centre to each element's edge nearest the pipe; that is the South
end (the paddle z planes and the run-3062 chip map put P4 = chip 0 = North on the far side, P5 nearest the pipe, as the epic file
labels them). z distances are measured upstream from the Baby BCAL east face (the back of its 31.75 mm Al support) to each
element's downstream face. Model frame: x = 0 at the AstroPix module centre, North negative; z = 0 at the SFIL1 upstream face.

| Element | v3 | v4 (survey) |
| --- | --- | --- |
| SFIL1-5 fibre length | 510 mm, all | 519.1 / 504.8 / 517.5 / 504.8 / 504.8 mm |
| SFIL1-5 x centre | 0 | -7.6 / -10.0 / -27.5 / -16.4 / -29.1 mm (South ends at +251.9 / +242.4 / +231.3 / +236.1 / +223.4) |
| SFIL1-5 upstream face z | 0 / 65.1 / 130.1 / 195.2 / 260.2 mm | 0 / 64.3 / 130.2 / 195.3 / 258.0 mm (pitch 64.3 / 65.9 / 65.1 / 62.7) |
| carbon frames in the SFIL slots | 2 x 2 mm per slot | none (not in the survey or the epic test-beam file) |
| AstroPix | x centre 0, centre z = -21.6 mm | x -90 .. +90 mm (South edge 41.91 cm from the pipe), downstream face z = -6.5 mm |
| Baby BCAL | front 61 mm behind SFIL5, x centre 0 | east face at the survey benchmark: front z = 331.1 mm (51.1 mm behind SFIL5), x centre -33.3 mm |
| paddles x (P3 P4 P5 P6) | -30 / -80 / +80 / +30 mm (data chip map) | +3.5 / -74.3 / +100.3 / +13.0 mm |
| paddles z centre | -180 / -60 / -180 / -60 mm | -166.0 / -66.0 / -183.5 / -64.4 mm |
| paddle size x / y / z | 23.8 / 137 / 14 mm | 23.8 / 149.2 / 14.3 mm |

The readout uses each layer's own centre and length for the North / South light paths (and the Baby BCAL its x centre).
`/geometry/test/run` finds no overlaps; a 20-event 4 GeV run gives the v3 shower profile. The `info` TNamed records the layout.

Kept from v3, and why:
- **Active height 4 x 19.05 mm.** The survey heights H1-H5 (11.43 / 11.27 / 11.27 / 11.43 / 10.48 cm) are the slab envelopes; the
  epic file, built from the same survey, keeps the Pb/ScFi at 76.2 mm. If the 4.5 in were active, the cells would be 28.6 mm.
- **Vertical positions.** The survey gives the AstroPix top 3.97 cm below the SFIL1 envelope top (Y1), but not where the active
  76.2 mm sits inside the envelope, so y stays from the data (beam at y = -5 mm, AstroPix centred).
- **Layer thickness.** W1-W5 (22.2-23.8 mm) are envelopes; the 17-row layer + edge skins is 21.96 mm and is centred in each.
- **Beam x** stays in the AstroPix frame from the run-3062 hits (P6 triggers at x = +30 mm, P5 at +80 mm).

**Open conflict between the survey and the run-3062 data.** The data put each paddle over its own AstroPix chips with P(paddle) =
0.96-1.0 there and ~0.05 elsewhere: P4 chip 0, P3 chips 2-3, P6 chips 5-6, P5 chip 8 (method 4.5), i.e. P3 and P6 ~60 mm apart.
The survey (and the epic file) put P3 and P6 only 9.5 mm apart, both near chip 4-5, and P5 at +100 mm, beyond chip 8 (which ends
at +90 mm). With the survey paddles, a beam at the data's P6 position (+30 mm) crosses P6 in only 5 of 20 test events. Either the
paddle x numbers 18 / 19 (49.37 / 48.42 cm) or the AstroPix x reference (21 = 41.91 cm; page 2 gives 42.23 cm via X1) are off by a
few cm, or the module's chip order along x is not 0-8. Nothing in the sim analysis selects on `paddle_e`, so this affects only
the upstream material the beam sees, but it should be settled with whoever made the survey before the paddle positions are used.

**v4 scan results (2026-10-30, `sim/v4/`, seeds 81-88, `sim_summary.py tag=_v4` -> `sim_summary_v4.*`, `stack_shares_pe.py sim/v4
_v4`, `stack_shares_calib.py sim/v4 _v4` -> `sfil_gain_calibration_v4.txt`, `stack_shares_sfil_pe_v4.*`):** SF 0.1022 (unchanged),
stack containment 0.77 (v3 0.75) with the Baby BCAL visible share 6-7 % higher and the SFIL share 1 % lower (the survey z and the
51 mm gap move light to the back), sigma/mu 4.2 / 3.9 / 3.6 / 3.4 % at 3-6 GeV. SFIL1-4 shares at 4 GeV 0.047 / 0.203 / 0.356 /
0.393 (v3 0.051 / 0.210 / 0.355 / 0.384). Run 3062 against v4: pe-equalised SFIL2-4 shares give E = 3.0 GeV (3.0-3.3) for all
paddles, P5 3.3 (3.0-4.0), as with v3; implied South gains at the 4 GeV / P6 anchor 8.6 / 24.2 / 18.4 / 14.7 ADC per MeV of South
light (v3 7.8 / 22.8 / 18.0 / 14.8), at 3 GeV 9.2 / 28.1 / 22.9 / 19.6; the scan now also has P5 (x = 80 mm) points at 3 and 4 GeV.

**v4b = v4 + the measured light model (2026-10-30, `sim/v4b/`, seeds 91-96; `EXTRA='/bic/geom/attenuationAlpha 0;/bic/geom/attenuationLength 1500 mm;/bic/geom/vEff 160'`):**
single-exponential attenuation 1500 mm (the SFIL2/3 AstroPix-x measurement) and v_eff 160 mm/ns instead of the FTBF two-component form and
170 mm/ns. Shares and shape energy identical to v4; implied 4 GeV / P6 gains 9.0 / 25.5 / 19.2 / 15.4 ADC per MeV of South light (+4 %);
tN-tS rms 0.18 ns. Outputs `sim_summary_v4b.*`, `stack_shares_sfil_pe_v4b.*`, `sfil_gain_calibration_v4b.txt`.

## Cosmic mode (2026-10-30)

`/bic/gun/cosmic 1` (macro `macros/cosmic.mac`, condor `condor/cosmic.{sh,job}`, output `/sphenix/tg/tg01/jets/tmengel/eic/sim/cosmic/`):
particles start on a horizontal plane at `cosmicY` (250 mm, above the paddles), uniform in x and z over `cosmicX0..X1` x
`cosmicZ0..Z1` (-400..400 x -350..750 mm, paddles to Baby BCAL), zenith angle from the sea-level cos^2 distribution weighted by
the flux through a horizontal plane (cos theta = u^1/4), uniform azimuth, kinetic energy log-flat between `cosmicEmin` and
`cosmicEmax` (0.5-20 GeV); `/bic/gun/particle mu-`. The tree gains `z0, ux, uy, uz` (vertex and direction). A 1000-muon test: 2-3 %
of the generated muons cross a given SFIL layer, a SFIL3 cell-3 crossing leaves 1.56 MeV visible (IQR 1.2-2.2 MeV) and ~1.3 MeV
of light at each end. Purpose: the two-end muon tag of the cosmic run 3088 (`calibration/cosmics/`) applied to the simulated light
per cell gives the path-length-weighted spectrum that turns the measured sqrt(A_S A_N) peaks into ADC per MeV, a data-only check
of the beam-anchored gains (18 ADC/MeV for SFIL3 predicts ~23 ADC per MIP).

**Cosmic target mode (2026-10-30 evening, `/bic/gun/cosmicTarget 1` with `cosmic 1`; `condor/cosmic_tgt.job` -> `sim/cosmic/tgt/`, 100 x 100k,
seeds 301-400).** Why: in run 3088 the SFIL3 two-end muon tags come 97-100 % with a paddle pulse (52 % of all triggers have one, and the tag
rate is 45 x higher in those triggers), spread over all four paddles, at ~1.4 tags per minute for the whole layer. The paddles are vertical
slabs upstream of the stack at beam height, so a muon through a paddle AND a SFIL cell is near-horizontal (zenith ~72 deg, 10-90 % 64-79),
and such a track starts 0.5-1 m upstream of any practical start plane: the plane productions (`sim/cosmic/`, `sim/cosmic/pad/`) hold
2-12 such crossings per cell in 2 M muons. The target mode draws the direction as before (cos^2 zenith, uniform azimuth, log-flat energy),
picks one of the four paddles as placed by the detector (`DetectorConstruction::paddleCentre/paddleHalf`), the entry face with
probability proportional to |d_i| x face area and the entry point uniformly on it (flux-correct for tracks through the paddle), and starts
the track `cosmicBack` (800 mm) before the entry point or at the plane y = cosmicY, whichever is closer. Test (20k): 94 % hit a paddle, 0.16 %
cross a paddle and a SFIL3 cell (edge cells 1 and 4 twice as often as 2 and 3); those tracks also cross SFIL1 (94 %) and SFIL2 (97 %), SFIL4 in
27 %; South light per SFIL3 crossing 1.08 MeV (plane production 1.13), i.e. the path through a cell is about the same along z as along y.
The run-3088 trigger rate (1/min for near-horizontal tracks through 33 cm^2 of paddle face into the 0.6 sr of SFIL3) matches this picture; a
paddle laid flat on the stack would give ~30 tags per minute.
