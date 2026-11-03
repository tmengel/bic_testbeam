"""BIC test beam analysis: calibration, simulation and performance update (November 3, 2026).
Builds the deck in the user's own PowerPoint style from context/bic_beamtest.pptx (white, Arial, title top-left, date / page
footer, CU logo on the title slide) and a pixel preview of every slide (preview/slideNN.png, 1920 x 1080) from the same layout,
since no PowerPoint renderer is available on the node.
    python3 build_deck.py  -> bic_beamtest_nov03.pptx, preview/slideNN.png
Figures: figs/ (trimmed copies of the analysis figures, see figs/README.txt). Details, cuts and caveats go to the speaker notes.
"""
import os, copy, re
from pptx import Presentation
from pptx.util import Inches, Pt, Emu
from pptx.dml.color import RGBColor
from pptx.enum.text import PP_ALIGN, MSO_ANCHOR
from pptx.oxml.ns import qn
from PIL import Image, ImageDraw, ImageFont

TEMPLATE = '/sphenix/user/tmengel/eic/context/bic_beamtest.pptx'
DATE = 'November 3, 2026'; OUT = 'bic_beamtest_nov03.pptx'
SW, SH = 13.3333, 7.5                       # inches
TAGCOL = {'DONE': '2B6B47', 'NEW': '2B6B47', 'IN PROGRESS': 'A9541C', 'NEXT': 'A9541C', 'TODO': 'A9541C', 'OPEN': '97243A', 'QUESTION': '97243A', 'NOT POSSIBLE': '97243A'}
BODY_PT, SUB_PT, TITLE_PT, FOOT_PT, TABLE_PT = 17, 15, 40, 14, 13
FONT_R = '/usr/share/fonts/urw-base35/NimbusSans-Regular.otf'; FONT_B = '/usr/share/fonts/urw-base35/NimbusSans-Bold.otf'
PX = 144.                                   # preview pixels per inch (1920 x 1080)

# ---------------------------------------------------------------- content
# bullet = (level, text, tag or None); **bold** segments allowed in text
S = []
def slide(title, bullets=(), fig=None, layout='left', notes='', table=None, figh=None):
    S.append(dict(title=title, bullets=list(bullets), fig=fig, layout=layout, notes=notes, table=table, figh=figh))

S.append(dict(title='BIC test beam analysis', layout='title', bullets=[], fig=None, notes='', table=None,
              subtitle='Calibration, simulation and performance update\nTanner Mengel\n' + DATE))

slide('Next steps from last time', [
    (0, '**1. Energy response and resolution**', 'DONE'), (1, 'template fits recover the saturated pulses; gains anchored to the model; linearity and resolution per paddle', None),
    (0, '**2. Shower shape with the matched AstroPix events**', 'DONE'), (1, '2.8 mm beam height, single-pixel clusters; attenuation and North / South gains from the AstroPix position', None),
    (0, '**3. Synchronization of runs 3063 and 3064**', 'DONE'), (1, 'all August runs on /gpfs01 are decoded, fitted and aligned (3064: 12 h, 18 x run 3062)', None),
    (0, '**4. SFIL5 South channel map**', 'OPEN'), (1, 'array order in height measured from the data; the SiPM assignment still needs the pin map', None),
    (0, '**Paddle to PS energy mapping?**  not available; the data now give their own energy estimate (slide 21)', 'OPEN'),
    (0, '**A-STEP latencies (19 / 88 ms), on-chip 10 MHz?**  no answer yet', 'OPEN'),
], layout='text', notes='Status of the four next steps and three questions of the October 20 talk.')

slide('What is new since October 20', [
    (0, '**Reconstruction and tools:** decoder fixes (PRs opened), raw format and per-trigger ROOT trees, trigger / DAQ loss budget, waveform anatomy, TOA / TOT timing, AstroPix ToT, template fits on every run, all new runs synchronized', None),
    (0, '**Calibration:** gain stability over 12 h, photostatistics, attenuation, cell structure, cosmic MIP scale, Baby BCAL and SFIL5 gains', None),
    (0, '**Simulation:** Geant4 model of the full stack, now at the August survey geometry; 3-6 GeV scans; cosmic mode', None),
    (0, '**Performance:** energy linearity, resolution and containment per paddle for runs 3062 and 3064; imaging layer', None),
    (0, 'Everything is in github.com/tmengel/bic_testbeam (docs/status.md is the overview); large outputs on SDCC', None),
], layout='text', notes='One line per thread; each has its own slides below.')

# ---- reconstruction and tools
slide('Decoders: three fixes, validated against the raw bytes', [
    (0, 'The H2GCROC format is understood bit by bit (docs/data_format.md); the stock decoder kept 6 of the 10 TOT bits, put out-of-order lines into the wrong sample and checked no CRC', None),
    (0, 'All three fixed in forks of online_distribution and ROC_pool, validated bit-exact against the raw bytes (75.5 M channel words per board, 0 mismatches); **pull requests opened October 29**', 'DONE'),
], table=[['item', 'stock decoder', 'now'],
          ['TOT', '6 of 10 bits (capped at 63)', 'full 10-bit TOT; scales with the saturated samples (Spearman 0.74 / 0.89)'],
          ['out-of-order lines', 'start a new, fragmented waveform', 'placed by timestamp in a reorder window'],
          ['CRC-32 per frame', 'not checked', 'checked: KCU0 0 failures on the live board; KCU1 failures only on the dead SFIL5 North board'],
          ['complete waveforms per packet', '61 % KCU0, 48 % KCU1', '75 % KCU0, 79 % KCU1 (88 / 83 % of recorded triggers across packets)']],
    layout='table', notes='decoder_tests/README.md. The central install (M. Purschke, Oct 21) has the TOT mask, the kept last sample and the Hamming flags; the CRC check and the reorder window are in the PRs. '
          'Hamming flags on 2.5 % of KCU0 triggers, mostly the paddle half and the SFIL North half. The remaining incomplete waveforms are UDP datagram loss in the DAQ, not decoding.')

slide('Raw format and per-trigger ROOT trees by subsystem', [
    (0, 'Begin-run packets: 100 us trigger busy, 20 samples per trigger; the ASIC analog-gain block is not recorded, so the KCU0 / KCU1 preamp gains cannot be read from the data', None),
    (0, 'Every stream converted to a ROOT tree, entry i = trigger i of the run, segment boundaries de-duplicated; the same converters ran on all new runs (3064 in parts); common/bic_streams.py reads the four streams together', 'DONE'),
], table=[['stream', 'run 3062 triggers', 'complete', 'size'],
          ['KCU0 (SFIL1-4, paddles)', '3.77 M', '88 %', '19 GB'],
          ['KCU1 (SFIL5 South)', '1.53 M', '83 %', '6 GB'],
          ['CAEN (Baby BCAL South, paddles)', '1.51 M', '-', '7 GB'],
          ['AstroPix', '2.86 M hit records', '-', '60 MB']],
    layout='table', notes='decoding/kcu_to_root.cc, kcu_merge.cc, condor jobs; checks in decoding/kcu0_root_check.py; files on /sphenix/tg/tg01/jets/tmengel/eic/{kcu0,kcu1,caen,astropix}/. '
          'Raw-format write-up: docs/data_format.md (48-word HGCROC line map, characterization mode words Tc Tp ADC TOT TOA, CAEN rollover word, AstroPix record layout). Ask for the DAQ YAML to get Global_Analog.')

slide('Where the triggers go: the DAQ loss budget', [
    (0, 'One trigger = 20 samples x 8 links on one timestamp; TRIGGER_OUT counts them: **4.10 kHz issued on both boards, 2.8 kHz recorded** (union of the three streams)', None),
    (0, 'KCU0 loses ~23 % to mid-packet UDP drops (about one socket buffer each); KCU1 is silent 80 % of the time; CAEN is dead-time limited at 1.0 kHz; 57 % of the issued triggers end as a complete KCU0 waveform', None),
    (0, 'Offline recovers ~95 % of what was received; the rest is DAQ: socket buffer / recvmmsg, common busy, dropouts, CAEN record length (docs/completeness_plan.md)', 'OPEN'),
], fig='figs/beam_rate.png', layout='full',
    notes='decoder_tests/kcu_event_anatomy.cc, docs/method.md 4.11. EVENT_COUNTER is a per-link sample counter and unusable; TRIGGER_IN always 0. 99.85 % of complete KCU0 waveforms carry a paddle pulse (one paddle in 88 %). '
          'Segment 0 starts with 7.14 s of unrecorded triggers until KCU1 started.')

slide('KCU0 waveforms: what the detector does', [
    (0, 'Decoding is bit-exact; the shape is the detector: rise at samples 5-6, peak at 7-8, tail still 55 % of the peak at sample 19', None),
    (0, '**Per-event baseline shift of 17-35 ADC shared by all SFIL channels of an ASIC** (not the paddles): subtract the per-event pedestal from samples 0-3', None),
    (0, '20-36 % of beam-row waveforms saturate at 1023; 2-2.5 % sit on the tail of an earlier particle; the SFIL1 North board is noisy (40-50 ADC)', None),
], fig='figs/kcu0_sfil3_pedsub.png', layout='full',
    notes='decoding/kcu_hists.cc + kcu_hists_plot.py, wf_overlay.py; docs/method.md 4.10. Shown: SFIL3, pedestal-subtracted ADC vs sample for the eight channels (300 k triggers). '
          'The baseline shift is a detector-side common mode (rho 0.82-0.99 within an ASIC, 0.01 with the paddles); the HGCROC common-mode words are not exposed by the decoder.')

slide('TOA / TOT timing: what the HGCROC TDC gives in run 3062', [
    (0, '**KCU1 single SiPMs: the TDC works**, ~190 ps per channel from a self-calibration (pair MAD 0.27 ns); a per-channel threshold scan in the next run, then ~200 ps per end North - South is realistic', None),
    (0, '**KCU0 summing boards (SFIL2 beam cell):** the TOA is a TDC after a circular unwrap and a baseline correction, but noise words on 14-37 % of samples share the pulse codes: position along the fibres 220-250 mm per shower from the TOA against 120-146 mm from the ADC crossing', None),
    (0, 't_N and t_S correlate at +45 deg: the trigger is latched on the 25 ns clock, so a common event time dominates (the CAEN paddle time removes it). TOT appears only on saturated pulses: an amplitude on KCU1, nothing for timing', None),
], fig='figs/fig_xres.png', layout='full',
    notes='timing/toa_timing.py, toa_selfcal.py (KCU1), timing/toa_clean/ (KCU0 SFIL2 cell 3, ADC-free: the ADC only finds the pulse, t = 25 s + code-density table); docs/method.md 4.9, 4.19. '
          'Strict TOA cuts (quiet before the pulse, one word, pure code band) keep 17 % of showers and give 145-156 mm. Do not quote sub-ns North-South from KCU0 TOA in 3062. '
          'The KCU0 TOA configuration of the run (why 21 of 32 channels emit words on pedestal samples and 10 never do) is still to be obtained from the DAQ group.')

slide('AstroPix ToT calibration along the paper method', [
    (0, 'Row + column pairing (|dToA| <= 1, matching ToT): 543 k pixel hits; **chips 0-4 show the MIP Landau peak** (MPV 331-724 ToT units, per-pixel spread 23-37 %); chips 5-8 read ToT = 0', None),
    (0, 'No source scans for this module on disk: one-point anchor (MIP 34 keV in 128 um Si) gives 21-46 keV per chip with a representative v3 per-pixel calibration; the module configuration (threshold 202 mV) was found', None),
    (0, 'Bench source scans of module 8 would replace the anchor by the per-pixel polynomial', 'QUESTION'),
], fig='figs/astro_tot_calib.png', layout='full',
    notes='astropix/astro_tot_calib.py following arXiv:2602.06084. ToT unit 10 ns from the firmware clock. Imaging-layer results (matched clusters, beam height) on the performance slide.')

# ---- reconstruction
slide('Template fits recover the saturated pulses', [
    (0, '**One pulse template per channel** from clean unsaturated pulses; amplitude, time and pedestal fitted per pulse', None),
    (0, '**KCU0 pulses at the 10-bit ceiling** (up to 66 % of the SFIL3 beam cell at P5) are recovered from the unsaturated samples, to ~1800 ADC', None),
    (0, '**KCU1 single-SiPM pulses clip at ~940 ADC:** amplitude rebuilt from the time over threshold (slope band 0-0.76 ADC per unit)', None),
], fig='figs/tplfit_sfil3.png', layout='full',
    notes='kcu0_sfil_fit.C / kcu1_sfil5_fit.C. Pulse finder threshold 60 ADC, peak at samples 7-9, chi2/ndf < 8 with a 5 ADC noise and 5 % shape error; '
          'saturated samples are masked in the fit. Templates from 100-1000 ADC pulses; channels with < 100 clean pulses take the template from another run. '
          'Run 3064: 137 M KCU0 and 229 M KCU1 fitted pulses (14 + 5 condor parts). The KCU1 TOT slope cannot be measured without the LED runs; '
          'b = 0.40 with the band 0-0.76 is adopted from the SFIL5/SFIL4 width.')

slide('All August runs decoded, fitted and synchronized', [
    (0, 'CAEN alignment had failed for the 12 h run (anchored on an empty slice); the anchor is now chosen in the middle of the run: **3064 matches like 3062** (54 % within 100 ns, MAD 6.4 ns)', None),
    (0, 'KCU0 - KCU1 is one constant offset per run (92.7 % at 0 ticks, 7.3 % at -4); the clock ratio to CAEN is 1.3653407 in every run', None),
], table=[['run', 'type', 'length', 'KCU0 triggers', 'KCU1 offset [ticks]', 'CAEN matched', 'use'],
          ['3062', 'beam', '25 min', '2.9 M', '-1 143 665 946', '50 %', 'reference'],
          ['3063', 'beam', '2 min', '236 k', '357 099', '66 %', 'extends 3062 by 6 %'],
          ['3064', 'beam', '12 h', '69 M', '306 738', '54 %', 'main sample, 18 x 3062'],
          ['3066-3069', 'beam, 0.3-0.8 nA', '13-21 min each', '17-34 k', '198 368 ...', '99.5 %', 'trigger unknown'],
          ['3088', 'cosmics', '15 h', '1.04 M', '-203 777', '9.5 %', 'MIP sample; CAEN on its own trigger'],
          ['3091', 'cosmics', '20 min', '13 k', '306 604', '99.9 %', 'AstroPix on']],
    layout='table', notes='sync/run_sync.py per run -> output/sync_runs/run<RUN>/sync.txt. The host-clock epochs wander by 1.3-1.6 s over 12 h, so the coarse '
          'offset is only good to ~1 s; the anchor is the most significant of 15 slices over the middle 90 % of the run (score (n-b)/sqrt(b+1), +-10 s search). '
          'AstroPix: 3066 locks strongly (218 sigma), 3064 only in a wide scan (21.9 sigma, no 10 MHz phase lock at 3.4 sigma); 3063 / 3091 weak. '
          '3088: 912 k of its 1.01 M CAEN triggers have no KCU partner, so the CAEN had its own cosmic trigger.')

# ---- calibration
slide('Gains stable to better than 1 % over 12 h (run 3064)', [
    (0, '10-minute medians of the fitted beam-cell amplitude: **0.8-1.3 % rms per layer**, slopes below 0.13 % per hour', None),
    (0, 'The four layers move together (+-1.5 %) and the North / South ratio moves with them: that is the beam, not the gains', None),
    (0, 'The run-3062 calibration carries over to 3064 with no time correction', None),
], fig='figs/gain_stability_kcu0.png', layout='full',
    notes='calibration/stability/gain_stability.C + gain_stability_figs.py. Run structure from the single-paddle shower rate (median 952 Hz): beam-off gaps at 1.5-2.1 h and 5.8-6.3 h, '
          'beam until 9.08 h, then a low-current TAC period at 3-6 Hz. Pedestals of the beam cells within +-2.5 ADC. KCU1 SFIL5 South unclipped pulses: 0.6-0.8 % rms. '
          'Fraction of South beam-cell pulses at the ceiling: SFIL2 0.18, SFIL3 0.53, SFIL4 0.34.')

slide('Photostatistics: the SFIL2 / SFIL3 gain ratio is 1.49', [
    (0, 'Both fibre ends see the same photoelectrons: the variance between the ends grows linearly with the amplitude and gives the **ADC per photoelectron**', None),
    (0, 'South beam cells: **SFIL2 14.1 +- 0.2, SFIL3 9.5 +- 0.1 ADC / pe**; "equal gains" was wrong by 50 % between these layers', None),
    (0, 'SFIL4 North is anomalous (0.35 of South) and SFIL1 North is masked: those two layers stay unmeasured by this route', None),
], fig='figs/photostat_ns.png', layout='full',
    notes='sim/analysis/photostat_ns.py, method 4.13: var(x) - cov(x, y) = c + a A with a = ADC per photoelectron times the excess noise factor. '
          'The SFIL2 N/S amplitude ratio 0.85 is a gain ratio, not a light asymmetry. Single-ended readouts (SFIL5 South single SiPMs, CAEN) cannot be calibrated this way.')

slide('Attenuation and N / S gains from the AstroPix position', [
    (0, 'ln(A_N / A_S) against the shower position along the fibres: **lambda_eff = 1329 mm (SFIL2), 1646 mm (SFIL3)**; g_N / g_S = 0.85, 1.03, 0.29 (SFIL4 North weak)', None),
    (0, 'The geometric mean of the two ends cancels the attenuation: the response rises 10-20 % from P3 to P5, so the beam energy changes across the paddles', None),
], fig='figs/attenuation_top.png', layout='full',
    notes='calibration/attenuation/astro_x_attenuation.py: 45.6 k showers with an AstroPix cluster of purity >= 0.3 (chips 3, 5, 6, 8). The epic two-component attenuation '
          'predicts lambda_eff 2341 mm; the measured single exponential (1500 mm, v_eff 160 mm/ns from the N-S timing) is the v4b light model of the simulation. '
          'The SFIL4 slope (619 mm) is not an attenuation: its North end is the anomalous one.')

slide('Cell structure: a shared SiPM pair and an empty cell', [
    (0, '**Cells 1 and 2 of each layer are read by one SiPM pair on both channels** (correlation 0.97-0.99): use their sum', None),
    (0, '**Cell 4 carries 0.5 % of the beam cell** against 5-11 % in the model; the off-beam cells see 3-10 x less lateral light than the model', None),
    (0, 'Run 3064 gives the same cell ratios as 3062 within 2 %', None),
], fig='figs/cell_equalisation.png', layout='full',
    notes='calibration/cells/cell_equalisation.py (raw peaks, triggers with an unsaturated beam cell). v2 SiPM boards (all South, North SFIL2-4): cell 3 = SiPM3 = beam row, cell 4 = SiPM4, '
          'cells 1 and 2 both read SiPM1 + SiPM2 (B. Kim deck, 10/8 version). Model at beam y = -5 mm; the bands are the beam-y scan. '
          'Whether the lateral deficit is the beam height, the cell-4 channels or the model is open.')

slide('Energy scale: gains anchored to the model at 4 GeV', [
    (0, 'Per layer: data South amplitude over simulated South light, P6 showers, survey geometry (v4)', None),
    (0, 'The anchor energy is the weak point: +-1 GeV moves every gain by -+25 %; the measured light model (v4b) adds +4 %', None),
], table=[['element', 'gain, 4 GeV anchor', 'note'],
          ['SFIL1 South', '8.6 ADC / MeV', 'half of the others under every hypothesis'],
          ['SFIL2 South', '24.2 ADC / MeV', 'SFIL2 / SFIL3 = 1.32 (photostatistics 1.49)'],
          ['SFIL3 South', '18.4 ADC / MeV', 'beam row; cosmic scale 27 +- 2 (slide 20)'],
          ['SFIL4 South', '14.7 ADC / MeV', ''],
          ['SFIL5 South (KCU1)', '14.0 x SFIL4', 'from the SFIL5 / SFIL4 ratio per shower, TOT-recovered'],
          ['Baby BCAL (CAEN)', '220 / 270 / 127 / 72 per array', 'CAEN units per MeV, rows 0-3; fall with depth open']],
    layout='table', notes='sim/analysis/stack_shares_calib.py (sfil_gain_calibration_v4.txt), sfil5_sfil4_gain.py, bbcal_gains.py. ADC per MeV of simulated South light (attenuation included). '
          'At 3 GeV the SFIL gains are 9.2 / 28.1 / 22.9 / 19.6. g5/g4 = 14.0 (13.7-14.3 by paddle, 14.6 / 13.0 at 3 / 6 GeV). '
          'Baby BCAL per channel: A0 226, A1 218, A2 223; A5 510, A6 537, A7 572; A8 356, A9 368, A11 418; A14 289 CAEN units per MeV.')

slide('Baby BCAL: per-channel gains and a fall with depth', [
    (0, 'DRS4 cell correction applied to every run; the block sits **5 mm below nominal** (row patterns); within a row the arrays agree to 2-6 %', None),
    (0, '**The gain falls by a factor 3 from row 0 to row 3**, which the FTBF MIP peaks do not show: either less light reaches the deep rows than the model says, or the array sums are not one-to-one', None),
], fig='figs/bbcal_gains_yfine.png', layout='full',
    notes='caen_calib/caen_apply.py (per-cell DRS4 pedestals of 3062 applied to 3064, 3088, 3091), sim/analysis/bbcal_gains.py with a bbcalYOffset scan. '
          'Channels A3 / A4 dead, North silent. Rows 0-1 gains are usable now; rows 2-3 carry the model shape as a systematic. The 3088 cosmic data (CAEN on its own trigger) '
          'should show MIP peaks per array and decide.')

slide('SFIL5: clipped pulses and the SFIL5 / SFIL4 gain', [
    (0, 'KCU1 channels saturate at ~940 ADC above pedestal; the amplitude is rebuilt from the time over threshold', None),
    (0, 'SFIL5 / SFIL4 per shower against the model light ratio: **g5 / g4 = 14.0** (13.7-14.3 by paddle)', None),
    (0, 'No MIP or LED calibration is possible from these data (cosmics: nothing above 30 ADC)', None),
], fig='figs/sfil5_sfil4_gain_rec.png', layout='full',
    notes='kcu1_tot/kcu1_sat_recovery.py, sim/analysis/sfil5_sfil4_gain.py: A = 943 + b (TOT - T_c), b = 0.40 adopted with the band 0-0.76 (+8 % on g5/g4 at b = 0.76). '
          'The North board was dead in 3062 (4 live channels only). KCU1 ADC per photoelectron from the North pairs: 24-34, i.e. 2.5-3.5 x the KCU0 summed-channel chain.')

slide('Cosmic run 3088: real muons in SFIL3 and SFIL4', [
    (0, '**Muon tag:** both ends of one cell (independent ASICs) at the beam timing with the pulse shape; common-mode subtraction brings the noise to 0.8-2.5 ADC', None),
    (0, 'SFIL3: 350-530 tags per cell with an accidental expectation below 36; spectra peak at 15-20 ADC with Landau tails (shown: cells 2-4)', None),
    (0, 'The tags come with a paddle pulse in 97-100 % (52 % of all triggers): **near-horizontal muons through the upstream paddles**, 1.4 per minute', None),
], fig='figs/cosmic_gm_3panels.png', layout='full',
    notes='calibration/cosmics/cosmic_mip.C, cosmic_mip_fit.py on the full run (984 k triggers). The "late pulses" of the first look were baseline ramps shared by the connected channels of an ASIC half; '
          'the per-half median over connected channels removes them. SFIL1 and SFIL2 cells 3-4 tags are accidentals (noisy ends); the in-board "vertical" tags on SFIL1 S, SFIL2 S and SFIL5 N '
          'are within-board bursts, not muons. The tag rate is 45 x higher in triggers with a paddle pulse; the rate matches near-horizontal tracks through the 33 cm2 paddle faces into SFIL3, '
          'not a paddle laid on the stack (~30 per minute).')

slide('Cosmic MIP scale from a trigger-matched simulation', [
    (0, '**Geant4 cosmic mode:** every muon sent through a trigger paddle (flux-correct sampling); 10 M muons, 5-10 k paddle + SFIL3 crossings per cell', None),
    (0, 'Readout model per end: Poisson(light x Y) x ADC / pe + noise; the light yield Y is fitted to the two-end spectrum', None),
    (0, '**SFIL3: 26.5 / 27.8 / 29.2 ADC per MeV for cells 3, 4 and 1+2 = 27 +- 2**, all cells consistent (shown: cells 3, 4, 1+2); SFIL1, 2 and 4 give no scale', None),
], fig='figs/cosmic_sfil3_3panels.png', layout='full',
    notes='calibration/cosmics/cosmic_mip_sim.py on sim/cosmic/tgt (condor/cosmic_tgt.job, /bic/gun/cosmicTarget 1). ADC per pe from the photostatistics (SFIL3 9.5, North = South x 1.03); the fit is repeated for 0.7 and 1.4 x that (24-28). '
          'Template choice (isotropic plane vs through-paddle) moves the result by < 5 %. Cells 1 and 2 must be summed (shared SiPM pair); their per-channel fits give 18-20. '
          'SFIL1: accidental tags; SFIL2: noisy cells 3-4, too few in 1-2; SFIL4: North end at 0.29 of South, tag efficiency 1 %.')

slide('Two data handles: the P6 showers carry about 3 GeV', [
    (0, '**Cosmic scale** of the SFIL3 beam row, 27 ADC / MeV, against the 4 GeV anchor of 18.4: the showers deposit 2.6-2.9 GeV-equivalent light', None),
    (0, '**Longitudinal shape** of the SFIL2-4 shares (photostatistics-equalised) against the model: 3.0 GeV (3.0-3.3) for P3, P4, P6; 3.3 (3.0-4.0) for P5', None),
    (0, 'The energy rises P4 < P3 < P6 < P5 (double ratios, geometric mean); the pair-spectrometer momentum map would settle the absolute value', 'OPEN'),
], fig='figs/stack_shares_sfil_pe_v4.png', layout='left',
    notes='stack_shares_pe.py sim/v4: SFIL2 scaled by 9.47/14.09, SFIL1 and SFIL4 assumed at the SFIL3 gain (outer bar +-17.5 %). Saturation removed on both sides by quantile matching gives the same residuals (SFIL3 +0.06, SFIL4 -0.03). '
          'Same result with v3 and v4b. The 4 GeV figure came from the May paddle map; nothing in the August data supports it.')

# ---- simulation
slide('Geant4 model at v4: the August survey geometry', [
    (0, 'Pb / ScFi SFIL1-5 with the ePIC fibre matrix (19.05 mm cells, 505-519 mm fibres), the AstroPix module, the Baby BCAL and the paddles at the surveyed positions', None),
    (0, 'Readout: Birks-saturated light, measured attenuation and v_eff, earliest arrival time per end, AstroPix pixels; cosmic mode', None),
    (0, 'Open: the survey puts P3 and P6 9.5 mm apart, the AstroPix chip map 60 mm; to be settled with the survey author', 'OPEN'),
], fig='figs/event_display_4GeV.png', layout='left',
    notes='sim/README.md: v1 pure Pb (2026-10-22), v2 ePIC materials, v3 epic BeamTest branch geometry (cell-assignment bug fixed), v4 survey (MeasureMap2), v4b = v4 + measured light model. '
          'Kept from v3: active height 76.2 mm, beam y from the data (-5 mm), beam x from the AstroPix hits (+30 mm P6, +80 mm P5). '
          'Known / assumed list in the README; the metrology replaced most assumed values.')

slide('Energy scan 3-6 GeV with the v4 model', [
    (0, 'Sampling fraction 0.102 in the SFILs; 78.6 MeV of visible fibre energy per GeV', None),
    (0, 'Resolution of the visible energy of the full stack: **7.3 % / sqrt(E) + 1.4 %**', None),
    (0, 'The stack holds 77 % of the energy, 15-16 % escapes behind the Baby BCAL; shower maximum in SFIL4 at 4 GeV', None),
], fig='figs/sim_summary_v4.png', layout='left',
    notes='sim/analysis/sim_summary.py tag=_v4 (5000 e+ per point, x = +30 mm; +80 mm points for P5). Against v3 the Baby BCAL share rises 6-7 % (survey z and the 51 mm gap); '
          'the SFIL1-4 layer shares at 4 GeV are 0.047 / 0.203 / 0.356 / 0.393. N-S timing rms 0.17 ns; 0.7 AstroPix pixels per event on average.')

# ---- performance
slide('Run 3062: linearity and resolution per paddle', [
    (0, 'E_vis = SFIL1-5 + Baby BCAL with the anchored gains and the attenuation correction; 100 k showers seen by all three boards', None),
    (0, '**Linear within +-3 %** over 3.5-4.5 GeV; **sigma / mu = 13.2 / 12.2 / 10.9 / 10.5 %** (P4 / P3 / P6 / P5) against 5 % in the model: the energy spread within a paddle (9-12 %) dominates', None),
    (0, 'Baby BCAL share 15-20 % below the model (the deep-row gains); SFIL1-4 alone: 19-30 %', None),
], fig='figs/performance_top.png', layout='full',
    notes='sim/analysis/performance.py. Paddle energies: P6 = 4.00 GeV anchor, P5 4.50 and P3 4.00 from the geometric-mean curve, P4 3.45 from the SFIL2-5 profile; '
          'the energy per paddle is known to 10-20 % and limits every statement. Core sigma/mu: iterated +-2 sigma. With the 3 GeV scale of slide 15 the visible energies scale by 3/4; '
          'the relative numbers do not change.')

slide('Run 3064: the same performance, 17 x the statistics', [
    (0, '34 M showers; SFIL1-4 alone: **sigma / mu = 29.6 / 26.7 / 20.6 / 19.1 %** (P4 / P3 / P6 / P5), equal to 3062 within 0.3 %', None),
    (0, 'Per 40-minute part the mean visible energy varies by 0.3-0.4 % and the width by 0.1-0.2 %: no drift over the run', None),
    (0, 'SFIL5 and the Baby BCAL join next (CAEN aligned, cells corrected) for the full-stack resolution', 'NEXT'),
], fig='figs/performance_3064.png', layout='full',
    notes='sim/analysis/performance_3064.py on the 14 layer-sum caches (data_fit_layers_3064_partPP.npz). Layer shares identical to 3062 within 0.001; '
          'shape energy 3.0 GeV (P5 3.9). The TAC period (~50 k triggers) sits at the end of part 13 and is not separated.')

slide('Imaging layer: single hits, 2.8 mm beam height', [
    (0, '5.7 % of the showers have a matched AstroPix cluster (3.4 % at purity >= 0.3, 1.6 % at >= 0.7)', None),
    (0, 'Clusters are **one pixel in 67 %, two in 31 %**: no charge sharing at 0.5 mm pitch for a 3-6 GeV positron', None),
    (0, 'Beam height 2.8 mm rms, the same on chips 3, 5, 6 and 8; used for the attenuation study and the shower position (140 mm per shower along the fibres)', None),
], fig='figs/astro_imaging.png', layout='full',
    notes='astropix/astro_imaging.py on 2.21 M showers. Matched fractions per paddle: P4 1.7 % (chip 0 never passes the purity cut), P3 8.0 %, P6 10.3 %, P5 3.0 %. '
          'Chip 0 records 625 k hits against 240-360 k for chips 1-8; the excess is not beam-correlated. Hot-pixel scan not started.')

slide('Plan versus status', [], table=[
    ['step', 'status'],
    ['1  event synchronization', 'done for every run; AstroPix statistical (80 % purity on three chips)'],
    ['2  TOA / TOT timing North - South', 'not possible in 3062 (KCU0 TOA not pulse-timed); recipe ready for a run with thresholds'],
    ['3  waveform calibration', 'pedestal, pile-up flag, template fits, saturation recovery done; KCU1 TOT slope needs the LED runs'],
    ['4  AstroPix hot-pixel scan', 'not started'],
    ['5  sampling fraction and energy scale', 'SF 0.102 from the model; gains anchored at 4 GeV; cosmic scale and shape say ~3 GeV'],
    ['6  shower reconstruction, position, leakage', 'layer shares, N/S position (140 mm), AstroPix clusters, containment done; no full reconstruction yet'],
    ['7  physics performance', 'linearity, resolution and containment per paddle for 3062 and 3064 at the present calibration']],
    layout='table', notes='docs/status.md section 8, updated.')

slide('Next steps', [
    (0, '**1. Full-stack performance of 3064:** join SFIL5 and the Baby BCAL (aligned, cell-corrected) to the 34 M showers', 'NEXT'),
    (0, '**2. Baby BCAL MIP peaks from the 3088 cosmics** per array: settles the factor-3 fall with depth', 'NEXT'),
    (0, '**3. Re-anchor the gains at 3 GeV** and recheck the SFIL5 and Baby BCAL double ratios across the paddles', 'NEXT'),
    (0, '**4. Decoding and synchronization note** (sync/note) with the per-run constants', 'NEXT'),
    (0, '**5. AstroPix hot-pixel scan** and the SFIL5 South map when the pin map arrives', 'TODO'),
], layout='text', notes='Order by leverage. Items 1-3 use data already on disk; 5 waits for external input.')

slide('Questions to the group', [
    (0, '**Beam energy:** a pair-spectrometer momentum map for August, or any run with a known energy; the data say ~3 GeV under P6', 'QUESTION'),
    (0, '**Survey versus data:** P3 and P6 are 9.5 mm apart in the survey and 60 mm apart in the AstroPix chip map', 'QUESTION'),
    (0, '**Trigger definitions** of the cosmic runs (paddle OR plus what?) and of the low-current runs 3066-3069', 'QUESTION'),
    (0, '**Missing runs:** cosmics 3089 (47 h), pedestal 3173-3178, LED 3179-3190, injection run 2804', 'QUESTION'),
    (0, '**SFIL5 South SiPM board pin map** and the Baby BCAL FTBF MIP constants per array', 'QUESTION'),
    (0, '**DAQ for the next run:** socket buffer / recvmmsg against the UDP loss, the common busy, the KCU0 TOA configuration and the DAQ YAML (analog gains)', 'QUESTION'),
], layout='text', notes='The LED runs would give the KCU1 TOT slope and the SFIL5 absolute scale; the cosmic run 3089 triples the MIP sample.')

# ---- backup
slide('Backup: synchronization per run', [], table=[
    ['run', 'KCU1 offset [ticks]', 'KCU0 - KCU1 at 0 / -4', 'CAEN clock ratio', 'CAEN < 100 ns', 'AstroPix'],
    ['3062', '-1 143 665 946', '92.6 / 7.4 %', '1.365340679', '49.7 %', '19.3 sigma, phase 189'],
    ['3063', '357 099', '92.7 / 7.3 %', '1.365340786', '66.0 %', 'weak (4-6 sigma)'],
    ['3064', '306 738', '92.7 / 7.3 %', '1.365340765', '54.5 %', '21.9 sigma wide scan, no phase lock'],
    ['3066', '198 368', '92.7 / 7.3 %', '1.365340654', '99.5 %', '218 sigma, phase 2'],
    ['3067-3069', '1 113 501 / 559 073 / 466 275', '92.7 / 7.3 %', '1.36534066-68', '99.5 %', ''],
    ['3088', '-203 777', '89 / 0 % (11 % elsewhere)', '1.365340927', '9.5 %', 'off'],
    ['3091', '306 604', '89 / 0 %', '1.365340859', '99.9 %', 'weak']],
    layout='table', notes='output/sync_runs/run<RUN>/sync.txt. CAEN residual MAD 6.4-6.6 ns in every run; drift +5.4-5.6 ppm.')

slide('Backup: the v4 geometry against v3', [], table=[
    ['element', 'v3 (epic BeamTest file)', 'v4 (August survey)'],
    ['SFIL1-5 fibre length', '510 mm, all', '519.1 / 504.8 / 517.5 / 504.8 / 504.8 mm'],
    ['SFIL1-5 x centre', '0', '-7.6 / -10.0 / -27.5 / -16.4 / -29.1 mm'],
    ['SFIL1-5 upstream face z', '0 / 65.1 / 130.1 / 195.2 / 260.2 mm', '0 / 64.3 / 130.2 / 195.3 / 258.0 mm'],
    ['carbon frames in the slots', '2 x 2 mm per slot', 'none'],
    ['AstroPix', 'x centre 0, z = -21.6 mm', 'x -90 .. +90 mm, downstream face z = -6.5 mm'],
    ['Baby BCAL', '61 mm behind SFIL5, x centre 0', '51.1 mm behind SFIL5, x centre -33.3 mm'],
    ['paddles x (P3 P4 P5 P6)', '-30 / -80 / +80 / +30 mm (chip map)', '+3.5 / -74.3 / +100.3 / +13.0 mm'],
    ['paddles z centre', '-180 / -60 / -180 / -60 mm', '-166.0 / -66.0 / -183.5 / -64.4 mm']],
    layout='table', notes='sim/README.md "Geometry v4". Kept: active height 76.2 mm, beam y -5 mm, beam x from the AstroPix hits. Energy scan results: v3 and v4 shares within 0.01; implied gains 7.8 / 22.8 / 18.0 / 14.8 (v3) -> 8.6 / 24.2 / 18.4 / 14.7 (v4) -> 9.0 / 25.5 / 19.2 / 15.4 (v4b).')

slide('Backup: why the cosmic muons are near-horizontal', [
    (0, 'The paddles are vertical slabs 60-180 mm upstream of the stack at beam height; a vertical muon never crosses a paddle and a SFIL cell', None),
    (0, 'In the data a paddle pulse sits on 52 % of all triggers but on 97-100 % of the SFIL3 muon tags, over all four paddles; the tag rate is 45 x higher with a paddle pulse', None),
    (0, 'SFIL3 tags come with SFIL2 South activity in ~60 % (1 % in random triggers): tracks through SFIL1-3 in sequence; the layers are side by side in z, so a downward muon crosses one layer only', None),
    (0, 'Rate: 1.4 SFIL3 tags per minute = near-horizontal flux through 33 cm2 of paddle face into 0.6 sr of SFIL3; a paddle laid flat on the stack would give ~30 per minute', None),
    (0, 'Simulation: 2 M muons from a plane above the stack gave 2-12 paddle + SFIL3 crossings per cell; the through-paddle mode gives 0.16 % of 10 M', None),
], layout='text', notes='Median zenith angle 72 deg (10-90 %: 64-79). The South light per SFIL3 crossing is 1.08 MeV against 1.13 MeV for the plane production: the template changes little.')

slide('Backup: saturation handled on both sides', [
    (0, 'The South beam-row channel sits at the 10-bit ceiling in 29 / 54 / 32 % of the showers for SFIL2 / 3 / 4 (rising P4 < P3 < P6 < P5); the template fit recovers those pulses', None),
    (0, 'Below the ceiling raw peaks and fitted amplitudes agree to 0.004 in every layer share: the fit is unbiased', None),
    (0, 'For the model comparison the same top fractions are dropped from the model light (quantile matching): the residuals are the same as for the full sample (SFIL3 +0.06, SFIL4 -0.03)', None),
    (0, 'The ceiling is an absolute reference: the gain that reproduces the saturation fraction is 1.3 x the gain from the mean amplitude, in every layer and energy', None),
], layout='text', notes='stack_shares_pe.py nosat, sat_gain.py. A threshold-based unsaturated sample is biased against large pulses: it kept 80 % of the model against 30 % of the data.')

slide('Backup: the setup as read out in run 3062', [], fig='figs/kcu0_channel_map.png', layout='full',
    notes='docs/kcu0_channel_map.png from the 8/13 hookup table, verified channel by channel in the data. Cells 1 and 2 share SiPM1 + SiPM2 on v2 boards.')

slide('Backup: code and data', [], table=[
    ['folder (github.com/tmengel/bic_testbeam)', 'content'],
    ['sync/', 'run_sync.py (any run), decoders, event tree, the decoding / sync note'],
    ['waveforms/', 'template fits kcu0_sfil_fit.C, kcu1_sfil5_fit.C, condor jobs for every run'],
    ['calibration/', 'stability/, cosmics/, cells/, attenuation/, caen_calib/, kcu1_tot/'],
    ['sim/', 'Geant4 model (BicGeometry.hh, macros, condor scans, cosmic mode), analysis/ (shares, gains, performance)'],
    ['decoding/, timing/, decoder_tests/', 'ROOT converters and checks, raw-format tools, TOA / TOT studies, decoder validation'],
    ['astropix/', 'imaging layer, ToT calibration, geometry checks'],
    ['docs/', 'status.md (overview, read first), method.md, data_format.md, completeness_plan.md, new_runs.md, this deck; collaborator write-up as a Claude Doc'],
    ['outputs (SDCC)', '/sphenix/tg/tg01/jets/tmengel/eic/{kcu0,kcu1,caen,sim}/ : fits, stability, cosmics, CAEN calibration, simulation']],
    layout='table', notes='')

# ---------------------------------------------------------------- helpers
def parse_bold(text):
    """split on ** : [(segment, bold)]"""
    out, bold = [], False
    for i, seg in enumerate(text.split('**')):
        if seg: out.append((seg, bold))
        bold = not bold
    return out

def title_pt(title):
    for pt in (TITLE_PT, 36, 32, 28):
        if pil_font(False, pt).getlength(title) <= 12.5 * PX: return pt
    return 28
font_cache = {}
def pil_font(bold, pt):
    key = (bold, pt)
    if key not in font_cache: font_cache[key] = ImageFont.truetype(FONT_B if bold else FONT_R, int(round(pt * 2)))
    return font_cache[key]
def text_width(segs, pt):
    return sum(pil_font(b, pt).getlength(s) for s, b in segs)
def wrap_segments(segs, pt, width_px):
    """wrap a list of (text, bold) segments into lines of segments, by word, measured with the preview font"""
    words = []
    for s, b in segs:
        parts = re.split(r'(\s+)', s)
        for p in parts:
            if p: words.append((p, b))
    lines, cur, curw = [], [], 0.
    for w, b in words:
        ww = pil_font(b, pt).getlength(w)
        if curw + ww > width_px and cur and not w.isspace():
            lines.append(cur); cur, curw = [], 0.
            if w.isspace(): continue
        if not cur and w.isspace(): continue
        cur.append((w, b)); curw += ww
    if cur: lines.append(cur)
    return lines

def bullet_lines(bullets, width_in):
    """number of wrapped text lines a bullet block takes, and its height in inches"""
    n, h = 0, 0.
    for lvl, text, tag in bullets:
        pt = BODY_PT if lvl == 0 else SUB_PT; indent = 0.3 + 0.35 * lvl
        segs = parse_bold(text) + ([('  ' + tag, True)] if tag else [])
        lines = wrap_segments(segs, pt, (width_in - indent) * PX); n += len(lines); h += len(lines) * pt * 1.25 / 72. + 0.08
    return n, h

# ---------------------------------------------------------------- pptx
prs = Presentation(TEMPLATE)
# drop the template's own slides, keep layouts and master
sldIdLst = prs.slides._sldIdLst
for sldId in list(sldIdLst):
    prs.part.drop_rel(sldId.rId); sldIdLst.remove(sldId)
L_TITLE, L_TITLEONLY = prs.slide_layouts[0], prs.slide_layouts[5]

def add_footer(sl, num):
    for x, txt, align in ((0.0, DATE, PP_ALIGN.LEFT), (10.33, str(num), PP_ALIGN.RIGHT)):
        tb = sl.shapes.add_textbox(Inches(x), Inches(7.09), Inches(3.0), Inches(0.4)); tf = tb.text_frame; tf.margin_left = tf.margin_right = Inches(0.1)
        p = tf.paragraphs[0]; p.alignment = align; r = p.add_run(); r.text = txt; r.font.size = Pt(FOOT_PT); r.font.name = 'Arial'; r.font.color.rgb = RGBColor(0x40, 0x40, 0x40)

def add_bullets(sl, bullets, x, y, w, h):
    tb = sl.shapes.add_textbox(Inches(x), Inches(y), Inches(w), Inches(h)); tf = tb.text_frame; tf.word_wrap = True
    tf.margin_left = tf.margin_right = Inches(0.05); tf.margin_top = Inches(0.02)
    first = True
    for lvl, text, tag in bullets:
        p = tf.paragraphs[0] if first else tf.add_paragraph(); first = False
        pt = BODY_PT if lvl == 0 else SUB_PT; p.level = lvl; p.space_after = Pt(6)
        pPr = p._p.get_or_add_pPr(); pPr.set('marL', str(int(Inches(0.3 + 0.35 * lvl)))); pPr.set('indent', str(int(-Inches(0.25))))
        bu = pPr.makeelement(qn('a:buChar'), {'char': '•' if lvl == 0 else '–'});
        for e in pPr.findall(qn('a:buNone')): pPr.remove(e)
        pPr.append(bu)
        for seg, bold in parse_bold(text):
            r = p.add_run(); r.text = seg; r.font.size = Pt(pt); r.font.name = 'Arial'; r.font.bold = bold; r.font.color.rgb = RGBColor(0, 0, 0)
        if tag:
            r = p.add_run(); r.text = '  ' + tag; r.font.size = Pt(pt); r.font.name = 'Arial'; r.font.bold = True; r.font.color.rgb = RGBColor.from_string(TAGCOL[tag])
    return tb

def add_picture_fit(sl, path, x, y, w, h):
    im = Image.open(path); iw, ih = im.size; s = min(w / iw, h / ih); pw, ph = iw * s, ih * s
    px = x + (w - pw) / 2; sl.shapes.add_picture(path, Inches(px), Inches(y), Inches(pw), Inches(ph)); return (px, y, pw, ph)

def add_table(sl, rows, x, y, w, pt=TABLE_PT):
    nr, nc = len(rows), len(rows[0]); rowh = pt * 2.2 / 72.
    shp = sl.shapes.add_table(nr, nc, Inches(x), Inches(y), Inches(w), Inches(rowh * nr)); tbl = shp.table
    # column widths from the longest entries
    lens = [max(len(r[c]) for r in rows) for c in range(nc)]; tot = sum(min(l, 60) for l in lens)
    for c in range(nc): tbl.columns[c].width = Inches(w * min(lens[c], 60) / tot)
    tblPr = tbl._tbl.tblPr; sid = tblPr.find(qn('a:tableStyleId'))
    if sid is None: sid = tblPr.makeelement(qn('a:tableStyleId'), {}); tblPr.append(sid)
    sid.text = '{5940675A-B579-460E-94D1-54222C63F5DA}'      # No Style, Table Grid
    for i, r in enumerate(rows):
        for c, v in enumerate(r):
            cell = tbl.cell(i, c); cell.margin_left = cell.margin_right = Inches(0.06); cell.margin_top = cell.margin_bottom = Inches(0.03)
            cell.fill.solid(); cell.fill.fore_color.rgb = RGBColor(0xF2, 0xF2, 0xF2) if i == 0 else RGBColor(0xFF, 0xFF, 0xFF)
            tf = cell.text_frame; tf.word_wrap = True; p = tf.paragraphs[0]; run = p.add_run(); run.text = v; run.font.size = Pt(pt); run.font.name = 'Arial'; run.font.bold = (i == 0); run.font.color.rgb = RGBColor(0, 0, 0)
    return shp

layout_record = []     # for the preview: per slide the boxes actually used
for i, s in enumerate(S):
    num = i + 1; rec = dict(num=num, title=s['title'], items=[])
    if s['layout'] == 'title':
        sl = prs.slides.add_slide(L_TITLE); sl.shapes.title.text = s['title']
        sub = sl.placeholders[1]; sub.text = s['subtitle']
        for p in sub.text_frame.paragraphs:
            for r in p.runs: r.font.size = Pt(24); r.font.name = 'Arial'
        rec['subtitle'] = s['subtitle']; layout_record.append(rec); continue
    sl = prs.slides.add_slide(L_TITLEONLY); t = sl.shapes.title; t.text = s['title']; tpt = title_pt(s['title']); rec['tpt'] = tpt
    for p in t.text_frame.paragraphs:
        for r in p.runs: r.font.size = Pt(tpt); r.font.name = 'Arial'; r.font.bold = False
    add_footer(sl, num)
    top = 1.3
    if s['layout'] == 'text':
        add_bullets(sl, s['bullets'], 0.35, top, 12.6, 5.6); rec['items'].append(('bullets', s['bullets'], 0.35, top, 12.6))
    elif s['layout'] == 'left':
        tw = 5.0; add_bullets(sl, s['bullets'], 0.35, top, tw, 5.6); rec['items'].append(('bullets', s['bullets'], 0.35, top, tw))
        box = add_picture_fit(sl, s['fig'], 5.55, top, 7.5, 5.6); rec['items'].append(('fig', s['fig']) + box)
    elif s['layout'] == 'full':
        n, h = bullet_lines(s['bullets'], 12.6); add_bullets(sl, s['bullets'], 0.35, top, 12.6, h + 0.1); rec['items'].append(('bullets', s['bullets'], 0.35, top, 12.6))
        fy = top + h + 0.15; box = add_picture_fit(sl, s['fig'], 0.35, fy, 12.65, 6.95 - fy); rec['items'].append(('fig', s['fig']) + box)
    elif s['layout'] == 'table':
        y = top
        if s['bullets']:
            n, h = bullet_lines(s['bullets'], 12.6); add_bullets(sl, s['bullets'], 0.35, top, 12.6, h + 0.1); rec['items'].append(('bullets', s['bullets'], 0.35, top, 12.6)); y = top + h + 0.2
        add_table(sl, s['table'], 0.35, y, 12.65); rec['items'].append(('table', s['table'], 0.35, y, 12.65))
    if s['notes']: sl.notes_slide.notes_text_frame.text = s['notes']
    layout_record.append(rec)
prs.save(OUT); print('wrote', OUT, len(S), 'slides')

# ---------------------------------------------------------------- preview (same geometry, Nimbus Sans for Arial)
os.makedirs('preview', exist_ok=True)
def draw_runs(d, x, y, segs, pt, color=(0, 0, 0)):
    for seg, bold in segs:
        f = pil_font(bold, pt); d.text((x, y), seg, font=f, fill=color); x += f.getlength(seg)
for rec in layout_record:
    im = Image.new('RGB', (1920, 1080), 'white'); d = ImageDraw.Draw(im)
    if 'subtitle' in rec:
        f = pil_font(False, 60); tw = f.getlength(rec['title']); d.text(((1920 - tw) / 2, 330), rec['title'], font=f, fill='black')
        for k, line in enumerate(rec['subtitle'].split('\n')):
            f2 = pil_font(False, 24); d.text(((1920 - f2.getlength(line)) / 2, 540 + 56 * k), line, font=f2, fill='black')
        d.text((40, 1000), '[CU Boulder logo from the template]', font=pil_font(False, 12), fill='gray')
        im.save('preview/slide%02d.png' % rec['num']); continue
    tp = rec.get('tpt', TITLE_PT); d.text((0.1 * PX + 14, 0.03 * PX + 30 + (TITLE_PT - tp)), rec['title'], font=pil_font(False, tp), fill='black')
    f = pil_font(False, FOOT_PT); d.text((0.1 * PX, 7.09 * PX + 6), DATE, font=f, fill=(64, 64, 64)); nt = str(rec['num']); d.text((13.23 * PX - f.getlength(nt), 7.09 * PX + 6), nt, font=f, fill=(64, 64, 64))
    for it in rec['items']:
        if it[0] == 'bullets':
            _, bullets, x, y, w = it; yy = y * PX + 4
            for lvl, text, tag in bullets:
                pt = BODY_PT if lvl == 0 else SUB_PT; indent = 0.3 + 0.35 * lvl; segs = parse_bold(text)
                if tag: segs = segs + [('  ' + tag, True)]
                lines = wrap_segments(segs, pt, (w - indent) * PX); lh = pt * 2 * 1.25
                d.text(((x + indent - 0.22) * PX, yy), '•' if lvl == 0 else '–', font=pil_font(False, pt), fill='black')
                for ln in lines:
                    xx = (x + indent) * PX
                    for seg, bold in ln:
                        col = tuple(int(TAGCOL[tag][k:k + 2], 16) for k in (0, 2, 4)) if (tag and bold and seg.strip() and seg.strip() in tag) else (0, 0, 0)
                        fnt = pil_font(bold, pt); d.text((xx, yy), seg, font=fnt, fill=col); xx += fnt.getlength(seg)
                    yy += lh
                yy += 6 * 2
        elif it[0] == 'fig':
            _, path, x, y, w, h = it; pic = Image.open(path).convert('RGB'); pic = pic.resize((max(1, int(w * PX)), max(1, int(h * PX))), Image.LANCZOS); im.paste(pic, (int(x * PX), int(y * PX)))
        elif it[0] == 'table':
            _, rows, x, y, w = it; nc = len(rows[0]); lens = [max(len(r[c]) for r in rows) for c in range(nc)]; tot = sum(min(l, 60) for l in lens)
            cw = [w * min(lens[c], 60) / tot * PX for c in range(nc)]; yy = y * PX; rowh = TABLE_PT * 2.2 / 72. * PX
            for i, r in enumerate(rows):
                xx = x * PX; hmax = rowh
                for c, v in enumerate(r):
                    lines = wrap_segments([(v, i == 0)], TABLE_PT, cw[c] - 18); hmax = max(hmax, len(lines) * TABLE_PT * 2 * 1.2 + 10)
                for c, v in enumerate(r):
                    d.rectangle([xx, yy, xx + cw[c], yy + hmax], outline=(180, 180, 180), fill=(242, 242, 242) if i == 0 else 'white')
                    lines = wrap_segments([(v, i == 0)], TABLE_PT, cw[c] - 18); ly = yy + 5
                    for ln in lines: draw_runs(d, xx + 9, ly, ln, TABLE_PT); ly += TABLE_PT * 2 * 1.2
                    xx += cw[c]
                yy += hmax
    im.save('preview/slide%02d.png' % rec['num'])
print('previews in preview/')
