# BIC test beam, run 3062: status of all studies (2026-10-23, new runs and to-do added 2026-10-28)

One page per topic would be `method.md` (synchronization and detector studies), `root_files.md` (the data
in ROOT trees), `subsystems.md`, `beam.md`, `channel_map.md`, `data_format.md`, `pitfalls.md`,
`completeness_plan.md`, and `sim/README.md` (the Geant4 model). This page is the overview: what
exists, what was found, what is open. Numbers are for run 3062 (Aug 13-14 2026, 4.5 SFILs + Baby BCAL,
3-6 GeV positrons) unless stated.

## 1. Where the data and tools are

| item | location |
| --- | --- |
| raw RCDAQ files (15, 36 GB) | `eic/data/` (symlink `bic_testbeam/data/raw`) |
| raw files of beam runs 3061, 3063-3069 and cosmic runs 3088, 3090, 3091 (0.73 TB, added 2026-10-27) | `/gpfs01/eic/data1/BIC_Beamtest/`, catalog and usage in `new_runs.md` |
| decoders (local forks with the TOT, CRC/Hamming and reorder fixes) | `eic/online_distribution` (branch tanner), `eic/ROC_pool`, installed in `eic/install`, `setup.sh` |
| every stream as per-trigger ROOT trees | `/sphenix/tg/tg01/jets/tmengel/eic/{kcu0,kcu1,caen,astropix}/`, symlinks `bic_testbeam/data/*_3062.root`, format in `root_files.md` |
| synchronization pipeline and event tree | `bic_testbeam/sync/` (`reproduce.sh`), `data/events_3062.root` (4.53 M common triggers), `method.md` 4.6 |
| analysis macros | `bic_testbeam/common/`: `bic_streams.{py,C}` (multi-stream reader + matching); `bic_testbeam/waveforms/`: `kcu0_sfil_plots.C` (waveform / TOA / TOT panels with quality flags), `kcu0_sfil_fit.C` (template pulse fits), `kcu1_sfil5_fit.C` + `kcu1_sfil5_plots.C` (the same for the 64 KCU1 SFIL5 South channels: fits with TOT/TOA per channel, and the figure set: templates, amplitudes, fit vs raw, chi2, TOT vs amplitude, array sums, waveform gallery), `kcu1_sfil5_raw.C` (raw galleries of BOTH SFIL5 ends: waveform, TOA code and TOT vs sample per channel, South 4 x 16 and North 2 x 72, plus the TOA pulse-time distribution; figures in `waveforms/kcu1_test_figs/`) |
| full-run pulse fits | `/sphenix/tg/tg01/jets/tmengel/eic/kcu0/fit_full/kcu0_sfil_fit_full.root` + figures |
| Geant4 model and energy scans | `eic/bic_testbeam/sim/` (README), outputs in `/sphenix/tg/tg01/jets/tmengel/eic/sim/` |
| decoder validation | `eic/decoder_tests/` (README with the numbers) |
| condor jobs | `bic_testbeam/sync/condor/`, `sim/condor/` (plain `condor_submit`, no DAGMan) |

External sources (ePIC geometry XML, Baby BCAL and AstroPix papers, two BIC talks) are kept as local copies in
`context/` and listed with URLs and what each is used for in `context/SOURCES.md` (2026-10-23).

## 2. Decoders and data integrity: done

- H2GCROC3 format understood bit by bit (`data_format.md`); the stock decoder kept only 6 of the 10 TOT
  bits, placed out-of-order lines into the wrong sample and checked no CRC: all three fixed in the local
  forks (stock decoder and ROC_pool), validated bit-exact against the raw bytes (75.5 M channel words per
  board, 0 mismatches) and against the old decoder where both are right (`decoder_tests/README.md`).
- Completeness: complete 20-sample waveforms 61 -> 75 % per packet (KCU0), 48 -> 79 % (KCU1); with the pool
  across packets 88 % (KCU0) and 83 % (KCU1) of recorded triggers. The remaining loss is UDP datagram loss
  in the DAQ (55.6 M datagrams lost on KCU0), not decoding; DAQ-side remedies in `completeness_plan.md`.
- Error counters per trigger and per link are in the ROOT files (`header.crc_ok`, `ham`, `present`) and
  in `kcu0_sfil_plots.C`: on KCU0 no CRC failures on the live board, Hamming flags on 2.5 % of triggers,
  mostly on the paddle half and the SFIL North half; KCU1 CRC failures confined to the dead SFIL5 North board.

**What the central install has (checked 2026-10-23 against `/eic/data/purschke/install`, built 21 Oct 17:20 / 18:46 from
`/eic/data/purschke/eic/{online_distribution,ROC_pool}`, commits 826f14b..a76c5c5 and d14a5e8..b48cf14):**

| fix | local fork (uncommitted on branch `tanner`, `ROC_pool` working tree) | central install |
| --- | --- | --- |
| TOT mask 0x3f -> 0x3ff | yes | **yes** (a41dd3e decoder, d14a5e8 pool, 21 Oct) |
| last sample of every Packet kept | yes (his 826f14b) | **yes** (20 Oct) |
| Hamming flags H1-H3 exposed | `iValue(wf,s,ch,"HAMMING")` per line + cumulative `"HAMMINGERRORS"` | **yes, as `iValue(wf,s,ch,"ERROR")`** per half-chip (035b6f3 decoder, cfb600f pool); different key name, same bits |
| CRC-32 of every frame recomputed and exposed (`"CRCERROR"`, `"CRCERRORS"`, `"BADFRAMES"`) | yes (`crc32_frame`) | **no**: the frame's CRC word is stored (`crc32`) but never checked |
| out-of-order lines placed by timestamp (`REORDER_WINDOW` = 4 x 164 clocks, `close_samples`, `"LATELINES"`) | yes, decoder and pool | **no**: a new sample opens whenever the timestamp jumps by more than 5 ticks in arrival order, so a late line still fragments the waveform (61 % instead of 75 % complete per packet on KCU0, 48 % instead of 79 % on KCU1) |
| cumulative counters `"BADLINES"`, `"RESETS"`, `"DECODEERRORS"` | yes | no |

So the TOT, last-sample and Hamming fixes are central; the CRC check and the reorder window are only here.

**Rebased 2026-10-29 (user asked):** both `tanner` branches are now committed on top of the current upstream heads
(`online_distribution` ePIC ce1cad2, `ROC_pool` main b48cf14), two commits each (build-product `.gitignore`s; the decoder
changes). The Hamming conflict is resolved to upstream: Martin's `error[]` field, `iValue(wf,s,ch,"ERROR")` and the ERR row
in `dump()`; the local `"HAMMING"` accessor is gone (`ts_dump`, `decoder_tests` renamed). Kept on top: the reorder window,
the marker test fix (`||`), UDP frame and line bounds checks, CRC-32 (`"CRCERROR"`, `"CRCERRORS"`) and the counters
(`"HAMMINGERRORS"`, `"BADFRAMES"`, `"BADLINES"`, `"LATELINES"`, `"RESETS"`, `"DECODEERRORS"`). Pre-rebase state: branch
`tanner-pre-rebase` in each repo. Split the same day into three commits per repo on top of the `.gitignore` commit: bug fixes (marker
test, half range, UDP frame and line bounds), reorder window, CRC-32 and counters (`online_distribution` a419f72, cc3c42f,
c5153ca; `ROC_pool` b671dce, 44540e4, 68f45c2); pre-split state in branch `tanner-presplit`. Not pushed; PRs to `sPHENIX-Collaboration/online_distribution` (ePIC) and
`mpurschke/ROC_pool` (main) still to be opened.

## 3. Synchronization: KCU0-KCU1 and CAEN final, AstroPix statistical

- KCU0 - KCU1: exact tick offset -1,143,665,946 (both at 160 MHz, no common reset), residual 0 or +4
  ticks; 927 k pairs. The 7.147 s is the difference of the two counter epochs, not a latency.
- CAEN - KCU0: group trigger time (8.533 ns) x 1.365340679 + an offset of 12.302 s, with a few-us wander
  followed by a running median; 59 % of CAEN triggers within 100 ns, MAD 6.4 ns (`bic_streams`).
- AstroPix - KCU0: the on-chip 8-bit timestamp is a 10 MHz counter phase-locked to KCU0 (phase 189, 19 sigma,
  same on all 9 chips); the readout stamp lags the trigger by 19 or 88 +- 5 ms, so per-event matching needs
  the paddle pattern and shower features: 188 k candidates, about 80 k true, 26 k at purity >= 0.8, all on
  chips 3, 5, 6 (a selection effect of the matcher, not of the beam).
- Trigger bookkeeping: `TRIGGER_OUT` is the board trigger number (+1 per trigger, 4.10 kHz, never reset);
  it counts missing triggers exactly; `EVENT_COUNTER` is per-link and unusable; `TRIGGER_IN` is 0.

## 4. The trigger and the beam (`beam.md`, `method.md` 4.11)

- Paddle OR at 4.10 kHz, 6.06 M triggers issued; KCU0 recorded 3.77 M (57 % complete of issued), KCU1 and
  CAEN about 1.0 kHz each. 99.85 % of triggers show a paddle; with the 8-16 sample dip exactly one paddle in
  practically all of them (P3 23 %, P4 22 %, P5 31 %, P6 24 %). 90 % of calorimeter showers come from the
  P4 region: the calorimeter acceptance ends near the North end.
- Positrons from the Hall D pair spectrometer arm at 1.662 T, nominally 3-6 GeV, no per-event energy (the
  August paddles were not calibrated against the PS counters). The gain-free SFIL5 and Baby BCAL double ratios
  (section 7, corrected 2026-10-23) show the shower depth rising P4 < P3 < P6 < P5, so the energy rises from North to
  South as in the May map; its size per paddle is unknown.

- **Trigger paddles as pulse-height devices (2026-10-24, `calibration/paddles/paddle_gains.{txt,png}`, `paddle_waveforms.png`,
  `paddle_unbiased.png`).** On the CAEN (channels 12-15, DRS4 cell-corrected) the paddle pulses are 10 ns wide (charge/peak
  11 ns) and reach the ADC floor (3750 counts) in 2 / 2 / 9 / 12 % of single-paddle triggers for P3 / P4 / P5 / P6. Each
  spectrum has a sharp threshold edge, at 1325 / 1275 / 1925 / 2150 counts, with its most-probable value only 10 % above the
  edge (modes 1460 / 1400 / 2090 / 2370): the discriminators sit on the MIP peak, so the trigger efficiency is paddle
  dependent and the true MPV is at or below the edge. 14-24 % of single-paddle triggers have the paddle pulse below 0.9 x its
  edge at the same time as the main population (threshold walk or a different signal copy at the discriminator). Edges and
  modes scale together, P3 : P4 : P5 : P6 = 1 : 0.96 : 1.45 : 1.62, so either the P5/P6 chains have 1.5 x the gain of P3/P4
  with thresholds set at the same point of the MIP spectrum, or the CAEN copies are attenuated differently; the two cannot be
  separated without a threshold scan or the discriminator settings. Triggers formed by another paddle do not give an unbiased
  MIP spectrum for a paddle: those coincident pulses are small (mode 175 counts) and spread in time (16 % more than 60 ns
  early), i.e. backsplash and noise, not a second crossing. On KCU0 (ASIC1, ch 84 / 86 / 81 / 97) the paddles are negative
  pulses on a 77 ADC pedestal that clip at ADC 0 in 28-38 % of triggers, carry no TOT word and a TOA in 83-93 %: timing only.

## 5. What the KCU0 waveforms are (`method.md` 4.9, 4.10; `kcu0_sfil_plots.C`, `kcu0_sfil_fit.C`)

- Pulse: rise over two samples (5-7), peak 7-9, long tail at 40-60 % of the peak to the end of the window;
  per-channel templates in `kcu0_fit_templates.png`. SFIL1 South is distinctly shorter.
- Noise: 3-5 ADC per sample on clean channels; the max of samples 4-12 over pedestal stays below ~40 ADC in
  99 % of empty waveforms; the SFIL1 North column (ch 146/150/154/158) has 19-50 ADC noise and is masked.
- Per-trigger common baseline shift of 20-45 ADC shared by all SFIL channels (not the paddles), removed by
  the per-event pedestal (samples 1-3). Pile-up tails on 2 % of beam-row waveforms.
- Saturation at 1023: of showers with a layer sum above 300 ADC, the SFIL3 South beam-row channel saturates in 53 / 58 / 62 /
  66 % of P4 / P3 / P6 / P5 triggers (SFIL2 28-45 %, SFIL4 22-53 %), rising with the paddle energy; the earlier statement of
  20-36 % under P4 and 0-2 % under P5/P6 was wrong. The template fit recovers amplitudes through saturation up to about
  1800 ADC; above that the fit pins (pile-up at the edge). Only 30 % of showers have no saturated channel at all.
  Full run: 6.98 M fitted amplitudes from 2.91 M triggers.
- TOA: on KCU0 only SFIL3 North cell 3 and the paddles carry pulse-timed TOA words; the other channels emit
  words on pedestal samples in fixed code bands. KCU1 TOA works (self-calibrated ~190 ps per channel). The LED
  macro's 6.37 ns window rule does not apply: in 3062 the forward convention t = 25 s + TOA_ns with no wrap is
  the continuous one. TOT fires only on saturated pulses (10-bit code, decompress v < 512 ? v : 8 (v - 512) + 4).
- North - South timing from the ADC crossing: 1.6-2.0 ns per end on the summing boards (too coarse for
  position); the ps goal needs TOA thresholds set on KCU0 in a future run.

## 6. Geant4 model (`sim/README.md`)

- Parameterised SFIL1-5 Pb/SciFi, AstroPix module, Baby BCAL, paddles. **Geometry v3 (2026-10-23 evening)** from the epic
  `BeamTest` branch (official DD4hep description of the August stack): 4 cells of 19.05 mm, 505-519 mm fibres, 65.1 mm layer
  pitch, 61 mm to the Baby BCAL, paddle slabs, lead/glue radiator (1.43 X0 per layer) with 0.61 mm edge skins; the AstroPix
  v3 module in front of SFIL1 (it is in the data: 1-2 pixel clusters; the official file has it commented out). Only its exact
  position is still assumed. A fibre-to-cell bookkeeping bug of v1/v2 was fixed with v3 (`pitfalls.md`).
- v3 scan 3-6 GeV (`sim/`, 5000 e+ per point; v1 in `sim/v1_purePb/`, v2 in `sim/v2_slots17mm/`): SFIL sampling fraction
  0.102, 76 MeV visible per GeV, sigma/mu 4.6 % (3 GeV) to 3.8 % (6 GeV), 74-76 % of the energy in SFILs + Baby BCAL, 15-16 %
  escaping, shower maximum in SFIL4. Six-element shares at 4 GeV: 0.024 / 0.100 / 0.169 / 0.183 / 0.151 / 0.372, within
  0.01 of v2: the longitudinal picture did not depend on the assumed spacing and cells; the lateral one did (cell 3 now
  holds 65 % of SFIL3 against the data's 72 %, cells 1+2 24 % against 28 %, cell 4 8 % against 0.5 %).
- Data versus v3 with equal gains: more in SFIL2 and less in SFIL4 than the model at any energy; with the measured
  SFIL2/SFIL3 gain ratio SFIL2 agrees and SFIL3 (+0.03) / SFIL4 (-0.015) are what an SFIL4 South gain of 0.8 x SFIL3 would give.

- **Rescan with the two-component attenuation (seeds 41-46, tag `_v3b`, 2026-10-24):** all numbers reproduce the
  v3 run within 0.003 in the shares (4 GeV: 0.024/0.099/0.170/0.185/0.154/0.369), SF 0.1021, 76 MeV/GeV.

## 7. Calibration attempts (`sim/analysis/stack_shares_calib6_v{1,2}.{pdf,txt}`, `sfil_gain_calibration_v{1,2}.txt`, `photostat_ns.{txt,pdf}`, `stack_shares_sfil_pe_v2.pdf`; description in `docs/stack_shares_calib6.md`, method 4.13)

- **Correction (2026-10-23).** The `pad` field of `events_3062.root` is a bitmask (P3 = 1, P4 = 2, P5 = 4, P6 = 8).
  The 2026-10-22 South-gain study selected `pad == 5` and `pad == 6`, i.e. two-paddle events, so its "P5 = P6 within
  0.3 GeV" was an artifact. With single-paddle masks the gain-free double ratios [element/SFIL3]_P5 / [...]_P6 are
  1.135 (SFIL5) and 1.201 (Baby BCAL), above even the scaled May map prediction (1.081, 1.103): P5 showers are
  deeper than P6 showers. The SFIL5-over-SFIL1-4 South ratio rises monotonically from the North paddle to the South
  one (P4 3.37, P3 3.53, P6 3.73, P5 3.91), as does the Baby BCAL ratio: the beam energy increases across the
  paddles in the May-map direction. Old caches/outputs are kept as `*_WRONGPAD_2026-10-22.*`.
- Improved data points: template-fitted (saturation-recovered) South-end amplitudes summed over the four cells of
  each layer, 1.97 M showers. SFIL1-4 shares 0.025 / 0.277 / 0.377 / 0.321 against 0.063 / 0.312 / 0.403 / 0.222
  before: the old SFIL4 deficit came from the beam-row-only, saturated, sqrt(SxN) measurement. SFIL2-4 sit between
  the 3 and 5 GeV pure-Pb curves; SFIL1 is less than half the simulated share.
- Preliminary per-layer gains (layer sums, ADC per MeV of simulated South light, pure-Pb model): at 4 GeV 4.8 /
  14.1 / 12.9 / 13.0, at 5 GeV 4.5 / 12.3 / 10.7 / 10.2. SFIL2-4 agree within 3-4 % at 3-4 GeV; SFIL1 needs 0.4 of
  their gain under every hypothesis. The earlier beam-row-cell gains (20 / 33 / 23.5 / 22.8) were about twice too
  high because the simulated beam straddled two cells. SFIL5 South 212 (KCU1 ADC per MeV) and Baby BCAL South 703
  (CAEN units per MeV) at a 4 GeV anchor, P6 showers.
- Six-element figure with uncertainties (`stack_shares_calib6_v1.pdf`): SFIL1-4 with one common gain (equal gains
  assumed), SFIL5 and Baby BCAL with their own simulation-anchored gains (so those two points are a normalisation,
  bar = anchor 3-6 GeV). Systematics on the SFIL1-4 shares: shower threshold (<= 0.007), saturation recovery
  (<= 0.006), noise bias (<= 0.001), statistics (negligible); the equal-gain assumption, taken from the North/South
  spread of the same layers (N/S = 0.86 / 0.97 / 0.64 for SFIL2-4, rms of ln = 0.175), gives 0.03 per layer and
  dominates: the layer shares cannot pin the energy without a real gain calibration. Paddle spread (physics): up to
  0.02.
- Equal-gain energy from the SFIL2-4 profile: 3.2 GeV all paddles (P4 3.0 < P3 3.0 < P6 3.4 < P5 4.0) with the
  pure-Pb geometry; to be redone with the corrected (v2) geometry, which moves the shower maximum later and lowers
  these.
- **Route 3, photostatistics from the two fibre ends (data only, method 4.13).** ADC per photoelectron of the South
  beam-row channels: SFIL2 14.09 +- 0.16, SFIL3 9.47 +- 0.06, ratio 1.49 +- 0.02: equal gains is wrong by 50 % between
  these layers. Both ends of SFIL2 and SFIL3 see the same photoelectrons (the SFIL2 N/S amplitude ratio 0.85 is gain).
  SFIL4's North end is anomalous (channels 186/187 carry a constant ~130 ADC, N cell 3 = 0.35 x S cell 3), so SFIL4 and
  SFIL1 South gains remain unmeasured. The v2 model implies SFIL2/SFIL3 = 1.29-1.41, v1 1.09-1.20: the measurement
  favours v2. Shares with the SFIL2/SFIL3 ratio applied (SFIL1, SFIL4 at the SFIL3 gain): 0.027 / 0.205 / 0.415 /
  0.353 against v2 at 4 GeV 0.049 / 0.208 / 0.356 / 0.387.
- v3 implied gains (layer sums, P6 showers, ADC per MeV of South light): 4 GeV 7.8 / 22.8 / 18.0 / 14.8, 5 GeV
  7.6 / 20.9 / 15.4 / 12.0 (v2 was within 10 % of these); SFIL2/SFIL3 = 1.27-1.40 across 4-6 GeV against the measured 1.49;
  SFIL5 200 KCU1 ADC/MeV and Baby BCAL 376 CAEN units/MeV at a 4 GeV anchor. The v3 gain-free double-ratio predictions
  (P6 5 GeV / P5 6 GeV: SFIL5 1.062, Baby BCAL 1.095) are still below the data (1.135, 1.201): the P5-P6 energy step is
  larger than 1 GeV or the model's depth dependence is too weak. Cell-level gains are now meaningful: beam-row (cell 3)
  gains agree between P5 and P6 showers to 2-9 %, but the data have far less light in cells 1+2 and 4 than the model
  (cell 4: 0.5 % against 8 %), so the beam y, the cell-4 channels or the lateral model still need attention.
- **Excluding saturated waveforms (2026-10-23 night, `stack_shares_pe.py <simdir> <tag> nosat|nosat_raw`).** Triggers with no
  saturated channel are 30 % of the showers (mean South total 1539 against 2287 ADC) and are biased against large SFIL3
  pulses by construction; the same threshold (930 ADC over the implied gains) keeps 80 % of the model's events, so the data
  saturate far more often than the model expects. On that sample the raw peaks and the fits agree within 0.004 in every
  share (the fit is unbiased below saturation) and the pe-equalised shares are 0.037 / 0.200 / 0.391 / 0.373 against the
  model under the same cut 0.050 / 0.197 / 0.353 / 0.400 at 4 GeV. The data have practically no light in cell 4 of any
  layer (0-1 % against 5-11 % in the model) and 88-98 % in the beam-row cell (model 65-75 %); dropping cell 4 from the model
  too gives 0.054 / 0.205 / 0.354 / 0.387 at 4 GeV: SFIL2 and SFIL4 agree, SFIL3 stays 0.035 high, SFIL1 0.02 low.
- **Saturation removed from both sides by quantile matching (2026-10-24, `stack_shares_pe.py ... nosat`).** Per layer the data lose
  the fraction f_L of showers whose South beam-row channel sits at the ceiling (0 / 0.29 / 0.54 / 0.32 for SFIL1-4, all paddles;
  rising P4 < P3 < P6 < P5); the model drops the same top fractions of its layer light, so the selection is equivalent whatever the
  gains (a 930 ADC / gain threshold had kept 80 % of the model against 30 % of the data: the data's amplitude distribution is
  1.3 x broader than the model's for the mean-implied gains). Kept fractions now 0.29-0.34 (model) against 0.31 (data). Under the
  matched cut the residuals are the same as for the full sample: data 0.037 / 0.200 / 0.391 / 0.373 against the model at 4 GeV
  0.061 / 0.202 / 0.331 / 0.406 (SFIL3 +0.06, SFIL4 -0.03, SFIL1 -0.02). The earlier "improvement" from the unsaturated sample came
  from the model not being cut equivalently.
- **Saturation-edge gains (`sat_gain.py`, `sat_gain.txt`).** The 10-bit ceiling is an absolute ADC reference: the gain at which the
  model reproduces the data's saturation fraction, g_sat, against the gain from the mean amplitude, g_mean. Per layer g_sat/g_mean
  is 1.3 and independent of the energy hypothesis (SFIL3 shows 1.6-1.8 only because its beam-row share was taken from unsaturated
  events, which are biased low), so the two never cross: the data have more events in the upper tail than the model for the
  mean-implied gains (partial showers and the energy spread within a paddle pull the mean down). Relative gains agree between the two
  estimators: SFIL2 : SFIL3 : SFIL4 = 1.31 : 1 : 0.87 (mean) and 1.40 : 1 : 0.87 (edge) for P5, 1.31 : 1 : 0.81 and 1.34 : 1 : 0.81
  for P6; photostatistics gives SFIL2/SFIL3 = 1.49 (excess noise included). Working relative gains: 1.35 : 1 : 0.85, +-0.07.
- **SFIL5 / SFIL4 South relative gain (`sfil5_sfil4_gain.py`, `.txt`).** The fit-tree triggers join the KCU1 stream through the
  KCU0 timestamp (100 % matched, 708 k with KCU1). <A5> / <A4> = 11.3-12.2 KCU1 ADC (arrays B+C+D, raw peaks, which saturate: see the KCU1 fit) per
  KCU0 ADC (fitted layer sum), nearly paddle-independent; the model's South-light ratio of the adjacent layers is 0.77-0.87 from
  3 to 6 GeV, so g5/g4 = 13.4-15.7 KCU1 ADC per KCU0 ADC for the same light (energy hypothesis +-8 %, paddle +-3 %). Dropping showers
  with a saturated SFIL4 beam row raises the ratio to 14.1-14.8 because that cut selects small SFIL4 pulses, not because the fit is off.
- **epic branches surveyed (2026-10-24, `context/SOURCES.md`).** The FTBF Baby BCAL macros in `epic/analysis/` supply a measured
  readout model: 1100 photoelectrons per GeV before attenuation, a two-component attenuation (18 % at 61 cm, 82 % at 418 cm),
  per-depth light-guide efficiencies, a 1 MeV cell threshold and the Birks range 0.063-0.189 mm/MeV; their data were MIP-calibrated.
  The attenuation is now in `bic_sim` (it moves South-end amplitudes by about 7 % between the P4 and P5 positions, the size of the
  paddle trend in g5/g4); the scan is being rerun with it (seeds 41-46; the single-exponential v3 files are in `sim/v3_singleexp/`).
- **SFIL5 South clipped pulses: what the KCU1 waveforms allow (2026-10-24, `waveforms/kcu1_tot/`).** The single-SiPM channels
  saturate at the ceiling like the summed ones (the earlier "never saturate" was a threshold mistake; 2.1 clipped channels per
  shower trigger on average). Raw-waveform study (`kcu1_tot/` scripts, 300k triggers): the KCU1 pulse shape narrows steadily
  with amplitude from 100 ADC up (sample after the peak 0.81 -> 0.72 of the peak, three samples after 0.75 -> 0.50 between
  100 and 940 ADC), whereas the KCU0 summed channels keep a fixed shape to 800 ADC; and after a clipped stretch the next
  KCU1 sample sits at ~750 ADC whatever the charge (1 to 7 clipped samples). So the template extrapolation (median 0.70 of the
  ceiling for one clipped sample) and every ADC-tail method are invalid for clipped KCU1 pulses: the only measure is the
  time over the ceiling. The decoded TOT is exactly that: median 502 for one clipped sample, +166 units per additional 25 ns
  sample (ceiling crossing T_c = 419), present on 63 % of 1-sample and ~90 % of >= 2-sample clipped pulses (missing ones get
  the median TOT of their nsat). The recovery is A = 943 + b (TOT - T_c) with the slope b in ADC per TOT unit **not
  determinable from run 3062**: (i) a neighbour predictor within an array fails (the 16 SiPMs see different fibre bundles,
  factor 2.5 spread, position-driven bias); (ii) the fall of <A5>/<A4> with A4 (1.31 -> 0.43 normalised to 300-400 ADC,
  SFIL4 unsaturated) is reproduced by the model's own adjacent-layer relation (<L5|L4> 1.13 -> 0.45 at g4 = 15 ADC/MeV): the
  "factor 4.5 nonlinearity" claimed earlier is shower physics (L5 rises only 0.18 per unit L4 at fixed energy), not a gain
  effect, and is insensitive to b; (iii) the width of the per-paddle SFIL5 sum (q90/q50 1.33-1.50 at b = 0) already matches
  the model's L5 width (1.31-1.45), so b = 0 unless the beam energy spread within a paddle is 20-30 % (b = 0.40-0.76). Adopted
  **b = 0.40 with the band 0-0.76** (the toe slope below the ceiling is 0.33): mean SFIL5 sums 7539 / 8236 / 8879 ADC at
  b = 0 / 0.40 / 0.76 (ceiling-only sums are 4 % below the old fit+raw sums). A proper slope needs the LED runs 3179-3190 or
  the HGCROC charge injection. KCU0 keeps the template recovery (its TOT is absent on saturated pulses; its shape is stable).
- **SFIL5 / SFIL4 relative gain with the recovered sums (`sfil5_sfil4_gain_rec.{txt,png}`):** <A5>/<A4> = 11.0-11.4 across the four
  paddles (394k showers with KCU1), g5/g4 = 14.0 at a 4 GeV anchor (13.7-14.3 by paddle; 14.6 / 13.0 at 3 / 6 GeV; +8 % at b = 0.76,
  -8 % at b = 0). The model enters only through the adjacent-layer light ratio 0.77-0.87.
- **KCU1 ADC per photoelectron from the four live SFIL5 North SiPMs (method 4.13):** 24-34 ADC/pe (x F) for ch 54/60/65/66, i.e. the
  KCU1 single-SiPM chain gives 2.5-3.5 x the ADC per photoelectron of the KCU0 summed channels (9.5-14). The South SFIL5 sum cannot be
  calibrated against them (bare SiPMs sample a third of the fibre face; correlation 0.2-0.3).
- **SFIL5 in the share plot (`stack_shares_sfil5_pe_rec_v3b.pdf`, recovered sums, v3b scan).** With g5/g4 = 14.0 the five South
  shares are 0.021 / 0.160 / 0.322 / 0.275 / 0.221 against the 4 GeV model 0.041 / 0.166 / 0.273 / 0.288 / 0.231; SFIL5 sits on the
  model at the anchor by construction (its bar is the 3-6 GeV range of the model light ratio), SFIL3 keeps its +0.05 excess,
  everything else is within 0.02. The b band moves the SFIL5 point by less than its bar because the anchor renormalises it.
- SFIL5 South and the Baby BCAL cannot be photoelectron-calibrated from run 3062: intra-array SiPM pairs share
  only 20 % of their variance (direct coupling, no light guide) and no dark-count peak is resolvable on KCU1, KCU0 or
  CAEN (method 4.13). They stay simulation-anchored in `stack_shares_calib6_pe_v2.pdf`.
- **Cosmic MIP data exist off-site (2026-10-24):** H. Klest's plot of run 3089 segment 0000 (`context/`) shows North MIP peaks of
  15-30 ADC with a South coincidence, on a noise core of 8-10 ADC; the converted file `kcu0_run3089_0000_baseline4_paddles_waveforms.root`
  is on his machine. August run 2804 was an H2GCROC calibration (injection) run for KCU0 ASIC0. Neither is on this disk.
- Not yet done: a data-only cross-calibration of the three readouts (KCU0, KCU1, CAEN); the cosmic runs 3088-3172
  (MIP, all readouts) and the LED runs 3179-3190 are the way to replace every simulation-anchored gain; the SFIL4
  North and SFIL1 South anomalies need the hookup table.

- **Baby BCAL calibration: what exists (2026-10-24).** Nothing usable as a gain table: (i) the 2024 FTBF
  constants (MIP method and electron method, `context/arXiv_2504.03079`) are not in any repository here; the paper
  prints the per-channel MIP peaks (South 12, 11, 6, 5 | 27, 23, 24, 23 | 25, 28, 85, 24 | 35, 38, 35, 35 fADC
  units for channels 0-15, about 10 units and 9-12 pe per array) for a JLab fADC250, not for the CAEN, and
  flags S2 and S3 (no MIP peak) and S10 (faulty amplifier); (ii) the epic `BabyBCal` macros carry only the
  digitization model (1100 pe/GeV, per-depth light-guide efficiencies 1.15-0.74, Birks); (iii) our own
  Baby BCAL gains are simulation-anchored totals (316-376 CAEN units/MeV). The August South side reads 12 of the 16
  FTBF channels (A0 A1 A2? A3 A4 | A5 A6 A7 | A8 A9 A11 | A14; `subsystems.md` 4, labels from the run-log sheet), the North side is silent, so
  neither the N-S geometric mean nor N-S photostatistics are available, and the two dominant channels A6/A7 are
  different columns (no shared light).
- **Baby BCAL after the DRS4 cell-pedestal correction (2026-10-24, `calibration/caen_calib/bbcal_classes.{txt,png}`,
  segments 1-2 of the CAEN stream so far).** With the per-cell pedestals (`pitfalls.md`) the baseline rms is 2
  counts in the pre-pulse region (9 further into the window). Two findings: **CAEN channels 2 and 3 (A3, A4) are
  dead**: zero shower charge, no correlation with the uncorrected values, and the same 9 % "above 5 sigma" rate as
  the noise; their earlier 1.4-1.5 % shares were the cell artefact. Corrected shower shares, all paddles, with the
  run-log sheet's channel labels (CAEN 4-10 = A2(?) A5 A6 A7 A8 A9 A14; `context/SOURCES.md`): A5 0.270, A6 0.281, A9 0.103,
  A1 0.085, A2 0.085, A7 0.071, A11 0.033, A8 0.028, A0 0.024, A14 0.019, stable within 1 point across paddles; every row
  peaks between columns 1 and 2, so the earlier "row 2 anomaly" was the label shift. **Run 3062 contains no usable MIP or cosmic sample for the Baby BCAL**:
  there is no 10 Hz pulser in the run (`method.md`), the pair-spectrometer arm delivers no muons, the random cosmic
  coincidence budget in the 256 ns CAEN window is ~8 events for the whole run, and in the 38k SFIL-quiet triggers
  (paddle fired, all four SFIL layers below 100 ADC) the live Baby BCAL channels are above threshold only 0.2-2 %
  more often than the dead ones; a one-column track pattern occurs 3 times in 820k triggers, and SFIL-quiet triggers
  with a Baby BCAL sum above 20k units (direct hits beside the SFIL stack) number 19, not the 9.9k of the uncorrected
  data. Routes left for a per-channel Baby BCAL calibration: the cosmic runs 3088-3172 and LED runs 3179-3190
  (raw files not on disk), or gains anchored to the v3 model by regrouping its 10 depth segments x 4 cells into the
  paper's rows (segment 0 | 1-2 | 3-5 | 6-9) and columns.

- **Interim Baby BCAL per-channel gains anchored to the v3 model (2026-10-24, `sim/analysis/bbcal_gains.py`,
  `bbcal_gains_yfine.{txt,png,npz}`).** The model's `bb_ls` (10 depth segments x 4 cells, South light) is regrouped into the
  FTBF rows (segment 0 | 1-2 | 3-5 | 6-9) and columns (cells); the vertical position of the Baby BCAL block (new
  `/bic/geom/bbcalYOffset`, condor scan -37.5 to +37.5 mm in 12.5 mm steps, then -7.5 / -5 / -2.5 mm, 3-6 GeV at -5 mm) is
  chosen so that the within-row patterns of rows 0 and 1 are reproduced with equal gain per array. **The block sits 5.0 mm
  below nominal**: at -5 mm the model gives row 0 A0 : A1 : A2 = 0.123 : 0.432 : 0.446 and row 1 A5 : A6 : A7 = 0.441 :
  0.455 : 0.104 against the data 0.123 : 0.440 : 0.437 and 0.437 : 0.448 : 0.115 (pattern chi2 0.03; 0.34 and 0.23 at -7.5 and
  -2.5 mm). Gains at a 4 GeV anchor in CAEN charge units per MeV of South light (per summed array in brackets): A0 226, A1 218,
  A2 223 (row 0, one array each); A5 510 (255), A6 537 (269), A7 572 (286); A8 356 (119), A9 368 (123), A11 418 (139); A14 289
  (72). The 3 -> 6 GeV anchor range scales them by +40 % / -33 %. Within a row the per-array gains now agree to 2 % (row 0) and
  6 % (rows 1 and 2), so the channel-to-channel equalisation is good to that level. What remains is the **fall with depth: 220
  (row 0), 270 (row 1), 127 (row 2), 72 (row 3) per array**, a factor 3 between the front and the back rows that the FTBF MIP
  peaks (10 units per array at every depth) do not show. Either the deep rows receive less light than the model's tail
  behind SFIL5 says, or the 3- and 4-array sums are not one-to-one; the cosmic data decide. Rows 0-1 gains are usable now;
  rows 2-3 carry the model's longitudinal shape as a systematic. Per paddle the gain-corrected data over the model are
  0.75 / 0.89 / 1.04 / 1.24 (P4 / P3 / P6 / P5), uniform over channels: the beam-energy ordering again.

- **AstroPix ToT calibration along the FTBF paper's method (2026-10-24, `astropix/astro_tot_calib.{txt,png}`,
  `astro_tot_mpv_map.png`, `astro_tot_pixels.npz`).** The paper (arXiv:2602.06084) pairs row and column records (|dToA| <= 1,
  |dToT|/ToT < 10 %), calibrates every pixel with 241Am / 133Ba photopeaks (second-order polynomial ToT -> keV), fits the
  120 GeV proton ToT spectrum per pixel with a Landau x Gauss and finds a MIP MPV of 34.0 +- 1.9 keV (128 um effective
  depletion). For run 3062: the same pairing gives 543k pixel hits (2.59 M clean records, upper ToT nibble 0). **Chips 5-8
  carry no ToT** (38-79 % of their column records read 0, medians 3-7), so only chips 0-4 can be calibrated. Chip-level
  Landau x Gauss MPVs in the beam band: 676 / 484 / 585 / 724 / 331 ToT units for chips 0-4 (in ToT clock units of **10 ns, settled by the firmware** (`clk_sample = clk_100` in the A-STEP top module of the
  external-clock build; Bobae Kim's CSV carries a 5 ns software default), so the MPVs are 6.8 / 4.8 / 5.9 / 7.2 / 3.3 us; the
  August configuration is `module8.yml` (`context/astropix_module8/`): threshold 202 mV on all nine chips, columns 0-2 and 18 hot
  pixels masked, injection off, identical settings on chips 5-8 (their missing ToT is not a configuration choice); against the representative v3 per-pixel calibration in `A-STEPdee` (6.3 keV/us median) they read 43 / 31 /
  37 / 46 / 21 keV, i.e. within +-35 % of a 34 keV MIP, to be settled by the module8 yaml; the module was configured with `set.py -y module8 -c 9`, the yaml is on the DAQ machine). Per pixel
  (1411 pixels with >= 150 hits, Landau fits, 2-3 % fit error) the MPV spread within a chip is 23-37 % rms, the raw
  pixel-to-pixel ToT gain variation the paper removes with its source curves. No source data exist for this module here, so
  the calibration is relative (per-pixel MPV map) with the absolute scale anchored to the MIP: the Landau-Vavilov MPV for
  3-6 GeV positrons in 128 um Si is 34.0 keV, the same as for the paper's protons (Fermi plateau), so E = 34 keV x
  ToT / MPV_pixel near the MIP (+-6 % from the paper's value, more where the depletion of these chips differs). Per chip that
  is 19.9 / 14.2 / 17.2 / 21.3 / 9.7 ToT units per keV. A proper curve needs the module's bench 241Am / 133Ba scans from the
  AstroPix group.

- **Attenuation and North/South gain ratios from the AstroPix x, model-free (2026-10-24, `calibration/attenuation/`).**
  45.6k showers carry an AstroPix cluster of purity >= 0.3 (chips 3, 5, 6, 8; stack-frame x = AstroPix x - 88.75 mm, from
  -32 to +88 mm; P4 has no matched chip). Beam cell, both ends unsaturated, ln(A_N/A_S) median per 10 mm against x:
  SFIL2 slope -0.00150/mm -> lambda_eff 1329 mm, g_N/g_S = 0.85; SFIL3 -0.00122/mm -> 1646 mm, g_N/g_S = 1.03; SFIL4
  -0.00323/mm -> 619 mm, g_N/g_S = 0.29 (SFIL4 North cell 3, ch 144, is on the v2-board side with the mixed arrays, so its
  slope is not an attenuation). The epic two-component attenuation predicts a slope of -0.00085/mm (lambda_eff 2341 mm): the
  data attenuate 1.4-1.8 times faster over the central 120 mm. Profile residuals 0.01-0.04.
- **Beam energy versus position from the geometric mean (same script).** sqrt(A_N A_S) cancels an exponential attenuation;
  normalised at x = 0 the SFIL2-4 geometric-mean layer sums (fits, 42.5k showers) are 1.000 at -30, 1.002 at +30 and 1.125 at
  +80 mm; the unsaturated beam cells give 1.02 / 0.99 / 1.05-1.15. So E(P5)/E(P6) = 1.12 (1.05 from the earlier share fits)
  and the energy is flat within 3 % between P3 and P6. P4 (x = -80) is not covered by the matched chips.
- **Position along the fibres from the North - South time (2026-10-24, `timing/toa_position/`).** The KCU0 TOA words on
  SFIL2 cell 3 (S ch 29 / N ch 173, both valid on 63 % of pulses) and SFIL4 cell 3 (S 0 / N 144, 21 %) are not usable for
  timing even after a per-channel code-to-time calibration against the ADC crossing: the code occupancy is not uniform, the
  code-time curves have plateaus, the TOA - crossing residual is 4-8 ns rms and the N - S difference is double-peaked with a
  7 ns rms. **The ADC 50 % crossing works:** the N - S crossing difference moves by 0.0128 ns/mm (SFIL2) and 0.0126 ns/mm
  (SFIL4) across the four paddles, i.e. v_eff = 2/slope = 157 and 159 mm/ns, the fibre value, with an rms of 1.8-2.6 ns per
  shower and layer (position rms about 140 mm per layer at 25 ns sampling). Fixed N - S offsets 0.66 ns (SFIL2) and 3.25 ns
  (SFIL4) at x = 0.
- **SFIL2 cell 3 TOA without the ADC timing (2026-10-28, `timing/toa_clean/toa_clean{,_figs}.py`, supersedes the 2026-10-27
  entry below).** The 10-27 code-to-time table was the ADC crossing phase per code and only words within 6 ns of the ADC crossing
  were kept, so that TOA was tuned to the ADC (its 2.1-2.5 ns N - S and the TOA/ADC/common split are not valid). Now: t = 25 s +
  T(code), s = sample of the word (k, else k+1; the ADC only finds the pulse), T from the code distribution assuming a uniform
  phase, even triggers for the calibration, odd for results, no time window. Noise words: 14 % (S) / 37 % (N) of samples before
  a pulse, on the same codes as pulse words, so no cut removes them. 104k odd-trigger showers with a word at both ends (59 %).
  corr(t_N, t_S) 0.85; 71 % within 5 ns of the diagonal. N - S core sigma 4.1-4.5 ns from codes only (116 mm/ns), 2.8-3.2 ns
  after an amplitude/baseline regression (ADC amplitudes and pedestals only; 160 mm/ns), core fraction 0.89-0.92; ADC crossing
  1.7-2.0 ns (146 mm/ns). sigma_x TOA 221-253 mm, ADC 125-146 mm; P4 vs P5 separated by 0.71 sigma (TOA), 1.24 sigma (ADC);
  against AstroPix 232 vs 128 mm. Within-paddle TOA/ADC residual correlation 0.38. Conclusion unchanged: the TOA times the
  shower but cannot place it; use the ADC crossing.
- **Dual-sided readout study, repeat of the 2017 sPHENIX EmCal note (2026-10-28, `dual_sided/`, note
  `note/dual_sided_note.pdf`).** SFIL2/3 cell 3, one paddle, no clipping (35.5 % of showers clipped): asymmetry (S-N)/(S+N) changes by
  0.116 / 0.094 over 160 mm (2017: 0.014), lambda 1382 / 1702 mm; N-S P5-P4 2.0-2.2 ns (2017 2.37 ns), v_eff ~150 mm/ns; per-shower
  width CFD (f 0.8, d 25 ns) 2.1-2.3 ns, ADC 50 % 1.7-2.0 ns, template 2.3-2.7 ns on SFIL2 (2017 ~9.3 ns rms).
- **Good-TOA selection, ADC-free (2026-10-28, `timing/toa_clean/toa_good{,_fig}.py`).** Cuts from the TOA words only, both ends:
  no word before the pulse (samples 1..k-1), not two words on k and k+1, code in a 16-code bin where pulse words outnumber the
  expected noise words (no-pulse code shape x pre-pulse word rate; keeps all South bands, 56 % of North pulse words). Kept
  29 / 27 / 17 % of showers, the same fraction per paddle and the same amplitudes. Final: TOA N - S 1.75-1.89 ns (166 mm/ns, core
  fraction 0.96), sigma_x 145-156 mm, P4-P5 1.10 sigma; ADC on the same showers 1.64-1.91 ns, 119-138 mm, 1.31 sigma. The noise
  words are what limit the TOA; clean words come within ~20 % of the ADC.
- **SFIL2 cell 3 TOA redone (2026-10-27, `timing/toa_sfil2/`, method 4.18).** With the pulse word taken on the crossing sample or
  the next, the code unwrapped circularly and the per-trigger baseline offset corrected (-0.05 ns/ADC), the KCU0 TOA is a working
  TDC (monotonic code, ~25 ps per code nominal), but its per-pulse error is 1.6-1.8 ns (in N - S) against 0.8-1.2 ns for the ADC
  50 % crossing; N - S sigma TOA 2.1-2.5 ns, ADC 1.6-2.3 ns, combined = ADC. A common 1.4-1.9 ns (light timing) is the floor:
  ~100 mm per shower for any readout of this channel. Sub-ns N - S needs single-SiPM channels (KCU1, ~190 ps) or a TOA threshold
  set above the baseline wander in the next run.
- **Rate and pile-up dependence: none at the 10^-3 level (2026-10-24, `calibration/rate/rate_dependence.{txt,png}`).** For 1.97 M
  showers, the SFIL1-4 South sum normalised to its paddle mean is 1.000 +- 0.001 in every bin of the time since the previous
  trigger from 0.10 ms (the busy) to 5 ms, identical to the control against the time to the next trigger, 1.000 +- 0.001 in
  every rate bin from 1 to 5 kHz (the 0.97 +- 0.02 point below 1 kHz is 354 showers at beam-off edges), and independent of
  the previous trigger's own size (1.001 +- 0.001 after showers above 4000 ADC within 0.5 ms). No layer deviates by more
  than 0.5 % where the statistics allow. The gains above therefore carry no rate systematic from this run; the 2 % pile-up
  seen in the waveforms is handled by the pile-up flag, not by a shift of the response.
- **Per-event position along the fibres for every shower (2026-10-24, `timing/toa_position/event_x_walk.{py,txt,npz,png}`).**
  Two handles on the beam cell of each layer. Timing: the template-fit t0 of both ends (`kcu0_sfil_fit_full.root`, 20 %
  crossing of the aligned template with clipped samples masked), dt = t0_N - t0_S calibrated on the four paddle positions:
  slopes 0.0124 / 0.0133 / 0.0146 ns/mm for SFIL2 / 3 / 4, i.e. v_eff = 162 / 150 / 137 mm/ns, with no amplitude walk left to
  correct (an explicit ln(A_N/A_S) term in the regression is degenerate with x through the attenuation and is not applied;
  clipped and unclipped pulses give the same rms). Amplitude: x_amp = -(lambda_eff/2)(ln(A_N/A_S) - ln g_N/g_S) with the measured
  1329 / 1646 mm and 0.85 / 1.03. Per layer and shower the rms against the paddle position is 199 / 151 / 258 mm for the timing
  (SFIL2 / 3 / 4) and 138 / 150 mm for the amplitude (SFIL2 / 3), medians within +-35 mm. The inverse-variance combination
  places 1.42 M showers at 132 / 137 / 149 / 175 mm rms for P4 / P3 / P6 / P5 with median residuals of -16 / +2 / -27 / -25 mm,
  and 145 mm against the AstroPix x. The five estimators do not average down (time-only 171, amplitude-only 143, combined
  132-175 mm), so the 140 mm is the per-shower floor of light sharing at 25 ns sampling, usable for statistical x
  dependences across the whole module but not as a per-event tag at the 10 mm level.

- **Cell-to-cell equalisation within a layer (2026-10-24, `calibration/cells/cell_equalisation.{py,txt,png}`).** From the raw
  trees (cells by name, first 800k triggers, raw peak minus the empty-trigger bias, beam cell 100-900 ADC): the off-beam
  cells carry, relative to the beam cell, South cells 1+2 = 0.067 / 0.034 / 0.284 / 0.098 and cell 4 = 0.027 / 0.010 / 0.008 /
  0.006 for SFIL1-4; North cells 1+2 = 0.22 (SFIL3) and 0.55 (SFIL4, the channels 187/186 with the constant pulse, not a
  shower share), SFIL2 North cells 1-2 read 0.003 (no response) and SFIL1 North is bias-dominated (its column noise). The v3
  model at the nominal beam height (y = -5 mm, 5 mm inside cell 3) puts 0.23 / 0.32 / 0.44 / 0.58 of the beam cell into cells
  1+2 and 0.02-0.20 into cell 4, so the data see 3-10 times less lateral light than the model in every off-beam cell, and
  the pattern across layers (SFIL3 highest, SFIL2 and SFIL4 lowest) is not the monotonic widening of a shower. Cells 1 and 2
  of the v2 boards are one summed pair (event-table docstring, Bobae's deck), so only the pair is a physical channel. A
  per-cell gain cannot be separated from the beam height with the data alone. The beam-height scan of the model (y = -14,
  -9.5, -5, -2, 0, +2 mm; cell 3 spans -19.05 to 0) gives for the cells-1+2 pair over the beam cell 0.06 / 0.13 / 0.23 / 0.35
  (SFIL1-4) with the beam at the cell centre (-9.5 mm), 0.23 / 0.32 / 0.44 / 0.58 at -5 mm and about 1 at the cell edge; cell 4
  is 0.06-0.35 at the centre and never below 0.02, against 0.006-0.027 in the data, so cell 4 is dead rather than unlit and
  cannot locate the beam. Reading the pair ratios with the beam at the cell centre: SFIL1 S 0.067 / 0.059 = 1.1, SFIL3 S
  0.284 / 0.229 = 1.2 and SFIL3 N 0.222 / 0.229 = 1.0 are consistent with equal gains, while **SFIL2 S 0.034 / 0.126 = 0.27 and
  SFIL4 S 0.098 / 0.353 = 0.28 respond at about 0.3 of their beam cell** (0.11-0.17 if the beam is 5 mm inside the cell, 0.5 if
  at -14 mm). So: the beam is near the centre of cell 3, the summed pair channels of SFIL1 and SFIL3 are on the beam-cell
  scale, those of SFIL2 and SFIL4 South are a factor 3-4 low (a relative gain, or the v2-board mixing), and SFIL2 North cells
  1-2 are off. Cell 4 stays empty (dropped from the model in comparisons).
## 10. Physics performance at the present calibration (2026-10-26)

- **Imaging layer (`astropix/astro_imaging.{txt,png}`).** 2.21 M showers; 5.7 % have any AstroPix cluster matched,
  3.4 % with purity >= 0.3 and 1.6 % with purity >= 0.7 (the statistical matching through the on-chip timestamp, `method.md`
  4.3). Matched clusters are single pixels in 67 % and two pixels in 31 % of cases (median 1, 90th percentile 2): a 3-6 GeV
  positron in front of SFIL1 is one hit, no charge sharing to speak of at 0.5 mm pitch. The beam height on the chips is
  2.8 mm rms (rows 14.5-16.5 median, rms 5.6-5.8 rows), the same on chips 3, 5, 6 and 8. Matched fractions per paddle:
  P4 1.7 % (chip 0 is under it but its clusters never pass the purity cut), P3 8.0 % (chip 3), P6 10.3 % (chips 5-6), P5 3.0 %
  (chip 8). Hit records over the run: chip 0 625k, chips 1-8 240-360k (80-120 row records per second); chip 0's excess is not
  beam-correlated.

- **Calorimeter (`sim/analysis/performance.{py,txt,png,npz}`).** 99.9k showers with SFIL1-4 (KCU0 fits), SFIL5 (rebuilt)
  and the Baby BCAL (cell-corrected, per-channel gains), each South amplitude corrected to x = 0 with the measured attenuation
  (lambda_eff 1329 / 1646 mm for SFIL2 / 3, 1500 mm assumed elsewhere; factors 1.055 to 0.948 across the paddles) and converted
  to visible MeV with the 4 GeV model-anchored gains. Paddle energies: P6 = 4.00 GeV is the anchor, P5 = 4.50 and P3 = 4.00 from
  the geometric-mean curve, P4 = 3.45 from the SFIL2-5 profile relative to P3 (the profile alone gives P5 / P3 / P4 over P6 =
  1.19 / 0.81 / 0.70, the geometric mean 1.12 / 1.00 / -, the earlier six-share fits 1.05 / 0.95 / 0.90: the energy per paddle is
  known to 10-20 %, and that limits every statement below).
  - **Linearity:** <E_vis>/E = 62.0 / 60.8 / 64.1 / 64.9 MeV per GeV for P4 / P3 / P6 / P5 against the model's 64.5: linear
    within +-3 % over the covered 3.5-4.5 GeV, after the attenuation correction (without it the South sums spread by 9 %).
  - **Energy resolution:** core sigma/mu of the total = 13.2 / 12.2 / 10.9 / 10.5 % (P4 / P3 / P6 / P5) against the
    monoenergetic model's 5.0 / 4.8 / 4.8 / 4.6 %. If the model resolution holds, the beam energy spread within a paddle is 9-12 %;
    the number is an upper limit on the detector resolution since the per-event SFIL5 TOT slope and the Baby BCAL depth gains
    also broaden it. SFIL1-4 alone: 29.9 / 26.5 / 20.4 / 18.8 % against the model's 17.7 % at 4 GeV, i.e. the four front layers
    reproduce the model's (containment-limited) resolution at the South paddles and are 50 % broader at the North ones; adding
    SFIL5 and the Baby BCAL halves the width, as containment should.
  - **Longitudinal profile:** the SFIL2-5 shares at the anchored energies agree with the model within 0.01-0.03 per element
    (SFIL3 still +0.02-0.03 high, SFIL4 and SFIL5 within 0.01), and the shower-maximum shift with energy is visible across the
    paddles. **Containment:** the Baby BCAL carries 0.31 / 0.32 / 0.34 / 0.36 of E_vis against 0.38-0.41 in the model, 15-20 % low,
    the same deficit as the fall of the Baby BCAL per-array gains with depth; the model puts 15-16 % of the energy behind the
    Baby BCAL.
  - **Position:** 140 mm per shower along the fibres from light sharing (section 7), 2.8 mm beam height from the AstroPix.

## 8. Plan versus status

| step | status |
| --- | --- |
| 1 event synchronization | KCU0-KCU1 and CAEN final; AstroPix statistical, per-event purity ~80 % on three chips; completeness limited by DAQ loss |
| 2 TOA/TOT ps timing N-S | not possible in 3062 (KCU0 TOA not pulse-timed); recipe and self-calibration ready for a run with thresholds |
| 3 waveform calibration (pedestal, saturation) | per-event pedestal, pile-up flag and template fits done; saturated KCU0 amplitudes recovered by template to ~1800 ADC; KCU1 clipped pulses rebuilt from the TOT (slope undetermined, band given); an LED or charge-injection calibration is needed for the slope |
| 4 AstroPix hot-pixel scan | not started (tools: `astropix_3062.root` hits tree, `bic_streams.astro_phase`) |
| 5 sampling fraction and energy scale | simulation gives SF 0.089 and the layer expectations; absolute scale blocked by the unknown beam energy; preliminary MC-anchored gains exist |
| 6 shower reconstruction, SFIL position, leakage | layer shares and N/S asymmetry done; AstroPix clustering exists; full reconstruction not started |
| 7 physics performance | not started |

## 9. Open questions, in order of leverage

1. Pair spectrometer momentum map and the module survey for August (energy per position).
2. Setup metrology: cell heights, fibre length, layer spacing, AstroPix slot (replaces the ASSUMED values).
3. KCU0 TOA configuration of run 3062 (why most channels emit untimed words) and thresholds for the next run.
4. Why ASIC3 (Baby BCAL North) and the SFIL5 North board were silent.
5. DAQ: UDP loss and the 100 us busy (`completeness_plan.md`).
6. Trigger definitions of the cosmic runs (3088-3172) and of the low-current runs 3066-3069: neither is the paddle OR.
7. Where are cosmic run 3089, the pedestal runs 3173-3178 and the LED runs 3179-3190 (not on `/gpfs01/eic/data1/BIC_Beamtest/`).

## 11. New runs on /gpfs01 (2026-10-28, `new_runs.md`)

- Beam runs 3061, 3063-3069 and cosmic runs 3088, 3090, 3091 are on `/gpfs01/eic/data1/BIC_Beamtest/`, same RCDAQ format
  and packets as 3062; `kcu_to_root` converts them unchanged (300-packet tests, 0 ADC mismatches).
- **3064** (12 h, ~69 M KCU0 / 28.5 M KCU1 / 28.7 M CAEN triggers) looks like 3062 in the sample (SFIL2 cell 3 and paddle
  fractions within a few %); AstroPix weak (1.51 M hits in 12 h vs 2.86 M in 3062's 25 min); its end is a low-current TAC period. 3063 extends 3062 by 6 %.
- **3066-3069** (low current, ~100 k triggers): beam-cell showers but a paddle on only 4-10 % of triggers (trigger unknown).
- **3088** (15 h cosmics, 1.04 M KCU0 triggers) and 3091 (20 min, AstroPix on): late pulses (samples 12-15), ~20 ADC per cell
  after per-ASIC common-mode removal, MIP peak not yet shown; **SFIL5 North summed boards now on KCU0 ch 153 / 156 / 163 / 167**
  (hookup change Aug 14 21:55 CT); SFIL5 South single SiPMs see nothing in cosmics.
- Still missing: cosmics 3089 (47 h) and 3092-3172, pedestal 3173-3178, LED 3179-3190, HGCROC injection run 2804, beam run 3060.

## 12. To-do (2026-10-28)

Data and tools for the new runs:
1. ~~Parameterise the condor converters~~ done 2026-10-28 (`run_to_root.{sh,job}`, `run_merge`, `run_watch.sh`); all new runs submitted (cluster 419386), merges follow automatically; check `output/root_logs/run_watch_419386.log`.
2. ~~Run-dependent KCU0 channel map~~ done: runs >= 3085 get the `sfil5n` tree (ch 153 / 156 / 163 / 167).
3. Per-run synchronization (KCU0-KCU1 offset, CAEN ratio/offset, AstroPix phase) and trigger table for 3063, 3064, 3088.
   Template pulse fits of all beam runs submitted 2026-10-29 (cluster 419662, `new_runs.md`).
4. ~~Decide full conversion vs skim for 3064~~ (converted in full, merged in parts of 10 segments; `new_runs.md`); find the TAC
   period from the rate against time.
5. Ask for the trigger definitions of the cosmic runs and of 3066-3069; ask for 3089, the LED and the pedestal runs.

Calibration chain (status deck order):
6. Step 3/5: MIP scale for SFIL1-4 from 3088 (2026-10-30: the pulses are at the beam timing, not late; two-end muon tag; jobs
   running, section 13.4); Baby BCAL South from the 3088 CAEN stream (cell correction running, 13.5); first look at SFIL5 North
   (KCU0) in cosmics (MPV 11-14 ADC from 10 k triggers, 13.4).
7. Step 3: off-beam cell gains and ~~gain stability over 12 h with 3064~~ (stability done 2026-10-30, section 13.3; off-beam
   cells still to do).
8. Step 4: SFIL5 TOT slope with the 3064 KCU1 statistics (LED runs still the real answer; nothing in 3064 changes the 3062
   argument, the slope stays 0.40 with the band 0-0.76 until an LED or injection run).
9. Step 5: Baby BCAL depth trend with 3064 CAEN statistics and the 3088 cosmics (CAEN cell correction of 3064 / 3088 / 3091
   submitted 2026-10-30, 13.5; the 3064 CAEN alignment is being redone, 13.1).
10. Step 6: joint fit of the paddle energies with the spectrometer optics (needs the PS momentum map).

Analysis follow-ups:
11. Dual-sided and good-TOA studies on 3064 (18 x the showers); the good-TOA cut keeps 17 %, so ~10 M showers.
12. ~~KCU1 single-SiPM SFIL5 North - South timing (sub-ns test)~~ not possible with the August data: in the beam runs the SFIL5
    North board on KCU1 was dead, and in the cosmic runs (SFIL5 North moved to KCU0) the SFIL5 South single SiPMs on KCU1 see nothing.

## 13. Threads, state and next steps (2026-10-30)

One paragraph per thread; the numbers are in the files named. Condor outputs of the day under `/sphenix/tg/tg01/jets/tmengel/eic/`:
`{kcu0,kcu1}/fit_runs/<tag>/` (template fits of every beam run), `{kcu0,kcu1}/stability/`, `kcu0/cosmics/`, `caen/calib/`, `sim/v4/`.

**13.1 Synchronization per run (`sync/run_sync.py`, `output/sync_runs/run<RUN>/sync.txt`).** KCU0-KCU1 is final for every run
(92.7 % of shared triggers at 0 ticks, 7.3 % at -4; exact offsets 357099 / 306738 / 198368 / -203777 / 306604 ticks for 3063 /
3064 / 3066 / 3088 / 3091; in the cosmic runs the -4 class is gone and 11 % sit elsewhere). CAEN-KCU0 works for the short runs
(3063 66 %, 3066-3069 99.5 %, 3091 99.9 % within 100 ns, MAD 6.4-6.6 ns, clock ratio 1.3653407 +- 0.0000001 in all of them) and
had **failed for 3064 (0.6 % matched) and 3088 (9.8 %)**. A direct scan of 3064 (23 slices of 2000 CAEN triggers, +-10 s,
scratch `caen_diag.py`) shows the offset is simply constant, -0.1414 s at 0.7 h to -0.1399 s at 12.3 h with the common clock
ratio (a 30 ppb residual drift), with 800-1900 of the 2000 triggers of every beam slice matched; `run_sync.py` had anchored on
the empty pre-beam slice (its score n / b with b = 1) and then walked a +-2 ms window from a wrong offset. Fixed 2026-10-30:
the anchor is the most significant of 15 slices over the middle 90 % of the run (score (n - b) / sqrt(b + 1), >= 50 pairs,
+-10 s search); 3064 and 3088 resubmitted (`run_sync_list_caenfix.txt`). The host-clock epochs wander by 1.3 s (KCU) / 1.6 s
(CAEN) over the 12 h (bursty host stamps on all three systems), so they fix the coarse offset only to about 1 s, which the
+-10 s search covers. AstroPix: 3066 (218 sigma, phase 2) and 3064 (21.9 sigma wide scan) lock; 3063 / 3091 are weak (4-6
sigma). The decoding / sync note (`sync/note/sync_note.tex`, `run_sync_figs.py`) still has its four placeholders; fill it after
the 3064 / 3088 rerun.

**13.2 Template fits of all beam runs: done** (2026-10-29/30, table in `new_runs.md`): 3064 gives 137 M KCU0 and 229 M KCU1 fitted
pulses, 18 x run 3062; the low-rate runs 3065 / 3069 use the 3064 templates (`templateFile` flag).

**13.3 Gain stability over the 12 h of run 3064 (`calibration/stability/gain_stability.{C,_figs.py}`, `gain_stability.txt`,
`gain_stability_{rate,kcu0,ns_ped,clip,kcu1}.png`).** Per 600 s bin the interpolated median of the fitted amplitude. Run structure
from the single-paddle shower rate (median 952 Hz): beam-off gaps at 1.5-2.1 h and 5.8-6.3 h, beam until 9.08 h, then the TAC
period at 3-6 Hz of showers to 12.3 h. Before 9.08 h: SFIL1-4 South beam cells 0.8 / 0.9 / 0.8 / 1.3 % rms of the 10-min medians
with slopes of 0.00 / -0.07 / -0.06 / -0.13 %/h; the SFIL1-4 South total 0.68 % rms, -0.02 %/h; SFIL3 per paddle 0.6-1.1 % rms.
The four layers move together (+1 % at 0-3 h, -1.5 % at 4-5.5 h, +1 % after 6.5 h) and ln(A_N/A_S) of the beam cells moves with
them (+0.015 at 4-5.5 h; rms 0.006-0.007, no slope), so the common drift is the beam (energy or x along the fibres), not the
gains: **the gains are stable to better than 1 % over the run**; no time-dependent correction is needed for 3064 and the 3062
calibration carries over. Pedestals of the beam cells stay within +-2.5 ADC (63-88 ADC). In the TAC period the medians are 2-3 %
lower (different beam). KCU1 SFIL5 South unclipped pulses, per array: 0.6-0.8 % rms, slopes <= 0.06 %/h. Fraction of the fitted
South beam-cell pulses at the 1023 ceiling over the run: SFIL2 0.18, SFIL3 0.53, SFIL4 0.34 (3062: 28-45 / 53-66 / 22-53 % by paddle).

**13.4 Cosmic MIP scale from run 3088 (`calibration/cosmics/cosmic_mip.C`, `cosmic_mip_fit.py`, condor `run_cosmic.{sh,job}`,
output `kcu0/cosmics/cosmic_3088_<start>.root` and `cosmic_3091_0.root`).** What the waveforms are (30 k-trigger study): the muon
pulses sit at the BEAM timing (rise at sample 6, peak 7-9, the usual shape, 40-70 ADC for a cell crossing); the "late pulses at
samples 12-15" of the first look were slow baseline ramps whose rms grows from 7 ADC at sample 2 to 25 ADC at sample 19. The ramp
is common to the connected channels of an ASIC half (corr 0.66-0.98, strongest within a summing board; the unused channels of the
half do not see it, corr -0.1): subtracting per sample the median over the connected channels of the half (SFIL1 North column
excluded) brings the noise from 5-9 to 0.8-2.5 ADC on the clean channels. Residual ramps leak into the pulse window on SFIL1,
SFIL2 c3-4, SFIL4 S c3-4 and SFIL5 N (2-12 % of triggers above 15 ADC against 0.02-0.6 % on the clean channels), rejected by a
pulse-shape cut on the stored waveform (sample peak+4 below 0.8 of the peak). The cosmic trigger is random in time (median
spacing 44 ms, 1 % periodic) and half of its triggers carry a ~90 ADC dip in exactly one paddle: the paddles are part of it. The
muon tag is the coincidence of the two ends of one cell (different ASICs, independent noise), both peaking at samples 7-9 within
one sample of each other, with the pulse shape at one end at least. **Full run (`cosmic_mip_3088.txt`, `cosmic_mip_{gm,ends,
vertical}_3088.png`):** per cell SFIL1 17-45 k and SFIL2 c3-4 7-11 k tags, all accidental (the SFIL1 North column and the SFIL2
c3-4 ramps); SFIL2 c1-2 550-920, SFIL3 350-530 and SFIL4 200-380 tags with accidental expectations of <= 36 (SFIL3) to a third
(SFIL4): real particles. Their sqrt(A_S A_N) spectra peak at 15-20 ADC with long tails (SFIL3: Landau MPV 16.5 +- 1.7 for c2,
medians 21-34; SFIL4 medians 10-13, its North end weak). What they are: 100 % carry a paddle dip (all triggers 52 %), 84 % exactly
one; a SFIL3 c1 tag lights c2 in 80-93 %, a c3 tag lights c2 / c4 in only 9-26 %; **no track through three layers exists** (same
cell at both ends in SFIL2, 3 and 4: 0-4 per cell in 984 k triggers), and the in-board "vertical" tags on SFIL1 S, SFIL2 S and
SFIL5 N (thousands) are within-board bursts, not muons (corr(A_c1, A_c2) = 0.97-0.99, ratio 0.91; SFIL3, with the same flux,
has 4-26) - the SFIL5 North "MPV 11-14 ADC" of the first look is therefore withdrawn. So the trigger of 3088 does not define the
track geometry, and without the path length per cell the 15-20 ADC peaks are not yet an ADC-per-MeV scale. Next step: a cosmic
mode of the Geant4 model (muons from above, cos^2 zenith, flat over the stack) with the same two-end tag, which gives the
path-length-weighted visible energy per tag and turns the SFIL3/4 spectra into gains (and says what fraction of the two-end
events are muons); and run 3089 (47 h) when it arrives. SFIL5 South (KCU1) cannot be MIP-calibrated from cosmics (nothing
above 30 ADC). The chain itself (common-mode subtraction, timing tag, spectra, fits) is in place.

**13.4b Off-beam cells from the fit trees: not possible** (`calibration/cells/cell_ratios_fits.py`, 3062 and 3064 parts 00 / 07):
the off-beam cells carry <= 0.07 of the beam cell (~40 ADC), below the fit's 60 ADC threshold, so they have a fitted pulse in
0-2 % of the showers (SFIL4 North c1-2 12-17 %, the constant ~130 ADC channels); the raw-peak method stays the tool:
`cell_equalisation.py` now takes a file and tag and runs on 3064 part 00 on condor (cluster 429056, `cell_equalisation_3064_part00.*`).

**13.5 Baby BCAL.** `caen_calib/caen_apply.py` now takes any CAEN file (`FILE TAG`); the per-cell DRS4 pedestals of 3062 are
applied to 3088, 3091 and the five 3064 parts (cluster 429051, `caen/calib/caen_calib_<tag>.npz`). Raw 3088 CAEN waveforms show
the cell pattern (26 ADC rms, all 12 channels above 100 together) and nothing else yet; 912 k of its 1.01 M triggers have no KCU
partner (to be rechecked with the repaired alignment), so the CAEN may have had its own cosmic trigger: look at the corrected
data for per-channel spectra with column patterns and MIP peaks for the rows 0-3 that the model-anchored gains need (the
factor-3 fall with depth, section 7).

**13.6 Simulation.** Geometry v4 (survey) is the default. **v4 energy scan done (cluster 429047, 8 points: 3-6 GeV at x = +30
and +80 mm, 5000 e+ each, seeds 81-88, `sim/v4/`; `scan.sh` takes an output sub-directory as 7th argument and an `EXTRA`
environment string of macro commands; `sim_summary.py tag=_v4 ...` -> `sim_summary_v4.*`, `stack_shares_pe.py sim/v4 _v4`,
`stack_shares_calib.py sim/v4 _v4` -> `sfil_gain_calibration_v4.txt`).** Against v3: SF unchanged (0.1022), the Baby BCAL visible
share rises 6-7 % (BBvis/E 0.030 vs 0.028 at 4 GeV) and the SFIL share falls 1 % (the 51 mm gap and the survey z moved more of the
shower into the back); stack containment 0.77 (0.75); sigma/mu of the visible energy 4.2 / 3.9 / 3.6 / 3.4 % at 3-6 GeV (v3 4.5 /
4.1 / 3.8 / 3.8). SFIL1-4 layer shares at 4 GeV 0.047 / 0.203 / 0.356 / 0.393 (v3 0.051 / 0.210 / 0.355 / 0.384): the survey
pushes a little light from SFIL1-2 to SFIL4. The data comparison is unchanged in its conclusion: **pe-equalised SFIL2-4 shares of
3062 give E = 3.0 GeV (3.0-3.3) for all paddles, P5 3.3 (3.0-4.0)**, as with v3. Implied South gains at the 4 GeV / P6 anchor
(ADC per MeV of South light): **SFIL1-4 = 8.6 / 24.2 / 18.4 / 14.7** (v3: 7.8 / 22.8 / 18.0 / 14.8), i.e. SFIL1 +10 %, SFIL2 +6 %,
SFIL3 +2 %, SFIL4 -1 %; at 3 GeV 9.2 / 28.1 / 22.9 / 19.6. The v4 scan also has P5 (x = 80) points at 3 and 4 GeV (P5 4 GeV:
8.8 / 25.4 / 19.4 / 16.7). The per-paddle figures of section 7 that quote v3 gains stay as they are (the differences are inside the
17.5 % gain band); the cosmic comparison below uses the v4 geometry. The paddle-position conflict with the AstroPix chip map
(sim/README.md) still needs the survey author.
Two more model updates from the available data (2026-10-30 afternoon):
- **v4b light model** (cluster 429059, `sim/v4b/`): the FTBF two-component attenuation replaced by the single exponential measured
  on SFIL2/3 from the AstroPix x (lambda_eff 1500 mm, between 1329 and 1646) and v_eff 160 mm/ns (N-S timing 157-162); the same 6
  energy points. **Done (`sim/v4b/`, `*_v4b` outputs):** the layer shares are identical to v4 (SFIL1-4 at 4 GeV 0.047 / 0.202 / 0.357 /
  0.393), so the shape energy stays 3.0 GeV (3.0-3.3; P5 3.2); only the absolute South light per MeV changes, and the implied 4 GeV / P6
  gains rise 4 % to 9.0 / 25.5 / 19.2 / 15.4 ADC per MeV of South light (v4 8.6 / 24.2 / 18.4 / 14.7); N-S timing rms 0.18 ns (v4 0.17).
  The light model is a 4 % term on the absolute gains and nothing else.
- **Cosmic mode** (`/bic/gun/cosmic 1`, `macros/cosmic.mac`, `condor/cosmic.{sh,job}`; sim/README.md): muons from a plane 250 mm
  above the setup, cos^2 zenith, log-flat 0.5-20 GeV; 2 M muons produced (`sim/cosmic/`, 20 files). 2-3 % of them cross a given
  layer; a SFIL3 cell crossing leaves 1.56 MeV visible and ~1.3 MeV of light per end, so the 18 ADC/MeV beam anchor predicts ~23 ADC
  per MIP, where the 3088 two-end spectra peak. `calibration/cosmics/cosmic_mip_sim.py` fits per cell the light yield Y (pe per MeV
  of end light) by the likelihood of the data sqrt(A_S A_N) spectrum against the simulated crossings read out as Poisson(L Y) x G +
  noise (G = ADC/pe from the 3062 photostatistics, North = South x the measured N/S ratio), with the same both-ends > 8 ADC tag, and
  gives Y x G in ADC per MeV of end light next to the beam-anchored gain; repeated for 0.7 G and 1.4 G. **Result with the 2 M
  production (`cosmic_mip_sim_3088.txt/.png`, 31 k simulated crossings per cell, isotropic flux): SFIL3 cells 1-4 = 19.7 / 19.7 /
  27.8 / 27.8 ADC per MeV of end light** (likelihood intervals +-5-10 %, 15-26 for 0.7-1.4 G) against the 4 GeV beam anchor's 18.0
  (15.4 at 5 GeV): the cosmic scale of SFIL3 is 10-55 % above the 4 GeV-anchored gain, i.e. the P6 showers deposit less than the
  model's 4 GeV, 2.6-3.6 GeV if the light yield is the same for muons and showers. SFIL1 / SFIL2 / SFIL4 fall to the fit floor: their
  data tags are accidentals (SFIL1, SFIL2 c3-4) or the weak North end (SFIL4, N/S 0.29: the simulated tag efficiency is 1 %), as
  13.4 said. **In the simulation only 1 % of the muons hit a paddle and almost none a paddle and a SFIL cell** (19-50 of 31 k
  crossings), while 100 % of the data tags carry a paddle dip; the paddle-tagged template of SFIL3 c1 gives 16.9 instead of 19.7
  ADC/MeV, so the track geometry matters at the 15 % level. **Trigger geometry resolved (2026-10-30 evening):** the paddle-footprint
  production (`sim/cosmic/pad/`, 2 M muons, 4 % hit a paddle) still has only 2-12 paddle + SFIL3 crossings per cell, because the
  paddles are vertical slabs upstream of the stack at beam height: a muon through a paddle AND a SFIL cell is near-horizontal (zenith
  ~72 deg) and starts 0.5-1 m upstream of any start plane. In the data the paddle pulse (positive, ~88 ADC, exactly one paddle) is on
  52 % of all 3088 triggers (51 % of the triggers with no SFIL activity at all) but on 97-100 % of the SFIL3 two-end tags, spread over all
  four paddles, and the tag rate per trigger is 45 x higher with a paddle pulse: the tagged muons did pass through a paddle, at 1.4 tags per
  minute for the whole layer, which is the rate expected for near-horizontal tracks through the 33 cm^2 paddle faces into SFIL3 (a paddle
  laid flat on the stack would give ~30 per minute). The data also show SFIL2 South activity in ~60 % of the SFIL3 tags (reference 1 %),
  as such tracks cross SFIL1-3 in sequence (SFIL4 cannot be tested: its South channels sit above 8 ADC in 30 % of random triggers from
  their own bursts). The simulation got a **cosmic target mode** (`/bic/gun/cosmicTarget 1`: every track passes through one of the
  paddles as placed, flux-correct face sampling; sim/README.md): 94 % paddle hits, 0.16 % paddle + SFIL3 crossings, median South light per
  crossing 1.08 MeV against 1.13 for the plane production, so the template changes little and the earlier 16.9 from 49 crossings was
  statistics. A 10 M-muon production (`condor/cosmic_tgt.job`, cluster 429062, `sim/cosmic/tgt/`, ~1000 crossings per SFIL3 cell) is
  queued; `cosmic_mip_sim.py ... 3088tgt pad` then gives the trigger-matched template.
  **Cells 1 and 2 share one SiPM pair** on the v2 boards (both channels read SiPM1 + SiPM2 mixed, r = 0.97-0.99; hookup note): a muon in
  physical cell 1 or 2 lights both channels (SFIL3 c1 and c2 tag together in 75 % of the cases, c2 and c3 in 6 %), so the per-channel fits of
  c1 and c2 see a share of the light and their lower gain (19.7 against 27.8 for c3 / c4) is expected. `cosmic_mip_sim.py` now adds the
  element "c1+2" (channel sums per end against the light of physical cells 1 + 2). **Result (isotropic template, likelihood scan now with
  common random numbers, `cosmic_mip_sim_3088.txt`): SFIL3 c3 = c4 = 27.8 ADC per MeV of end light (24-28 for 0.7-1.4 G), c1 / c2
  17.8 / 19.7, c1+2 summed 18.7 (21.5 for the G variants).** The summed pair stays at ~0.7 of the c3 / c4 gain, so the channel sharing
  does not explain the whole difference: either the passive network loses charge or the SiPM1+2 pair has a lower gain; the beam-row cell
  c3 and the independent c4 agree with each other and are the values to use. **Against the v4 4 GeV / P6 anchor of 18.4 the beam row's
  cosmic gain is 1.5 x higher: the P6 showers deposit about 2 / 3 of the 4 GeV light, i.e. ~2.7 GeV (2.6-3.1 over the G variants), in
  agreement with the 3.0 GeV (3.0-3.3) from the layer-share shape of section 7 / 13.6.** Two independent handles (absolute MIP scale and
  longitudinal shape) now say the same thing. **Trigger-matched template done (10 M target-mode muons, cluster 429062, 94 % through a
  paddle, 5.5-9.7 k paddle + SFIL3 crossings per cell; `cosmic_mip_sim_3088tgt.txt/.png`): SFIL3 c3 26.5, c4 27.8 (25-27 for 0.7-1.4 G),
  c1 / c2 per channel 18.7 / 19.7, and the c1+2 channel sum 29.2 (30.4 for the G variants).** With the right track geometry (a near-
  horizontal track crosses only one of the two physical cells of the shared pair) the summed pair agrees with cells 3 and 4: the SFIL3
  board gain is 27-29 ADC per MeV of end light in all four cells, and the per-channel c1 / c2 values are what a channel that carries
  ~2 / 3 of the pair's charge shows. The template choice moves c3 / c4 by < 5 % (isotropic 27.8 / 27.8). **Final cosmic scale for SFIL3:
  27 +- 2 ADC/MeV (fit 4 %, ADC/pe 7 %, template 3 %), 1.45-1.55 x the v4 4 GeV / P6 anchor of 18.4 (1.4-1.5 x the v4b 19.2), so the
  P6 showers deposit 2.6-2.9 GeV-equivalent light.**

**13.7 Performance.** The 3062 chain (`sim/analysis/performance.py`) needs for 3064: the KCU0 fits (done), the per-event layer
sums (`extract_fit_layers.py`, now `[fit file] [output]`, submitted for the 14 parts, cluster 429060 -> `data_fit_layers_3064_partPP.npz`;
the parts are slices of ~3.85 M triggers, 0.65 h each; the TAC period after 9.08 h is ~50k triggers at the end of part 13, not a part of its own). **Done for SFIL1-4 (`performance_3064.py` -> `performance_3064.txt/.png`, 53.6 M triggers, 34.1 M showers):** with the v4 4 GeV gains and the 3062 definitions, per paddle P4 / P3 / P6 / P5: <E_vis> 114.2 / 126.2 / 130.1 / 139.6 MeV, core sigma/mu 29.6 / 26.7 / 20.6 / 19.1 % (3062: 29.9 / 26.5 / 20.4 / 18.8 %; model SFIL1-4 alone 17.6-18.0 %), pe-equalised shape energy 3.00 / 3.00 / 3.27 / 3.88 GeV (all 3.02, 3.0-3.7), layer shares identical to 3062 within 0.001. Over the 14 parts <E_vis> per paddle has an rms of 0.3-0.4 % and sigma/mu of 0.001-0.002: 3064 is 3062 with 17 x the statistics and no drift, so the 3062 performance numbers (section 11) stand for both runs. Still missing for the full stack: the KCU1 rebuilt SFIL5 sums (`kcu1_tot/kcu1_sat_recovery.py` on
the 3064 KCU1 fits, still 3062-only), the CAEN cell correction (queued) and the CAEN alignment (13.1, queued). The first 3064
result to aim for is the energy resolution per paddle with the TAC period excluded and the gains of 3062 unchanged; the 0.7 %
stability of 13.3 is the systematic floor for that. (`calibration/stability/resolution_3064.py` gives the width of the stability
sums per paddle, 34-39 %, but those include showers with missing layers and are not the 3062 definition; use the layer-sum files.)

**13.8 Not possible with this data:** KCU1 SFIL5 N-S sub-ns timing (item 12); the SFIL5 South TOT slope (LED runs); the absolute
beam energy per paddle (PS momentum map).

## 14. Documentation index, plan forward and requests (2026-11-03)

**Where each result is written up.** Every number in the artifacts below traces to a file named in this document.

| document | date | covers | state |
| --- | --- | --- | --- |
| `docs/status.md` (this file) | 2026-11-03 | every study, result and open item | the overview; read first |
| Deck "BIC test beam update, November 3" (https://claude.ai/artifact/NSsUaqzfPWswYWRokddbfQ; pptx + preview PDF in `docs/slides_nov03/`) | 2026-11-03 | reconstruction and tools, calibration chain (merged from the October 26 deck: sync figure, chain table, three-handle relative gains, CAEN cells, five-element shares, energy per paddle), simulation (geometry, readout model and versions, scan, model against data), performance for 3062-3091, plan and requests; 44 slides in the artifact (the pptx in the folder is the earlier 36-slide build) | current |
| Claude Doc "BIC Hall D test beam 2026: run 3062 analysis summary" (https://claude.ai/artifact/R45o42e7GoZhS2vi9WgVFK) | 2026-10-23, update section 2026-11-03 | the long-form write-up with figures: setup, sync, decoders, waveforms, trigger, v3 model, calibration, performance | sections 7-9 describe the v3 model and the 4 GeV-anchored gains; the update section at the top carries the current numbers |
| Deck "BIC Run 3062 Status" (https://claude.ai/artifact/E1KkKtTPNqQwKk5F2pwoEu) | 2026-10-26, content updated 2026-11-03 | the calibration chain step by step at the 3062 level | current (v4 model, cosmic scale, ~3 GeV energy, PRs); its content is merged into the November 3 deck |
| Deck "SFIL2 TOA Calibration" (https://claude.ai/artifact/VRuZZSDjD3ZQueEXeQ4Ttv) | 2026-10-28, event-time slide added 2026-11-03 | the ADC-free KCU0 TOA study | current |
| Claude Doc "BIC Raw Data Format: KCU, CAEN and AstroPix Streams" (https://claude.ai/artifact/HyG8TDatu9MYd2v6F78c3n) | 2026-10-21, update 2026-11-03 | raw packet, line and word layouts; TOA / TOT fields; clocks | current; formats verified bit-exact |
| `docs/method.md`, `data_format.md`, `root_files.md`, `completeness_plan.md`, `new_runs.md`, `subsystems.md`, `channel_map.md`, `pitfalls.md`, `beam.md`, `stack_shares_calib6.md` | Oct 21-30 | methods per stream, formats, trees, DAQ loss, new runs, subsystems, hookup, pitfalls, beam | current |
| `sim/README.md` | 2026-10-30 | the Geant4 model v1-v4b, scans, cosmic mode | current |
| `sync/note/sync_note.tex` | draft | decoding and synchronization of runs 3061-3091 | four placeholders left (%%RUNS%%, %%DECODE%%, %%SYNC%%, %%SUMMARY%%); fill with `run_sync_figs.py` |
| `dual_sided/note/dual_sided_note.pdf` | 2026-10-28 | the 2017 dual-sided readout study repeated on 3062 | current |

**Studies in the repository and where they are documented** (checked 2026-11-03 against every folder with results): decoding/ and decoder_tests/ (section 2, raw-format doc), sync/ (3, 13.1), waveforms/ and waveforms/kcu1_tot/ (5, 7), timing/toa_sfil2 + toa_clean (section 7 and the TOA deck; timing/toa_position is the withdrawn ADC-tuned version, kept for the record), astropix/ (10, 7: imaging, ToT, geometry checks, SFIL5 array order), calibration/stability (13.3), cosmics (13.4, 13.6), cells (7, 13.4b), attenuation (7), caen_calib (7, 13.5), paddles (4), rate (7), dual_sided/ (7), sim/ (6, 13.6), event_displays/ (method 4.x, figures only). Nothing with a result is undocumented.

**Plan forward, in order of leverage.**
1. Full-stack performance of 3064: join the rebuilt SFIL5 sums (`kcu1_tot/kcu1_sat_recovery.py` generalised to the 3064 KCU1 fits) and the cell-corrected, aligned CAEN to the 34 M KCU0 showers; resolution and linearity per paddle with the 3062 definitions.
2. Baby BCAL MIP peaks per array from the 3088 CAEN data (own trigger): decides the factor-3 gain fall with depth; then the Baby BCAL share and containment.
3. Energy scale: re-anchor the gains at 3 GeV (cosmic scale and shape agree), recheck the SFIL5 and Baby BCAL double ratios across the paddles, and redo section 10 at that anchor.
4. Decoding and synchronization note: fill the four placeholders from `output/sync_runs/run*/sync.txt` and `run_sync_figs.py`.
5. AstroPix: hot-pixel scan; the SFIL5 South map when the pin map arrives.
6. Next run (hardware / DAQ): TOA threshold scan per channel and a delay scan, socket-buffer / recvmmsg against the UDP loss, the common busy, the DAQ YAML in the begin-run record.

**Requests to the group.** (1) The pair-spectrometer momentum map for August, or any run with a known energy: the data say ~3 GeV under P6 and the whole energy scale hangs on it. (2) The survey versus the AstroPix chip map: P3 and P6 are 9.5 mm apart in the survey and 60 mm apart in the data. (3) Trigger definitions of the cosmic runs (paddle OR plus what?) and of the low-current runs 3066-3069. (4) The missing raw runs: cosmics 3089 (47 h), pedestal 3173-3178, LED 3179-3190, injection run 2804. (5) The SFIL5 South SiPM board pin map and the Baby BCAL FTBF MIP constants per array. (6) The AstroPix module-8 bench source scans. (7) The KCU0 TOA configuration of run 3062 and the DAQ YAML (analog gains).

