# New raw data on `/gpfs01/eic/data1/BIC_Beamtest` (survey 2026-10-28)

Added on 2026-10-27 (owner purschke): **beam runs 3061, 3063-3069 and cosmic runs 3088, 3090, 3091**, about 0.73 TB of
RCDAQ files, every run split into the four streams `{beam,cosmics}_{ROC_0,ROC_1,CAEN,AstroPix}-0000RRRR-SSSS.evt` (ROC_0 =
KCU0 packet 12001, ROC_1 = KCU1 packet 12002, CAEN packet 2001, AstroPix packet 13001; clock packets 1010-1013), exactly the
format of run 3062. Run 3062 itself is not in this directory (it stays in `eic/data/`).

How this was checked: run-log sheet (`context/runlog_sheet/Runlog.csv`, `ASIC_Mapping.md`), `dlist -i` on the first event of
each stream, file sizes scaled to the run-3062 record counts, and a 300-packet conversion of five streams with the existing
`kcu_to_root` (0 ADC mismatches against the stock pool on every one) read with `decoding/newrun_qa.py`. Nothing beyond 300
packets per stream was decoded; the CAEN streams of the new runs were not opened.

## Catalog

Records: KCU = triggers (scaled from the file size with run 3062's 6.4 kB/trigger for KCU0), CAEN = triggers, AstroPix = hit records
(AstroPix column from the converted files, 2026-10-28; the first survey quoted packet counts there).

| run | type | date (CT), length | KCU0 | KCU1 | CAEN | AstroPix | size | conditions (run log) | verdict |
| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |
| 3061 | beam | Aug 13 20:31, 1 min | ~1 k | off | ~1 k | off | 0.01 GB | KCU1 out, AstroPix off | test run, skip |
| 3063 | beam | Aug 14 00:09, 2 min | 235 k | 100 k | 104 k | 165 k | 2.5 GB | same as 3062 (TPOL B, JD80 radiator, 146 nA) | **extend 3062 (+6 %)** |
| 3064 | beam | Aug 14 00:12-12:30, 12 h | **~69 M** | 28.5 M | 28.7 M | 1.51 M | 707 GB | TPOL C (Be ~750 um), radiator undefined, 128 nA average; "low-current TAC run began during this run" | **the main new sample: 18 x run 3062** |
| 3065 | beam | Aug 14 12:30, 1 min | 1.7 k | 1.5 k | 1.7 k | off | 0.03 GB | low current | skip |
| 3066 | beam | 12:31, 20 min | 36 k | 34 k | 34 k | 20 k | 0.6 GB | 0.78 nA, Al radiator 1.5 um, "possible TAC run?" | low rate, see below |
| 3067 | beam | 12:52, 20 min | 34 k | 33 k | 33 k | 18 k | 0.5 GB | 0.73 nA, same | low rate |
| 3068 | beam | 13:12, 21 min | 30 k | 28 k | 29 k | 13 k | 0.5 GB | 0.59 nA, same | low rate |
| 3069 | beam | 13:33, 13 min | 17 k | 16 k | 16 k | 4.8 k | 0.3 GB | 0.34 nA, same | low rate |
| 3088 | cosmics | Aug 14 23:11 - Aug 15 14:07, 15 h | 1.1 M | 0.84 M | 0.97 M | off (0-byte file) | 15.6 GB | new hookup (below) | **MIP candidate sample** |
| 3090 | cosmics | Aug 17 13:03, 1 min | 0.6 k | 0.5 k | 0.6 k | 0 | 0.01 GB | | skip |
| 3091 | cosmics | Aug 17 13:04, 20 min | 13 k | 10 k | 12 k | 5.4 k | 0.2 GB | AstroPix on | small MIP sample |

**Still not on any disk we can read:** run 3089 (47 h cosmics, 1.67 M CAEN events, the largest cosmic run), cosmic runs
3092-3172 (AstroPix on), pedestal runs 3173-3178, LED runs 3179-3190, the HGCROC injection run 2804, and beam run 3060
(SFIL1-4 only, KCU1 out, 1.24 M KCU0 packets).

## What each new run is good for

**Run 3064 (beam, 12 h).** The 300-packet sample looks like run 3062: SFIL2 cell 3 has a pulse on 63 % of triggers (3062:
65-68 %) with median 500-600 ADC (3062: 515-651), the paddle fractions are P3 0.25, P4 0.23, P5 0.14, P6 0.10 (3062: 0.24,
0.23, 0.13, 0.10), the KCU0 lines are 99 % complete, 0.1-0.3 of the beam-row pulses clip as before. Same PS field (1.662 T),
so the paddle energies should be the same; the converter target changed (TPOL B -> C) and the radiator is "undefined".
Useful for every calibration step that is statistics-limited in 3062:
- per-cell and off-beam-cell gains (step 3 needs the cell equalisation that 3062 cannot separate from beam height);
- SFIL5 TOT slope and the KCU1 saturation rebuild (28.5 M KCU1 triggers);
- gain stability over 12 h (SiPM current 343 uA in 3062 vs 164 uA average here; no temperatures logged);
- Baby BCAL per-array gains and depth trend with 19 x the CAEN statistics;
- dual-sided timing and good-TOA studies with ~18 x the showers.
Weak for AstroPix matching: 1.51 M AstroPix hits in 12 h, about half of what 3062 has in 25 min (2.86 M), so AstroPix was
off or not triggered for most of the run (when it ran still has to be found from the record times). Attenuation from AstroPix x therefore stays with 3062; attenuation from the paddles works.
The run-log note says a low-current TAC run began during 3064: the end of the run has to be found (trigger rate vs time)
and treated separately.

**Runs 3066-3069 (beam, ~100 k triggers, 0.3-0.8 nA).** Real showers in the beam cell (SFIL2 cell 3 on 35-38 % of
triggers, median 300-430 ADC; SFIL3 cell 3 ~830 ADC), but **a paddle fires on only 4-10 % of triggers** (3062: 99.85 %), so
the trigger is not the paddle OR. AstroPix recorded 5-20 k hits per run (0.3-0.6 per trigger; 3062: 0.76). Possible use:
pile-up-free pulse templates and AstroPix matching at low rate. Not usable for energy or position per paddle until the
trigger of these runs is known.

**Runs 3088 and 3091 (cosmics).** Mapping after Aug 14 21:55 CT (`ASIC_Mapping.md`, "After 14/08"): unchanged for SFIL1-4
and paddles, but **the SFIL5 North summed boards moved onto KCU0 ASIC2, channels 153 / 156 / 163 / 167 (SFIL5 cells 1-4
North)**. In our channel map these are "unused", so they land in the `unused` tree; they do carry pulses (5-28 % of triggers
above 15 ADC). The sample shows:
- pulses arrive late, peak at samples 12-15 (beam: 7-8), so the pulse window must move;
- the per-trigger common baseline shift is large (pedestal spread ~10 ADC raw, 4.8 ADC after removing the per-ASIC median);
- after that, candidate vertical-muon signals are about 20 ADC per cell in SFIL1 and SFIL4 South (4-5 x the noise),
  while some cells show almost no hits; with 2.6 k triggers this does not yet show a MIP peak;
- paddles fire on 2.5-8 % of triggers, so the cosmic trigger is not the paddles either (definition unknown);
- SFIL1 North cells 1, 2, 4 are noisy (pulses on 78-99 % of triggers in 3066 and 3091), as in 3062;
- SFIL5 South single SiPMs (KCU1) show nothing above 30 ADC (0.1-0.2 % of triggers): no MIP calibration for them from
  cosmics;
- Baby BCAL North (KCU0 ASIC3) still silent / railed, as in 3062.
If the MIP peak is confirmed on the full run, 3088 gives the first data-only (model-free) absolute scale for SFIL1-4 and the
first SFIL5 North data. The CAEN (Baby BCAL South) cosmic data were not opened yet.

## What is needed to use these runs

1. **Run-generic condor jobs.** The converters work unchanged, but the job wrappers in `sync/condor/` hard-code run 3062,
   the `beam_` prefix and `eic/data/`. Parameterise run, prefix and input directory.
2. **Output space and time for 3064.** At run-3062 rates (5 kB per KCU0 trigger, 0.13 s per packet) the full conversion is
   ~0.35 TB (KCU0) + ~0.12 TB (KCU1) + ~0.13 TB (CAEN) on `/sphenix/tg/tg01/jets/tmengel/eic/` (14 PB volume, 60 % used),
   and ~270 CPU hours for KCU0 (140 segment jobs of ~2 h). Alternative: a skim that keeps only the SFIL, paddle and KCU1
   channels with a pulse, or run the template fit directly from the raw files.
3. **Run-dependent channel map** in `kcu_to_root` for runs >= 3085 (SFIL5 North on KCU0 ch 153 / 156 / 163 / 167).
4. **Per-run synchronization**: KCU0-KCU1 tick offset, CAEN clock ratio and offset, AstroPix phase, and the common trigger
   table / event tree are run-specific (`sync/reproduce.sh` assumes 3062).
5. **Trigger definitions** for the cosmic runs and for 3066-3069 (neither is the paddle OR): from the RCDAQ setup packet 900
   (if written) or the shift crew.
6. **Cosmic pulse finding**: per-ASIC common-mode subtraction, pulse window at samples 10-17, then per-channel MIP spectra on
   tracks (all four cells of a layer for vertical muons, or matched hits across layers for inclined ones).
7. **Find the TAC part of 3064** from the trigger rate against time before using the run as one sample.

## Conversion (submitted 2026-10-28)

All 269 non-empty segment files of the runs above are being converted to the run-3062 tree format by condor
(`sync/condor/run_to_root.{sh,job}`, list `run_to_root_list_all.txt`; test job cluster 419385, main cluster 419386).
Outputs: `/sphenix/tg/tg01/jets/tmengel/eic/<stream>/run<RUN>/<stream>_<RUN>_seg<k>.root` for stream = kcu0, kcu1,
caen, astropix; merged per run by `run_merge.{sh,job}` into `<stream>_<RUN>.root`, symlinked as
`bic_testbeam/data/<stream>_<RUN>.root`. `decoding/condor/run_watch.sh` (nohup on the login node, log
`output/root_logs/run_watch_419386.log`) waits for the cluster to drain, checks every segment log and submits the merges only
for stream/run pairs whose segments all exited 0 (plain condor_submit, no DAGMan); failed segments are listed in the log.

Changes for these runs: `kcu_to_root` takes the run number as an 8th argument and, for runs >= 3085, writes the SFIL5 North
summed boards (KCU0 ch 153 / 156 / 163 / 167) to a new tree `sfil5n` (`Nc1..Nc4_*`, map labels "SFIL5 cell c North end"),
leaving 176 channels in `unused`; runs < 3085 are bit-identical to before (checked on 300 packets of 3064). `kcu_merge`
carries `sfil5n`. `astro_to_root` titles the file with its input name instead of "run 3062".
Expected size: ~0.6 TB, nearly all from run 3064. Not done by these jobs: synchronization / trigger table per run.

**Result (2026-10-28).** All 269 segment jobs exited 0. Merged files (triggers in `header`; AstroPix hit records):

| run | KCU0 | KCU1 | CAEN | AstroPix | note |
| --- | --- | --- | --- | --- | --- |
| 3061 | 1 097 | 0 (KCU1 out) | 1 075 | 0 | |
| 3063 | 236 051 | 96 505 | 103 592 | 165 280 | KCU0 complete 89.7 % |
| 3064 | parts | parts | parts | 1 506 973 | merged in parts, below |
| 3065 | 1 680 | 1 541 | 1 711 | 0 | |
| 3066 | 34 451 | 34 467 | 34 591 | 20 221 | |
| 3067 | 33 104 | 33 173 | 33 325 | 18 204 | |
| 3068 | 29 391 | 29 303 | 29 427 | 13 066 | |
| 3069 | 16 567 | 16 566 | 16 593 | 4 816 | |
| 3088 | 1 035 395 | 1 027 992 | 1 007 790 | off | KCU0 has the `sfil5n` tree |
| 3090 | 593 | 593 | 599 | 0 | |
| 3091 | 11 947 | 11 947 | 11 945 | 5 394 | |

Run 3064 does not fit one file (~350 GB for KCU0, and ROOT's 100 GB tree limit aborted the first `hadd`), so it is merged in
parts of 10 segments: `run_merge.sh <stream> 3064 <part>` (list `run_merge_list_3064.txt`, cluster 419519) writes
`<stream>_3064_partPP.root` (KCU0 parts 00-13, KCU1 and CAEN 00-04, ~25-30 GB each), symlinked under the same name in
`bic_testbeam/data/`. Every part starts on a de-duplicated head file, so chaining the parts in order gives each trigger
once. Each KCU file has exactly one backward `ts0` step, flagged in `header.backward`: the first few triggers carry a stale
clock from before the reset at run start (3063: 4 triggers at 1486 s, then the clock restarts at 0); drop them in the
per-run synchronization.

## Template pulse fits (submitted 2026-10-29)

The next stage after the conversion: `kcu0_sfil_fit.C` (SFIL1-4) and `kcu1_sfil5_fit.C` (SFIL5 South) over every beam run, as for
3062 (`waveforms/condor/run_fit.{sh,job}`, `Queue board, tag from <list>`; list `run_fit_list_beam.txt`, cluster 419662, 31 jobs:
KCU0 3063, 3065-3069, 3064 parts 00-13; KCU1 3063, 3065-3069, 3064 parts 00-04). Each job runs in its output directory
`/sphenix/tg/tg01/jets/tmengel/eic/<board>/fit_runs/<tag>/` (the macros write their figures to the working directory), logs
`output/root_logs/fit_<board>_<tag>.out`. Each 3064 part builds its own templates from its first triggers. The ACLiC libraries
were rebuilt (`.C++`) after the thread reorganisation, because the KCU1 one still pointed at the old `analysis/` path.
- 3061 (run locally, 1 097 triggers): no beam; 13 927 of 14 056 channel waveforms below `minAmp`, 0 fitted pulses.
- Not fitted: the cosmic runs 3088 / 3090 / 3091 (pulses at samples 12-15 and ~20 ADC, outside the beam fit's window; they need the
  cosmic pulse finding of to-do 6 in `status.md`); KCU1 3061 (no triggers).

Results of cluster 419662 (all exit 0; fitted pulses kept / rejected by the chi2 cut):

| run | KCU0 triggers | KCU0 kept | KCU0 chi2 | KCU1 triggers | KCU1 kept | KCU1 chi2 |
|---|---|---|---|---|---|---|
| 3063 | 236 k | 446 217 | 207 959 | 96.5 k | 811 658 | 201 385 |
| 3064 (14 / 5 parts) | 69.1 M | 137.2 M | 54.2 M | 28.4 M | 228.9 M | 55.2 M |
| 3066 | 34.5 k | 29 544 | 13 192 | 34.5 k | 129 644 | 27 308 |
| 3067 | 33.1 k | 25 631 | 11 161 | 33.2 k | 111 902 | 22 386 |
| 3068 | 29.4 k | 14 679 | 6 353 | 29.3 k | 58 848 | 11 185 |
| 3069 | 16.6 k | 4 024 (3064 templates) | 2 005 | 16.6 k | 13 721 | 2 202 |
| 3065 | 1.7 k | 1 471 (3064 templates) | 566 | 1.5 k | 5 710 (3064 templates) | 1 325 |

KCU0 3065 / 3069 and KCU1 3065 have < 100 clean pulses per channel in the template pass, so no channel got a template and every
fit failed. Both macros now take `templateFile` (flag; 3rd argument of `run_fit.sh`, list `run_fit_list_lowrate.txt`): channels
with too few pulses use the `hTemplate_k` of another run's fit output (KCU0: 3064 part 00, KCU1: 3064 part 04). Local check on
KCU0 3065: 1 471 kept, 566 chi2-rejected (28 %, as 3064's 28 %). Rerun as cluster 419797 (all exit 0, values in the table); the failed first pass is kept in
`/sphenix/tg/tg01/jets/tmengel/eic/fit_runs_firstpass/` and `root_logs/fit_*_firstpass.out`. The KCU0 figure captions now take the
run from the input file name (they said "run 3062" for every run; the figures of the first jobs still carry it). KCU1 3064
parts 00-03 end with a bash syntax error after ROOT exited (`run_fit.sh` was edited while they ran); their output files are
complete (not recovered, tree entries = kept triggers).

**Run 3064 performance, SFIL1-4 (2026-10-30 evening, `sim/analysis/performance_3064.py`):** 53.6 M triggers, 34.1 M showers; per paddle
P4 / P3 / P6 / P5 core sigma/mu 29.6 / 26.7 / 20.6 / 19.1 % and <E_vis> 114 / 126 / 130 / 140 MeV with the v4 4 GeV gains, the same as
3062 within 0.3 %; per-part (0.65 h) rms of <E_vis> 0.3-0.4 %. Shape energy 3.0 GeV (P5 3.9). Figure `performance_3064.png`.

## Cosmic run 3088: what the waveforms show (2026-10-30, `calibration/cosmics/`)

The sample study above (2.6 k triggers) said the cosmic pulses peak at samples 12-15. A 30 k-trigger study with the per-sample
common-mode subtraction corrects that:
- **the muon pulses sit at the beam timing** (rise at sample 6, peak 7-9, the usual shape; 40-70 ADC for a cell crossing);
  the late population is a slow baseline ramp, rms 7 ADC at sample 2 growing to 25 ADC at sample 19, shared by the connected
  channels of an ASIC half (0.66-0.98 correlated, strongest within one summing board), absent on the unused channels (-0.1);
- subtracting per sample the median over the connected channels of the half (SFIL1 North column excluded) leaves 0.8-2.5 ADC
  on the clean channels (SFIL2 c1-2, SFIL3, SFIL4 N, SFIL4 S c1-2); SFIL1, SFIL2 c3-4, SFIL4 S c3-4 and SFIL5 N keep residual
  ramps in the pulse window (2-12 % above 15 ADC against 0.02-0.6 %), which the pulse-shape cut of `cosmic_mip_fit.py` removes;
- the trigger is random in time (median spacing 44 ms) and half of the triggers have a ~90 ADC dip in exactly one paddle, so
  the paddles are in the cosmic trigger (the 2.5-8 % above used a dip > 40 over samples 8-16; the cosmic dips are at 6-11);
- muons: the two ends of a cell fire together (different ASICs, independent noise), 8-27 per cell per 29 k triggers in SFIL3-4,
  about 1e-3 per channel per trigger; no near-horizontal track through three layers in 30 k triggers; vertical patterns (three
  cells of one board) appear on SFIL1 S, SFIL2 N and SFIL5 N but are contaminated by the within-board noise;
- the in-board vertical tags (three cells of one board) on SFIL1 S, SFIL2 S and SFIL5 N are within-board bursts, not muons
  (corr(A_c1, A_c2) 0.97-0.99; SFIL3 has 4-26 such events against thousands): no SFIL5 North MIP value from them;
- no track through three layers in 984 k triggers; the SFIL3/4 two-end samples (200-500 per cell, 100 % with a paddle dip)
  peak at 15-20 ADC with long tails, path length undefined: the ADC-per-MeV scale needs a cosmic mode of the simulation
  (`status.md` 13.4).
Full-run jobs: `calibration/condor/run_cosmic.job` (cluster 429052, 3088 in four chunks + 3091), output
`/sphenix/tg/tg01/jets/tmengel/eic/kcu0/cosmics/`; then `cosmic_mip_fit.py` on the hadd. The tree stores the corrected
waveforms of triggers with any candidate (57 % of triggers, ~1 GB for the run) for the shape cuts and template fits.

**CAEN alignment of 3064** (section 13.1 of `status.md`): a direct scan shows the CAEN-KCU0 offset is constant at -0.140 s
over the whole run with the common clock ratio (800-1900 of 2000 CAEN triggers matched per beam slice); `run_sync.py` had
anchored on the empty pre-beam slice. Anchor choice fixed 2026-10-30, 3064 and 3088 resubmitted.

**Outcome (2026-10-30 evening, details in status.md 13.4 / 13.6):** the two-end muon tags are near-horizontal tracks through the upstream
trigger paddles (paddle pulse on 97-100 % of the tags, all four paddles, 1.4 tags per minute in SFIL3); the simulation's cosmic target
mode reproduces that population. SFIL3 cosmic scale 27 +- 2 ADC per MeV of end light in all four cells (cells 1 and 2 share one SiPM
pair and must be summed), 1.5 x the 4 GeV model anchor: the P6 showers of 3062 carry 2.6-2.9 GeV of light, in line with the 3.0 GeV
layer-share shape. SFIL1, SFIL2 and SFIL4 give no cosmic scale from 3088 (accidental tags, noisy channels, SFIL4 North at 0.29).
