# BIC test beam, run 3062: stream synchronization

Analysis of the Barrel Imaging Calorimeter (BIC) test beam at JLab Hall D (August 2026), run 3062.
The four DAQ streams are put on one clock, the KCU0 FPGA timestamp (6.25 ns ticks), and combined
into one per-trigger event tree.

| Stream | Link to KCU0 | Precision | Matched | Null control |
| --- | --- | --- | --- | --- |
| KCU1 (H2GCROC, SFIL5 South) | fixed offset -1,143,665,946 ticks | 0 or +25 ns | 60% of KCU1 triggers (841k of 1.39M) | 0.0% (+10 ms shift) |
| CAEN DT5742 | clock ratio 1.365340684 + running-median correction | residual MAD 6.4 ns | 57.1% of 1.51M | 0.1% (+10 ms shift) |
| AstroPix 9-chip module | on-chip ts = 10 MHz counter locked to the KCU clock | 100 ns (mod 25.6 us) | 53.4k clusters at >= 70% per-match purity (43k true) | 20.3 sigma vs 2.8 sigma shuffled |

Numbers are from the 2026-10-21 re-run with the corrected decoders (10-bit TOT, out-of-order lines,
ROC_pool stitching: 3.45M complete KCU0 waveforms instead of 3.00M, 1.39M KCU1 instead of 1.16M). The
pre-fix outputs are kept in `../output/pre_fix_2026-10-21/`.

Only the stock `online_distribution` decoders are used, plus Purschke's `ROC_pool` for the KCU streams
(both from the forks in `../online_distribution` and `../ROC_pool`, branch `tanner`, built into
`$ONLINE_MAIN`; see `../decoder_tests/README.md` for what was changed there and why).

Synchronization review (2026-10-21): the KCU1 and CAEN links are final. For AstroPix the on-chip
`ts` phase is the right and only sub-ms handle, it is identical on all nine chips (no per-chip
alignment exists or is needed), and `fpga_ts` is a readout stamp that cannot pick the trigger on
its own. Details in `docs/method.md` 4.3; raw formats of all three streams, including the
HGCROC line bit map and the TOA/TOT fields, in `docs/data_format.md`.

## Status

`docs/status.md` is the overview of every study, result and open item (2026-11-03; section 14 indexes every deck and write-up and gives the plan and requests); start there. The current talk is `docs/slides_nov03/` (November 3, 2026 update deck, pptx + preview PDF, also an artifact).

## Layout

One folder per analysis thread. Each thread folder has an `out -> /sphenix/user/tmengel/eic/output` link,
so scripts that write `out/...` work from any of them; heavy jobs sit in the thread's own `condor/`.

```
setup.sh            environment (sPHENIX stack + ONLINE_MAIN=/eic/data/purschke/install)
common/             shared by several threads
  plotstyle.py      sPHENIX-style matplotlib setup (imported via sys.path.insert(0, '<this folder>'))
  astro_geo_common.py  AstroPix / SFIL geometry helpers
  bic_streams.{py,C}   multi-stream reader + matching for the ROOT trees; cache/ = its npz cache
  channel_map_fig.py   KCU0 channel-map figure
decoding/           raw .evt -> per-trigger ROOT trees (all runs) and their checks
  {kcu,caen,astro}_to_root.cc, kcu_merge.cc, kcu_hists.cc   (make all)
  kcu0_root_check.py, kcu_hists_plot.py, newrun_qa.py
  rawcheck/         raw-byte format checks (diagnostics only; see its README)
  condor/           run_to_root / run_merge jobs (+ the run-3062 kcu/caen/astro jobs), run_watch.sh
sync/               stream synchronization: the run-3062 pipeline and the run-generic sync
  reproduce.sh      run-3062 stages: build dump prep caen astro events plots diag
  ts_dump.cc, wf_acc.cc   stock-decoder dumps and weighted waveform sums (make all)
  *.py              alignment, matching, event tree, figures (see the stage table below)
  run_sync.py, run_decode_stats.py, run_sync_figs.py   any run, from the decoding/ trees; condor/run_sync*
  note/             LaTeX note on the sync of runs 3061-3091
  out/              all outputs (symlink to /sphenix/user/tmengel/eic/output, 5.7 GB; reproduce.sh recreates it)
waveforms/          KCU0 SFIL / KCU1 SFIL5 pulse fits, templates and galleries
  kcu0_sfil_{fit,plots}.C, kcu1_sfil5_{fit,plots,raw}.C, wf_overlay.py; kcu1_tot/ = SFIL5 saturation and TOT
  condor/           full-run fits (kcu0_fit, kcu1_fit)
timing/             TOA / TDC timing: toa_timing.py, toa_selfcal.py, toa_sfil2/, toa_clean/, toa_position/
calibration/        caen_calib/ (Baby BCAL), cells/, paddles/, attenuation/, rate/, shower.py; condor/caen_apply
astropix/           imaging-layer performance and ToT, AstroPix-based geometry (astro_geo_*, sfil5_map.py)
calibration/stability, cosmics/   run-3064 gain stability; cosmic MIP scale (3088) with the Geant4 cosmic mode
sim/                Geant4 model of the stack (BicGeometry.hh, v1-v4b), condor scans, analysis/ (shares, gains, performance 3062 / 3064)
dual_sided/         the 2017 dual-sided readout study repeated on run 3062 (note/)
event_displays/     synchronized event displays (evdisp.py, evdisp3d.py)
docs/slides_nov03/  the November 3, 2026 update deck (build_deck.py, figs/, artifact/)
dual_sided/         two-ended (N/S) timing and asymmetry study and its LaTeX note
event_displays/     evdisp.py (2D), evdisp3d.py + template (3D), and dated display sets
sim/                Geant4 model of the test-beam setup (SFIL1-5, AstroPix, Baby BCAL, paddles); see sim/README.md
  bic_sim.cc, src/, include/, macros/, CMakeLists.txt   (cd build; cmake ..; make -> build/bic_sim)
  condor/           energy / position scans (scan.job and variants); outputs in /sphenix/tg/tg01/jets/tmengel/eic/sim
  analysis/         sim summaries, data-vs-sim shares, per-layer and Baby BCAL gain calibration, photostatistics
maps/               channel map: BicMap.h, map_run3062.json
docs/
  status.md         overview of every study, result and open item; start there
  method.md         alignment and matching method, cuts and results, stream by stream
  channel_map.md    detector to KCU/ASIC/channel assignment and its sources
  pitfalls.md       decoder and timing pitfalls
  data_format.md    raw formats: KCU line bit map, TOA/TOT, CAEN words, AstroPix records, clocks
  new_runs.md       catalog of runs 3061-3091
data/
  raw/              .evt files (symlink to /sphenix/user/tmengel/eic/data, not in git)
  events_3062.root  full per-trigger event tree (925 MB, not in git)
  events_3062_skim.root  reduced tree in git: three-stream or AstroPix-matched triggers
```

## Running

```
source setup.sh
cd sync
./reproduce.sh                 # all stages, about 2 h on 8 cores
./reproduce.sh plots diag      # or selected stages
```

`reproduce.sh` reads the raw files from `data/raw/` (-> `/sphenix/user/tmengel/eic/data`) and writes to
`sync/out/` (-> `/sphenix/user/tmengel/eic/output`). Override with `RAW=... OUT=... ./reproduce.sh`.

| Stage | Scripts | Main output |
| --- | --- | --- |
| build | `make ts_dump wf_acc` | binaries |
| dump | `ts_dump` (kcu, kcuwf, kcuadc, kcu1adc, caen, caenadc, astro) | per-stream timestamps and amplitudes |
| prep | `prep_kcu0.py`, `prep_astro.py` | KCU0 timestamps, unwrapped AstroPix stamps |
| caen | `caen_kcu0.py` | CAEN on the KCU0 clock |
| astro | `astro_clusters.py`, `prep_firsthit.py`, `astro_ts10.py`, `astro_latency.py`, `astro_pred.py`, `astro_physics.py`, `trigger_table.py`, `astro_matcher.py`, `astro_final.py` | common trigger table, AstroPix matches |
| events | `event_tree.py`, `skim_events.py` | `out/events_3062.root`, `out/events_3062_skim.root` |
| plots | `sync_banger.py`, `astropix/astro_geo_kcu0.py`, `astropix/astro_geo_kcu1_caen.py`, `astro_summary_plot.py`, `astropix/sfil5_map.py`, `kcu1_align.py`, `paddle_trigger.py`, `sumcheck.py` | summary and geometry figures |
| diag | `sync_diag.py`, `astro_phase_split.py`, `calibration/shower.py`, `wf_keys.py`, `wf_acc`, `wf_bkgsub.py`, `event_displays/evdisp.py`, `event_displays/evdisp3d.py`, `slide_panels.py` | diagnostics, waveforms, event displays |

### What is needed for what

| Purpose | Files | Can be skipped if you only want the sync? |
| --- | --- | --- |
| Synchronization proper (KCU0 clock, KCU1 offset, CAEN clock fit, AstroPix phase + matching, trigger table) | `Makefile`, `ts_dump.cc` (modes kcu, caen, astro, kcuadc; links `libROC_pool` + `libEvent` from `$ONLINE_MAIN`), `prep_kcu0.py`, `prep_astro.py`, `caen_kcu0.py`, `astro_clusters.py`, `prep_firsthit.py`, `astro_pred.py`, `astro_physics.py` (writes `kcu0_amp.npz` for the likelihood), `trigger_table.py`, `astro_matcher.py`, `astro_final.py`, plus the modules `plotstyle.py`, `astro_geo_common.py` | no |
| Event tree (the usable output) | `event_tree.py`, `skim_events.py`, `ts_dump` modes kcu1adc, caenadc (amplitude branches) | no, unless you only need timestamps |
| Measuring the run constants for a new run | `kcu1_align.py` (KCU1 offset), `astro_onoff.py` (-0.43 s), `astro_ts10.py` (phase 189), `astro_latency.py` (latency windows), `astro_phase_split.py` (per-chip check); `caen_kcu0.py` measures its own | for run 3062, yes (constants are hard-coded) |
| Figures and physics checks | `sync_banger.py`, `sync_diag.py`, `astropix/astro_geo_*.py`, `astro_summary_plot.py`, `astropix/sfil5_map.py`, `paddle_trigger.py` (+ `ts_dump kcuwf`), `sumcheck.py`, `calibration/shower.py`, `wf_acc.cc`, `wf_keys.py`, `wf_bkgsub.py`, `event_displays/evdisp.py`, `event_displays/evdisp3d.py`, `slide_panels.py`, `k1empty_wf.py` (raw-waveform check of the KCU1-empty band, method.md 4.7) | yes |
| Format diagnostics | `decoding/rawcheck/`, `ts_dump kcutot` | yes |

Every file in `sync/` is used by one of these; there is no dead code left after the 2026-10-21 cleanup
(two unread outputs, `astro_best.npz` and `astro_cand_amp_idx.npz`, were dropped).

Shared modules: `plotstyle.py` (sPHENIX ROOT style, label and legend helpers) and
`astro_geo_common.py` (AstroPix cluster geometry). `astro_onoff.py` measures the -0.43 s beam on/off
correction to the AstroPix coarse offset.

### Waveform figures and the KCU0 channel map (2026-10-21)

`decoding/kcu_hists` fills the histograms of the LED macro `PlotEventBICSingle.C` for all 288 KCU0
channels (stock decoder + ROC_pool); `decoding/kcu_hists_plot.py` draws them as 4x4 / 2x2 panels per
detector group (`sync/out/kcuhists/{SFIL_South,SFIL_North,paddles}_{ADC,TOA,TOT,ADCcorr,ADCpulse}.png`,
`ASIC_<q>.png`); `waveforms/wf_overlay.py` overlays 30 individual waveforms per channel
(`wf_overlay_{pulses,random}.png`, from the `kcutime` file); `common/channel_map_fig.py` draws the KCU0
channel map along the beam line (`docs/kcu0_channel_map.png`). What they show: `docs/method.md` 4.10.

### Raw data in ROOT trees (2026-10-22)

`/sphenix/tg/tg01/jets/tmengel/eic/{kcu0,kcu1,caen,astropix}/` (symlinks `data/kcu0_3062.root`,
`kcu1_3062.root`, `caen_3062.root`, `astropix_3062.root`): every trigger of every stream, one entry per
trigger in stream order, organised by subsystem with every raw field, error flag and counter kept.
Makers `decoding/kcu_to_root.cc`, `caen_to_root.cc`, `astro_to_root.cc`, boundary step `kcu_merge.cc`;
condor jobs in `sync/condor/`; format and conventions in `docs/root_files.md`.

### Playing with the streams (`common/bic_streams.py`, 2026-10-22)

Importable module over the four stream files: `Streams().index()` loads the per-trigger clocks of the
whole run (cached in `common/cache/`), `match_kcu1()` pairs KCU0 and KCU1 by the exact tick offset
(927k pairs), `match_caen()` pairs CAEN triggers with the linear clock model plus a two-pass running
median (59% of CAEN triggers within 100 ns, MAD 6.4 ns), `astro_phase()` is the AstroPix on-chip phase
test, `event(i)` / `sfil_waveforms(i, layer)` fetch one trigger with its partners, and `pedestal`,
`pulse_features`, `fit_pulse` are the waveform tools (pedestal from samples 0-3, crossing time at a
fraction, saturation mask, a rise/decay pulse model fit). `python3 common/bic_streams.py` runs a demo and
writes `common/demo_event.png`. The same module as a ROOT macro: `common/bic_streams.C`
(`.L bic_streams.C+`, class `BicStreams`, `PulseFeatures`, `FitPulse`, `ToaTime`; `bic_streams()` runs the
demo and writes `demo_event_root.png`; the clock index is cached in `common/cache/index.root`).

### Subsystems

`docs/subsystems.md`: what each detector is and what it reads out (paddles, SFIL1-4, SFIL5, Baby BCAL,
AstroPix, DAQ and clocks).

### Beam conditions

`docs/beam.md`: species, energy and rate of the Hall D positron beam in run 3062 (pair spectrometer
arm, 3-6 GeV, no per-event energy in August, 4.10 kHz paddle-OR trigger), with `docs/beam_rate.png`.

## Event tree

`events_3062.root`, tree `events`: 4,533,542 entries (2026-10-21 re-run; 4,189,768 pre-fix), one per
trigger of the common paddle trigger.

| Branch | Content |
| --- | --- |
| `trig_ts`, `t_s` | trigger time (KCU0 ticks, seconds) |
| `has_kcu0`, `has_kcu1`, `has_caen` | which streams recorded the trigger |
| `pad`, `pad_kcu0`, `pad_caen` | paddle pattern (bits P3, P4, P5, P6) |
| `sfil_s_w`, `sfil_n_w` [4][4] | SFIL1-4 South / North, [cell][layer], pedestal and common-mode subtracted |
| `pad_dip`[4] | KCU0 paddle dip depth |
| `sfil5_s_w`[72] | KCU1 ch 144-215 (SFIL5 South arrays A-D) |
| `caen_q`[16] | CAEN pulse integrals: Baby BCAL South ch 0-11, paddles ch 12-15 |
| `sfil_s`, `sfil_n`, `sfil5_s`, `caen` | max minus pedestal amplitudes (for saturation flags) |
| `astro_x`, `astro_y`, `astro_chip`, `astro_log10LR`, `astro_purity`, `astro_latency_ms`, ... | best AstroPix match |

A stream that missed the trigger has `has_<stream> = 0` and zeros in its branches.

`events_3062_skim.root` has the same tree and branches, restricted to the 327,207 triggers recorded by
KCU0, KCU1 and CAEN (214k) or carrying an AstroPix match with log10 LR > 1 (189k); 393,799 entries
after the 2026-10-21 re-run (155k / 178k pre-fix). It is the copy
tracked in git; the full tree is too large for GitHub.

```
root -l data/events_3062_skim.root
events->Draw("sfil_s_w[2][1]", "has_kcu0")
```

## Run constants

These are measured for run 3062 and must be re-measured for other runs (3063, 3064):

- KCU1 offset 1,143,665,946 ticks, ts(KCU0) = ts(KCU1) + offset: measured by `kcu1_align.py`, hard-coded in
  `trigger_table.py`, `sync_banger.py`, `sync_diag.py`, `sumcheck.py`
- KCU0 RTCLOCK epoch 804631.17857 s (`caen_kcu0.py`)
- AstroPix coarse offset correction -0.43 s (measured by `astro_onoff.py`, hard-coded in the `astro_*.py`
  scripts and `sync_diag.py`), 10 MHz phase 189, latency windows 7-26 / 76-97 ms
