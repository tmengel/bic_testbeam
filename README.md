# BIC test beam, run 3062: stream synchronization

Analysis of the Barrel Imaging Calorimeter (BIC) test beam at JLab Hall D (August 2026), run 3062.
The four DAQ streams are put on one clock, the KCU0 FPGA timestamp (6.25 ns ticks), and combined
into one per-trigger event tree.

| Stream | Link to KCU0 | Precision | Matched | Null control |
| --- | --- | --- | --- | --- |
| KCU1 (H2GCROC, SFIL5 South) | fixed offset -1,143,665,946 ticks | 0 or +25 ns | 53% of KCU1 triggers | 0.0% (+10 ms shift) |
| CAEN DT5742 | clock ratio 1.365340679 + running-median correction | residual MAD 6.4 ns | 49.7% of 1.51M | 0.0% (+10 ms shift) |
| AstroPix 9-chip module | on-chip ts = 10 MHz counter locked to the KCU clock | 100 ns (mod 25.6 us) | 50.8k clusters at 81% purity | 19.3 sigma vs 3.0 sigma shuffled |

Only the stock `online_distribution` decoders are used.

## Layout

```
setup.sh            environment (sPHENIX stack + ONLINE_MAIN=/eic/data/purschke/install)
sync/               the pipeline (run everything from here with reproduce.sh)
  reproduce.sh      stages: build dump prep caen astro events plots diag
  ts_dump.cc        stock-decoder dump of timestamps and amplitudes (KCU, CAEN, AstroPix)
  wf_acc.cc         weighted waveform sums for background subtraction
  *.py              alignment, matching, event tree, figures (see the stage table below)
  out/              all outputs (created by reproduce.sh)
maps/               channel map: BicMap.h, map_run3062.json
docs/
  method.md         alignment and matching method, cuts and results, stream by stream
  channel_map.md    detector to KCU/ASIC/channel assignment and its sources
  pitfalls.md       decoder and timing pitfalls
data/
  raw/              run 3062 .evt files (link to the raw data area, 36 GB, not in git)
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

`reproduce.sh` reads the raw files from `data/raw/`. Point that link (or `RAW` in the script) at
another copy of the run 3062 files if needed.

| Stage | Scripts | Main output |
| --- | --- | --- |
| build | `make ts_dump wf_acc` | binaries |
| dump | `ts_dump` (kcu, kcuwf, kcuadc, kcu1adc, caen, caenadc, astro) | per-stream timestamps and amplitudes |
| prep | `prep_kcu0.py`, `prep_astro.py` | KCU0 timestamps, unwrapped AstroPix stamps |
| caen | `caen_kcu0.py` | CAEN on the KCU0 clock |
| astro | `astro_clusters.py`, `prep_firsthit.py`, `astro_ts10.py`, `astro_latency.py`, `astro_pred.py`, `astro_physics.py`, `trigger_table.py`, `astro_matcher.py`, `astro_final.py` | common trigger table, AstroPix matches |
| events | `event_tree.py`, `skim_events.py` | `out/events_3062.root`, `out/events_3062_skim.root` |
| plots | `sync_banger.py`, `astro_geo_kcu0.py`, `astro_geo_kcu1_caen.py`, `astro_summary_plot.py`, `sfil5_map.py`, `kcu1_align.py`, `paddle_trigger.py`, `sumcheck.py` | summary and geometry figures |
| diag | `sync_diag.py`, `shower.py`, `wf_keys.py`, `wf_acc`, `wf_bkgsub.py`, `evdisp.py`, `evdisp3d.py`, `slide_panels.py` | diagnostics, waveforms, event displays |

Shared modules: `plotstyle.py` (sPHENIX ROOT style, label and legend helpers) and
`astro_geo_common.py` (AstroPix cluster geometry). `astro_onoff.py` measures the -0.43 s beam on/off
correction to the AstroPix coarse offset.

## Event tree

`events_3062.root`, tree `events`: 4,189,768 entries, one per trigger of the common paddle trigger.

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
KCU0, KCU1 and CAEN (155k) or carrying an AstroPix match with log10 LR > 1 (178k). It is the copy
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
