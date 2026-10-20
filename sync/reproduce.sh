#!/bin/bash
# Reproduce the run-3062 stream synchronization (KCU0 reference; KCU1, CAEN, AstroPix).
# Stock online_distribution decoders only.
#
#   ./reproduce.sh [stage ...]      stages: build dump prep caen astro events plots diag   (default: all)
#
# Wall time on an 8-core node: dump ~25 min (one process per segment), rest ~1 h.
# Run-specific constants (re-measure for another run, see the write-up):
#   caen_kcu0.py      Lk = KCU0 RTCLOCK epoch (804631.17857 s)
#   astro_*.py        O  = envelope offset - 0.43 s (beam on/off correction, astro_onoff.py)
#                     phase 189, latency windows 12-24 / 82-94 ms
#   sync_banger.py, astro_geo_kcu1_caen.py   KCU1 offset -1143665946 ticks
set -e
cd "$(dirname "$0")"
source ../setup.sh >/dev/null 2>&1
RAW=../data/raw
RUN=00003062
STAGES=${@:-build dump prep caen astro events plots diag}
mkdir -p out

waitall() { for p in $(jobs -p); do wait $p; done; }

for s in $STAGES; do
  echo "=== stage $s  ($(date +%T))"
  case $s in
  build)
    make ts_dump wf_acc ;;
  dump)
    # timestamps: one row per trigger (KCU, CAEN) or per hit (AstroPix), plus RTCLOCK
    for f in $RAW/beam_ROC_0-$RUN-*.evt $RAW/beam_ROC_1-$RUN-*.evt; do
      ./ts_dump kcu out/$(basename $f .evt).root $f > out/$(basename $f .evt).log 2>&1 &
    done
    ./ts_dump caen  out/caen_3062.root  $RAW/beam_CAEN-$RUN-000{0,1,2}.evt > out/caen.log 2>&1 &
    ./ts_dump astro out/astro_3062.root $RAW/beam_AstroPix-$RUN-0000.evt   > out/astro.log 2>&1 &
    ./ts_dump kcuwf out/kcuwf_0001.root $RAW/beam_ROC_0-$RUN-0001.evt > out/kcuwf.log 2>&1 &   # raw paddle waveforms
    waitall
    # amplitudes: KCU0 SFIL1-4 + paddles, KCU1 SFIL5 South (72 ch), CAEN ch 0-15
    for f in $RAW/beam_ROC_0-$RUN-*.evt; do ./ts_dump kcuadc  out/adc_$(basename $f .evt).root $f > out/adc_$(basename $f .evt).log 2>&1 & done
    for f in $RAW/beam_ROC_1-$RUN-*.evt; do ./ts_dump kcu1adc out/adc_$(basename $f .evt).root $f > out/adc_$(basename $f .evt).log 2>&1 & done
    for f in $RAW/beam_CAEN-$RUN-*.evt;  do ./ts_dump caenadc out/adc_$(basename $f .evt).root $f > out/adc_$(basename $f .evt).log 2>&1 & done
    waitall ;;
  prep)
    python3 prep_kcu0.py          # out/kcu0_ts.npy
    python3 prep_astro.py ;;      # out/astro_unwrapped.npz
  caen)
    python3 caen_kcu0.py ;;       # out/caen_kcu0_match.npz, out/caen_kcu0.png
  astro)
    python3 astro_clusters.py     # out/astro_clusters.npz (clusters + envelope offset)
    python3 prep_firsthit.py      # out/astro_cluster_firsthit.npy
    python3 astro_ts10.py         # out/astro_ts10_H.npy: 10 MHz phase lock vs shuffled ts
    python3 astro_latency.py      # out/astro_latency.npz: K - readout stamp, phase-matched vs control
    python3 astro_pred.py         # out/astro_attr.npz, out/astro_cand_bands.npz
    python3 astro_physics.py      # out/kcu0_amp.npz + true-vs-random observables
    python3 trigger_table.py      # out/triggers_3062.npz: common paddle-OR trigger table (KCU0+KCU1+CAEN)
    python3 astro_matcher.py      # out/astro_matcher_rows.npz (TRIG=kcu0: KCU0 complete waveforms only)
    python3 astro_final.py ;;     # out/astro_kcu0_matches.npz (purity table printed)
  events)
    python3 event_tree.py         # out/events_3062.root: one entry per trigger, all streams + AstroPix
    python3 skim_events.py ;;     # out/events_3062_skim.root: three-stream or AstroPix-matched triggers (in git)
  plots)
    python3 sync_banger.py        # out/sync_banger.{png,pdf}
    python3 astro_geo_kcu0.py     # out/astro_geo_kcu0.png
    python3 astro_geo_kcu1_caen.py     # out/astro_geo_kcu1_caen.png
    python3 astro_summary_plot.py      # out/astro_summary.png
    python3 sfil5_map.py               # out/sfil5_map.png: SFIL5 South vertical array order
    python3 kcu1_align.py              # out/kcu1_align.png: KCU1 -> KCU0 offset and residual
    python3 paddle_trigger.py          # out/paddle_trigger.png: paddle signal and KCU0/CAEN agreement
    python3 sumcheck.py ;;             # out/sumcheck.png: ADC sum check, max - ped vs window - ped - CM
  diag)
    python3 sync_diag.py          # out/sync_diag.png: every link vs run time
    python3 shower.py             # out/shower.png: shower shape, relative energy response per paddle
    python3 wf_keys.py            # out/wfkeys_*.txt: AstroPix phase-band triggers + event-display triggers
    mkdir -p out/wfacc            # waveform sums (stock decoders), ~25 min on 8 cores
    for f in $RAW/beam_ROC_0-$RUN-*.evt; do b=$(basename $f .evt); ./wf_acc kcu0 out/wfkeys_kcu0.txt out/wfacc/kcu0_${b##*-}.root $f > out/wfacc/kcu0_${b##*-}.log 2>&1 & done
    for f in $RAW/beam_ROC_1-$RUN-*.evt; do b=$(basename $f .evt); ./wf_acc kcu1 out/wfkeys_kcu1.txt out/wfacc/kcu1_${b##*-}.root $f > out/wfacc/kcu1_${b##*-}.log 2>&1 & done
    ./wf_acc caen out/wfkeys_caen.txt out/wfacc/caen.root $RAW/beam_CAEN-$RUN-000{0,1,2}.evt > out/wfacc/caen.log 2>&1 &
    waitall
    python3 wf_bkgsub.py          # out/wf_bkgsub.png: waveform-level background subtraction
    python3 evdisp.py             # out/evdisp/ev<row>.png, events.json
    python3 evdisp3d.py           # out/evdisp/bic_event_viewer.html (interactive, schematic geometry)
    python3 slide_panels.py ;;    # out/*_slide.png: reduced figures for the slides
  *) echo "unknown stage $s"; exit 1 ;;
  esac
done
echo "=== done ($(date +%T))"
