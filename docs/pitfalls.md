# Pitfalls


- **Stock KCU decoder fragments waveforms.** Use only `SAMPLESIZE == 20` waveforms for timing; a
  fragment's first timestamp is shifted by k x 164 ticks. Drop the first stale `ts~0` buffer flush.
- **RTCLOCK nearest-neighbour matching is meaningless at ~kHz rates.** A neighbour within 1 ms is
  always there; a 10 s shift control gives 63% against 69%.
- **CAEN `EVENTTIME` is the wrong counter.** Use `GROUPTRIGGERTIME`.
- **AstroPix fpga_ts coincidence searches fail at every bin size,** because of the 19/88 +- 5 ms
  latency. Use the on-chip `ts` phase instead.
- **The AstroPix decoder's sort is a no-op.** `oncsSub_idastrov1::CompareByTimestamp` never compares
  fpga_ts, so hits come back in arrival order.
- **Beam-off seconds must come from CAEN.** "KCU0 silent" seconds are mostly KCU0 DAQ dropouts.
- **Summed-amplitude correlations are inflated by on/off.** KCU0-vs-KCU1 r = 0.78 inclusive, but only
  0.15 with both boards active.
- **Paddles on KCU0 saturate negative to ADC 0.** Do not use max - pedestal for them. Use the
  zero-run hit flag in 3.2.
- **Pedestal 1023 with MAD 0 marks a railed 36-channel group** (Kim deck slide 23's "high pedestal").
  Check per-group pedestal medians first on any new run.
- **An earlier timing study (`roc_timing/`) used a custom decoder.** Its timing results agree with the stock-decoder pipeline,
  but the analyst asked for stock decoders only. It is not part of this repository.

