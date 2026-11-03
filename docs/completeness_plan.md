# Maximizing waveform / event completeness per stream (plan, 2026-10-21)

Where the losses are, measured on run 3062 (300-1500 RCDAQ packets per KCU, `decoder_tests/`), and
what to do about each. Rates are relative to the triggers the KCU FPGA counted (`trig_out`, ~4.15 kHz
in both boards), not to the 2.8 kHz union of recorded triggers.

## Loss budget

| Loss | KCU0 | KCU1 | CAEN | AstroPix | Recoverable offline? |
| --- | --- | --- | --- | --- | --- |
| Out-of-order lines split waveforms in the decoder | 14% of waveforms | 31% | - | - | **yes, done** (local decoder) |
| Waveforms cut at RCDAQ Packet boundaries | ~19% | ~19% | 0 (1 trigger = 1 packet) | - | **yes, done**: ROC_pool with the reorder fix, used by `ts_dump` / `wf_acc` (complete 93.4% / 92.8% of received) |
| UDP datagrams dropped mid-packet (gaps of 26-272 frames, median 136 = 186 kB) | ~23% of triggers, 0.52 gaps per packet over the segment | 0.53 gaps per packet over the segment (none in the first 300 packets) | - | - | no |
| Multi-second dropouts (whole stream silent) | occasional (1.7 s gap seen in 300 packets; t ~ 550 s outage) | dominant: duty cycle ~19% | 0 | readout windows of ~40 ms with 70-100 ms gaps | no |
| Single-sample / 2-sample readouts | ~5% of waveforms | ~2% | - | - | no (no pulse in them) |
| Dead time / rate limit | - | - | 1014 Hz recorded of ~4.15 kHz issued (DT5742 conversion + USB) | self-triggered, no dead time, but 8-bit timestamp | no |
| Trigger not stamped | - | - | - | latency 19 / 88 +- 5 ms: no per-trigger match | no (80% purity ceiling) |

Result for run 3062 after the offline steps: ~95% of the triggers each KCU *received* become complete
waveforms (from 60% / 48%); the triggers the DAQ never received stay lost, so the three-stream
coincidence stays at the few-% level set by the DAQ.

## Offline steps for run 3062 (in order)

1. **Re-dump with the local decoder and the pool** (both done 2026-10-21; `ONLINE_MAIN` is the local
   install): `cd bic_testbeam/sync && ./reproduce.sh dump prep caen astro events` (launched 2026-10-21
   17:58, log `output/reproduce_2026-10-21.log`; the pre-fix outputs are kept in
   `output/pre_fix_2026-10-21/`). Complete waveforms per received trigger: KCU0 61% -> 93%, KCU1
   48% -> 93%; more rows in `triggers_3062.npz`.
2. **ROC_pool reorder port** done (`h2gcroc3_10G_pool`: `_open` map, `close_samples()`, reset guard,
   CRC/hamming and error counters); `pool_validate` PASSED on both boards; `ts_dump` and `wf_acc` read
   the KCU streams through it (`h2gcroc3_10G_pool`, `addPacket`, `depth_ok()` loop, `drain()` at the end).
3. **Use partial waveforms where they are enough.** A fragment cut at its end still has pedestal and
   peak (samples 1-3 and 6-9) if it has >= 10 samples; a fragment cut at its start has an unknown
   offset (k x 164) and is useless for timing but its timestamp can still be matched to a known trigger
   from another stream (`kcu0_frag` in the trigger table does this). Lift `SAMPLESIZE == 20` to
   `>= 10 and first sample is a pedestal sample` in `kcuadc` for amplitudes only.
4. **Cut on integrity, not on luck:** require `CRCERROR == 0` for every channel used (only the dead
   KCU1 board fails), and record `HAMMINGERRORS` / `DECODEERRORS` per packet in the dumps.
5. **AstroPix:** nothing more to gain offline; the on-chip 8-bit 10 MHz timestamp plus the two
   latency populations is the information there is (81% purity at 50.8k clusters).

## DAQ steps for the next run (these are where the real completeness is)

1. **Common busy / veto.** The paddle OR fans out to the KCU trigger card and both CAENs with no shared
   busy, so every stream records a different subset (KCU0+KCU1+CAEN coincidence = 3.7% of the union).
   OR the busy outputs of both KCUs and the CAENs into the trigger veto, so that a trigger is issued only
   when every stream can take it. This alone makes the recorded trigger sets identical and the matching
   1:1.
2. **Stop the mid-packet UDP loss on KCU0.** The gaps are 26-272 datagrams (36-372 kB), i.e. about one
   default Linux UDP socket buffer (`net.core.rmem_default` = 212 kB), and they occur inside a buffer,
   0.6 times per 273 kB packet, at 126 MB/s of line data. Raise `net.core.rmem_max` and have the rcdaq
   plugin set `SO_RCVBUF` to tens of MB; read with `recvmmsg`; pin the KCU0 rcdaq process to its own
   core; consider jumbo frames (1366-byte datagrams of 7 lines are 92k datagrams/s per KCU). KCU1 on the
   same host shows no such loss while alive, so compare the two processes' settings first.
3. **Find the long dropouts.** KCU1 is silent ~80% of the time, in stretches of seconds; KCU0 has them
   too. Correlate with the rcdaq logs (buffer switches, file rollover on GPFS at 3 GB segments) and the
   host's disk writes; a local SSD spool or larger/more ring buffers is the usual cure. The frame
   sequence number (datagram bytes 8-11) is the monitor: a pmonitor plot of sequence gaps per second
   shows loss live (the raw pass in `decoder_tests/kcu_decoder_check.cc` is the reference code).
4. **CAEN dead time.** 1014 Hz recorded of ~4.15 kHz: shorten the record (256 -> 128 samples, or
   fewer groups) or move from USB to the optical link; with a common busy its dead time then gates all
   streams the same way instead of silently dropping its own.
5. **AstroPix trigger stamp.** Feed the paddle OR into a spare A-STEP input so the FPGA stamps every
   trigger with `fpga_ts` (or a sync pulse into both the KCU trigger card and the A-STEP, or a common
   counter reset). Then hits match triggers at 25 ns instead of through an 8-bit counter and a 19/88 ms
   latency model. Log the A-STEP timestamp-clock setting.
6. **Shorter RCDAQ event buffers or a trigger-aligned cut** would reduce the boundary fragments at
   the source, but with ROC_pool in place this is not needed.
7. **Keep `trig_out`.** It is the only count of issued triggers; use it to quote per-stream efficiency
   (recorded / issued) in every run summary.

## Expected effect

| Stage | KCU0 complete waveforms per issued trigger | KCU1 | three-stream coincidence |
| --- | --- | --- | --- |
| run 3062, Purschke's decoder | 0.53 x 0.61 = 32% | 0.19 x 0.48 = 9% | 3.7% of union |
| + local decoder (now) | 0.53 x 0.75 = 40% | 0.19 x 0.79 = 15% | up to ~5% |
| + ROC_pool with reorder fix (measured on the re-run: KCU0 complete 2322 Hz, KCU1 937 Hz, three-stream 213,754 = 4.7% of a 4.53M union) | 0.56 x 0.92 = 51% | 0.23 x 0.91 = 21% | 4.7% |
| next run with common busy, UDP loss and dropouts fixed | ~95% of issued | ~95% | ~90% (CAEN-rate limited) |
