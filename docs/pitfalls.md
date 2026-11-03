# Pitfalls


- **Stock KCU decoder fragments waveforms.** Use only `SAMPLESIZE == 20` waveforms for timing; a
  fragment's first timestamp is shifted by k x 164 ticks. Drop the first stale `ts~0` buffer flush.
- **RTCLOCK nearest-neighbour matching is meaningless at ~kHz rates.** A neighbour within 1 ms is
  always there; a 10 s shift control gives 63% against 69%.
- **CAEN `EVENTTIME` is the wrong counter.** Use `GROUPTRIGGERTIME`.
- **AstroPix fpga_ts coincidence searches fail at every bin size,** because of the 19/88 +- 5 ms
  latency. Use the on-chip `ts` phase instead. The phase (189) is the same on all nine chips, for
  row and column records and over the whole run (`docs/astro_phase_split_3062.txt`), so there is no
  per-chip offset to apply.
- **The stock KCU decoder keeps only 6 of the 10 TOT bits** (`(w >> 10) & 0x3f`). The full code is
  bits 19:10 in the HGCROC3 compressed form and does track pulse size; the 6-bit value is a sawtooth.
  See `docs/data_format.md` and `decoding/rawcheck/kcu3.cc`. **Fixed 2026-10-21 in the local build**
  (`online_distribution` fork, `$ONLINE_MAIN=/sphenix/user/tmengel/eic/install`); `output/` was dumped
  with the old library. Verification and the full decoder assessment: `../decoder_tests/README.md`.
- **Out-of-order lines fragment waveforms** in Purschke's decoder and in ROC_pool (unsigned timestamp
  differences; a link runs a sample ahead now and then, on both boards). **Fixed in the local build
  2026-10-21** (timestamp-keyed sample assembly): complete waveforms per Packet 61% -> 75% on KCU0,
  48% -> 79% on KCU1, every well-formed old waveform reproduced identically. ROC_pool got the same
  fix and `ts_dump` / `wf_acc` read the KCU streams through it, which also removes the Packet-boundary
  cuts: 93% complete on both boards. See `../decoder_tests/README.md`.
- **TOA has no fine time without per-channel calibration** (MAD 9.6 ns against the ADC 50% crossing;
  uniform random would give 6.25 ns). Use the ADC 50% crossing.
- **CAEN group 1's trigger-time word is the high bits of group 0's**, not a second clock; after the
  stock decoder, `lValue(1,"GROUPTRIGGERTIME")` is 0 by design.
- **The AstroPix decoder's sort is a no-op.** `oncsSub_idastrov1::CompareByTimestamp` never compares
  fpga_ts, so hits come back in arrival order.
- **Beam-off seconds must come from CAEN.** "KCU0 silent" seconds are mostly KCU0 DAQ dropouts.
- **Summed-amplitude correlations are inflated by on/off.** KCU0-vs-KCU1 r = 0.78 inclusive, but only
  0.15 with both boards active.
- **Paddles on KCU0 saturate negative to ADC 0.** Do not use max - pedestal for them. Use the
  zero-run hit flag in 3.2.
- **Pedestal 1023 with MAD 0 marks a railed 36-channel group** (Kim deck slide 23's "high pedestal").
  Check per-group pedestal medians first on any new run.
- **An earlier timing study (`roc_timing/`, deleted 2026-10-21) used a custom decoder.** Its timing results agree with the stock-decoder pipeline,
  but the analyst asked for stock decoders only. It is not part of this repository.


- **`pad` in `events_3062.root` is a bitmask, not a paddle number** (P3 = 1, P4 = 2, P5 = 4, P6 = 8; `evdisp.py` PBIT,
  `shower.py` PADS). `pad == 5` is P3+P5, `pad == 6` is P4+P5. The 2026-10-22 South-gain calibration and the
  `data_shares.npz` / `data_p5_south.npz` per-paddle entries used `pad == 3..6` and were wrong; corrected 2026-10-23
  (old files kept as `*_WRONGPAD_2026-10-22.*`). The fit tree `fits` of `kcu0_sfil_fit.C` stores the paddle NUMBER
  (3-6) in `pad`: the two conventions differ.

- **SFIL4 North cells 1 and 2 (KCU0 ch 187, 186) carry a constant ~130 ADC "pulse" in every shower trigger**, independent
  of the shower size (SFIL3's neighbours scale with the shower), and SFIL4 North cell 3 (ch 144) sees 0.35 of the
  South cell-3 amplitude. Do not use the SFIL4 North end for energy or N/S studies until the hookup is understood;
  the layer-sum N/S ratio 0.64 of SFIL4 is made of those two channels. Found by the photostatistics study (method 4.13).

- **bic_sim cell assignment bug (fixed 2026-10-23 evening).** The fibre copy number packed `row * 1000 + fibre` under the cell
  digit; rows 10-16 of every layer overflowed into it, so 41 % of the fibres were booked one cell too low (towards cell 4).
  Every cell-level simulation result before the fix is wrong: the "beam straddles two cells", the "simulation puts more light
  in the off-row cells than the data, so the real cells are larger", the per-cell gains of 2026-10-22 and the cell-level
  double-ratio predictions. Layer sums (all shares, layer gains, the photostatistics comparison) were not affected.
  v1 and v2 outputs in `sim/v1_purePb/`, `sim/v2_slots17mm/` carry the bug; the v3 scan is the first clean one.

- **Cell 4 of every SFIL layer is empty in the data** (0-1 % of the layer sum at both ends, all paddles; channels S 14/33/48/4,
  N 158/177/192/148), while the model puts 5-11 % of a layer's light there for a beam 5 mm inside cell 3. Either those channels
  are dead or the beam sits higher in cell 3 than the official offset. Cells 1+2 carry one signal on the v2 summing boards, so
  summing both channels may double count that row. Lateral (cell-level) comparisons are therefore not yet meaningful; layer
  sums are affected at the few-percent level (drop cell 4 from the model when comparing).

- **SFIL5 South single-SiPM channels saturate too** (2026-10-24): the pedestal-subtracted peak of KCU1 channels pins at ~943 ADC
  (1023 minus a pedestal of ~80); the 90th percentile of pulses above 60 ADC sits there. The earlier statement that they never saturate
  came from testing `sfil5_s >= 1000`, which a pedestal-subtracted peak can never reach. All SFIL5 sums from raw peaks
  (`events_3062.root`, `data_sfil5_bbcal.npz`, the first `sfil5_sfil4_gain.py` numbers) are clipped; use the KCU1 template fits
  (`kcu1_sfil5_fit.C`). KCU1 pulses peak one sample earlier than KCU0 (samples 6-7): the KCU0 peak window 7-9 must not be used for KCU1.
- **KCU0 TOT is unusable for saturation recovery** (2026-10-24): on SFIL3 South cell 3 a TOT word accompanies only 0.5-1.6 % of saturated
  pulses; on KCU1 it accompanies 67-89 % of them and 17 % of unsaturated pulses above 900 ADC (calibration overlap). Use TOT on KCU1 only.
- **The CAEN DT5742 waveforms are not cell-calibrated: the "noise" is the DRS4 cell pedestal pattern** (2026-10-24,
  `calibration/caen_calib/`). Sample s of an event with group start cell c0 sits in physical cell (c0 + s) mod 1024, and
  each cell has its own offset (cell-to-cell spread 32-34 counts on every channel, within-cell rms 9). The raw
  pedestal-subtracted waveforms therefore show 26-28 counts rms "noise", a start-cell dependent baseline step
  between the first and last samples (+30 counts for start cells 256-511, -30 for 768-1023) and, because the pattern is
  the same on all 16 channels of the board (corr 0.95 between Baby BCAL and idle paddle channels, 0.96 between the two
  DRS4 groups), an apparent common-mode excursion: in SFIL-quiet triggers all twelve Baby BCAL channels "fire" at once
  with a flat ~2000 units of charge, which is just a 15-30 count offset integrated over 140 samples. A per-channel
  per-cell pedestal table from the pre-pulse samples of 150k events (`drs4_pedestals.py`) reduces the rms from 28 to
  9 counts on the Baby BCAL channels and removes the step; `caen_apply.py` (condor, one job per segment) rewrites the
  per-event pulse quantities with it. The event-table `caen_q` / `caen` branches are the uncorrected values: the small
  Baby BCAL channels (A3, A4: 530-560 units per shower) and their 0.6-0.8 correlations are mostly this artefact.
- **KCU1 single-SiPM pulses change shape with amplitude and flat-top when clipped** (2026-10-24, `waveforms/kcu1_tot/`): the
  tail/peak ratios fall steadily from 100 to 940 ADC and the sample after a clipped stretch is ~750 ADC for any charge, so a
  template fit with masked samples returns LESS than the clipped peak and no ADC-tail method works; use the decoded TOT
  (time over the ceiling, 166 units per 25 ns) with A = 943 + b (TOT - 419), b = 0.40 (0-0.76). KCU0 summed channels keep
  a fixed shape below 800 ADC. The earlier "SFIL5 nonlinear by a factor 4.5" reading of A5/A4 versus A4 was shower physics:
  the model's <L5|L4> falls the same way.
- **CAEN paddle spectra are threshold-truncated at the MIP peak** (2026-10-24): the single-paddle pulse-height spectra rise
  sharply at 1325 / 1275 / 1925 / 2150 counts (P3 / P4 / P5 / P6) and peak 10 % above, with 14-24 % of triggers below the
  edge at the same time; a Landau fit "MPV" is just the edge. Do not read relative paddle gains from the peaks alone, and do
  not use other-paddle-triggered events as an unbiased sample (they are backsplash). P5/P6 pulses are 1.45-1.62 x P3/P4 in the
  CAEN, thresholds included.
- **AstroPix chips 5-8 have no usable ToT in run 3062** (2026-10-24): 38-79 % of their column records carry ToT = 0 and the
  rest are a few units; chips 0-4 are fine (median 400-830 units). Hit positions and timestamps of chips 5-8 are unaffected.
  Records with a non-zero upper ToT nibble (8 %) have shifted ToT medians and are excluded until the A-STEP hit format
  settles what the nibble means.
- **Baby BCAL CAEN labels: use the run-log sheet, not the PDF snapshot** (2026-10-24): the sheet's "Changed on Aug 10" table
  maps CAEN 0-3 = A0 A1 A3 A4, 4 = unlabeled (A2 by its pattern), 5 = A5, 6-9 = A6-A9, 10 = A14, 11 = A11. The PDF text
  extraction dropped the "4" and shifted channels 4-10 one label too low, which produced an apparent row-2 anomaly.
