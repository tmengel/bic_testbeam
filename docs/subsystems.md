# The subsystems of the August 2026 BIC test beam and what each one reads out (run 3062)

Beam enters from the left: trigger paddles, AstroPix module, SFIL1-4, SFIL5, Baby BCAL (`docs/kcu0_channel_map.png`).
"North" and "South" are the two fibre ends of the Pb/SciFi layers; chip 0 of the AstroPix module is the
North end. Sources: the 8/13 hookup table (`maps/channel_mapping.txt`), the TIC deck, the Kim SiPM
deck (10/8 update), the run log, and the data (`docs/channel_map.md`, `docs/root_files.md`,
`docs/method.md`).

## 1. Trigger paddles (P3, P4, P5, P6)

- **What:** four plastic scintillator counters on Hamamatsu R580 photomultipliers, mounted along the
  AstroPix module in front of the calorimeter stack. From North to South: P4 over chip 0, P3 over
  chips 2-3, P6 over chips 5-6, P5 over chip 8. Their sizes are not in the notes (the simulation
  assumes 20 x 60 x 10 mm).
- **Role:** the trigger. The four PMT signals go through NIM modules; their OR is the trigger sent
  to both KCUs and to the CAEN (no common busy, each board applies its own dead time), at 4.10 kHz
  in run 3062. In May the paddles were calibrated against the pair-spectrometer counters and tagged
  the positron energy (P4 < P3 < P6 < P5); in August they were not calibrated, so they only tag the
  position along the module.
- **Read out twice:** analog copies go to KCU0 ASIC1 channels 84 (P3), 86 (P4), 81 (P5), 97 (P6),
  20 samples of 25 ns, negative pulses that rail at ADC 0 for about 9 samples from sample 8-9; and to
  the CAEN DT5742 channels 12-15 (P3-P6), 256 samples of 1 ns, negative pulses on a pedestal of
  3720-3770 that reach 0 on large hits. 99.85% of KCU0 triggers show a paddle pulse, one paddle in
  88%, two in 9.5%; KCU0 and CAEN agree on the pattern in 96-99% of shared triggers.

## 2. SFIL1-4: Pb/SciFi intermediate layers with summing-board readout

- **What:** four "SciFi intermediate layers", each a slab of the GlueX-style Pb/SciFi matrix (1 mm
  scintillating fibres in corrugated lead, 17 fibre rows, about 1.4 X0 thick), fibres running
  North-South, stacked one behind the other along the beam at a 65 mm pitch (epic BeamTest branch), the AstroPix module in front of
  SFIL1. Light from each fibre end goes through a light guide to a
  4 x 4 array of 3 x 3 mm SiPMs (Hamamatsu S14161-3050); each layer end has four such arrays stacked
  vertically, the readout "cells" 1-4. The beam runs through cell 3 of every layer.
- **Summing boards:** on every end a passive board sums the 16 SiPMs of an array into one signal, so
  each layer gives 4 channels per end, 32 channels in all. The v2 boards (all South ends, North
  SFIL2-4) mix the SiPM1 and SiPM2 arrays, so cells 1 and 2 carry the same signal and count as one
  row; North SFIL1 is a v1 board with one-to-one cells. SFIL4 has a large North/South gain offset.
- **Readout:** H2GCROC3 ASICs on the KCU0 FPGA board, in characterization mode: per channel and per
  trigger 20 samples of 25 ns, each with a 10-bit ADC, a 10-bit TOA code (TDC, nominal 25 ps LSB,
  fires once per pulse crossing its threshold) and a 10-bit TOT code (fires only on pulses that
  saturate the ADC). South ends are on ASIC0 (ch 2/6/10/14, 21/25/29/33, 63/59/52/48, 43/37/0/4 for
  SFIL1-4 cells 1-4), North ends on ASIC2 (146/150/154/158, 165/169/173/177, 207/203/196/192,
  187/186/144/148). Pulses are positive, peak at sample 7-8, have a long tail and clip at 1023 in
  20-36% of all beam-row waveforms (53-66% of those with a shower, rising with the paddle energy P4 < P3 < P6 < P5); pedestals 60-105 ADC with a per-trigger common baseline shift of
  about 26 ADC. Of the TOA, only the paddles and one SFIL channel (196) emit words on their pulses
  in this run; the TOT appears only on saturated pulses. (`docs/method.md` 4.9, 4.10)

## 3. SFIL5: the same layer with single-SiPM readout

- **What:** a fifth Pb/SciFi layer of the same build, downstream of SFIL4, instrumented with the v1
  single-SiPM boards: each of the 64 SiPMs of the South end is its own channel (four 16-SiPM arrays
  A, B, C, D on KCU1 ASIC2, ch 144-160, 162-178, 180-196, 198-214 with the gaps 152, 161, ... not
  connected; 53 of 64 fire; vertical order from the data A < B < D < C, bottom to top). The North
  end board sits on KCU1 ASIC0/ASIC1 and is dead in this run: half 0 of both ASICs (ch 0-35, 72-107) is railed at 1023 in
  every trigger (inputs pegged, a bias or board fault, nothing to recover offline); half 1 has normal pedestals but only four
  SiPMs see light (ch 54, 60, 65, 66: pulses of 75-200 ADC in 22-40 % of South showers, each 1-2.5 % of the South total,
  correlation 0.3-0.5 with it; ch 60 has a pedestal at 0). No North layer sum exists; the four channels can only serve the
  TOA timing against the South end. In the run plan SFIL5 was meant for a Cherenkov comparison on the same layer.
- **Readout:** H2GCROC3 on the KCU1 FPGA, same format as KCU0: 20 samples of ADC / TOA / TOT per
  channel per trigger. The TOA works on all 16 live array channels tested (self-calibrated to about
  190 ps per channel), which is what the single-SiPM boards were for. KCU1 recorded about 1.0 kHz
  of the 4.1 kHz triggers (links and UDP loss).

## 4. Baby BCAL: the containment calorimeter

- **What:** the Baby BCAL of the 2024 FTBF test (arXiv:2504.03079, `context/`): a GlueX-type Pb/SciFi block,
  58 cm long fibres read out from both ends, SiPM face 11.8 x 22.2 cm, 40 SiPM arrays per side in 4 lateral
  columns x 10 rows in depth. The FTBF summing boards concatenate each column's 10 arrays into 4 channels of
  1, 2, 3 and 4 arrays (rows 0-3, upstream to downstream), so column c holds channels c, c+4, c+8, c+12 and
  the "A-number" labels in the run log are these 16 channel numbers.
- **Readout (run 3062):** the South end on the CAEN DT5742 desktop digitizer. The run-log sheet ("Changed on Aug 10",
  `context/runlog_sheet/`) maps CAEN ch 0 = A0, 1 = A1, 2 = A3, 3 = A4, 4 = unlabeled, 5 = A5, 6-9 = A6-A9, 10 = A14,
  11 = A11 (an earlier reading of the PDF snapshot had channels 4-10 one label too low). Channel 4 carries 8.5 % of the shower
  charge with the same pattern as A1 and sits where A2 (row 0, beam column) belongs, so it is taken as A2. A10, A12, A13 and
  A15 are not read; A3 and A4 (CAEN 2, 3) are dead. 12-bit DRS4 samples at 1 GS/s, 256 samples, negative pulses from
  sample ~20 on pedestals of 3450-3530; no channel ever reaches the ADC floor.
  The CAEN reads about 1 kHz of the triggers and stamps each with its group trigger time (8.533 ns ticks, the
  clock used for the synchronization). The North end is read as 40 single arrays through the active summing
  boards (BBS01-40 on J1-J40) on KCU0 ASIC3, ch 216-287, which carries no signal in run 3062 (ch 217-251 at a
  clean pedestal, 216 and 252-287 railed): no North-South geometric mean and no N-S photostatistics are possible.
- **Calibration state (2026-10-24, `status.md` 7):** no gain table exists for the August Baby BCAL. The FTBF
  paper gives the 2024 MIP peaks per channel for a different ADC (JLab fADC250), the epic `BabyBCal` macros only
  the digitization model; our own Baby BCAL gains are simulation-anchored totals. The CAEN waveforms need the
  per-cell DRS4 pedestal correction (`pitfalls.md`, `calibration/caen_calib/`) before any channel-level work: with
  it, **CAEN channels 2 and 3 (A3, A4) are dead** and the ten live channels share the shower charge as A5 27 %,
  A6 28 %, A9 10 %, A1 and A2 8.5 % each, A7 7 %, A11 3 %, A8 3 %, A0 2 %, A14 2 %, nearly independent of paddle: the
  beam sits between columns 1 and 2 in every row.
  No in-run MIP sample exists (no pulser, no muons, ~8 random cosmic coincidences in the whole run).

## 5. AstroPix: the imaging layer

- **What:** one 9-chip module of AstroPix v3 monolithic silicon pixel sensors (the ePIC BIC imaging
  layer technology), 35 x 35 pixels of 0.5 mm per chip, chips 20 mm apart along the module, placed
  in front of SFIL1 (the matched clusters are 1-2 pixels, a single track, which a module behind 1.4 X0 would not give). Chip 0 is the
  North end; rows run vertically and columns along the module, x = 20 mm x chip + 0.5 mm x column.
  The beam band sits at rows 12-22 and is densest on chips 3, 5, 6 and 8; columns 0-2 are nearly
  empty.
- **Readout:** self-triggered streaming through the A-STEP FPGA, independent of the paddle trigger.
  Every pixel over threshold produces a row record and a column record, each an 11-byte SPI record
  with the chip id, row-or-column flag, location, an 8-bit on-chip timestamp (10 MHz counter), a
  12-bit time-over-threshold and the 32-bit FPGA readout stamp (40 MHz). A pixel hit is a row record
  plus a column record with the same on-chip timestamp; the readout stamp lags the trigger by 19 or
  88 ms, so matching to the calorimeter goes through the on-chip timestamp phase and the paddle
  pattern (`method.md` 4.3). 2.86 M hit records in run 3062, about 1.9 kHz, 99% from layer 1.
- **ToT (2026-10-24, `status.md` 7):** 12-bit column and row ToT agree to 3 %; chips 0-4 give MIP spectra with
  Landau x Gauss MPVs of 676 / 484 / 585 / 724 / 331 ToT clock units (5 ns in Bobae Kim's decoded CSV, 10 ns in the public A-STEP decoder; the yaml decides), chips 5-8 read ToT 0 on most records. One-point energy
  scale E = 34 keV x ToT / MPV_pixel (`astropix/`); the module's source calibration is not on disk.

## 6. The DAQ and the clocks

- **RCDAQ** records the four streams as separate files (beam_ROC_0, beam_ROC_1, beam_CAEN,
  beam_AstroPix), each RCDAQ event stamped with the host RTCLOCK (ns). One KCU packet holds 200 UDP
  datagrams, about 8-9 triggers; UDP loss and the 100 us trigger busy are the main sources of
  missing triggers (`completeness_plan.md`).
- **Clocks:** one 40 MHz card feeds both KCUs (160 MHz timestamps, 6.25 ns ticks), the AstroPix FPGA
  (40 MHz readout stamp, 10 MHz on-chip counter) and the CAEN (117.1875 MHz group trigger time)
  without a common reset, so each stream has its own epoch: KCU1 - KCU0 = -1,143,665,946 ticks,
  CAEN 1.365340679 KCU ticks per CAEN tick plus an offset of 12.302 s with a few-us wander, AstroPix
  on-chip phase 189 (`docs/root_files.md`, `common/bic_streams.py`).
- **Not read out:** the pair-spectrometer hodoscope (so no per-event beam energy), the beam current
  and the SiPM bias and temperature (only the run-log values, HV 41.70 V, 343 uA total).
