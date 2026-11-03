# Synchronization: method and results


Every result is checked against a null control: the streams shifted by +10 ms, or the pairing
shuffled. Do not accept any correlation peak without one; several early "results" in this project
were artifacts that a control would have caught.

### 4.1 KCU1 to KCU0 (done)

- **Offset:** ts(KCU1) = ts(KCU0) - 1143665946 ticks. The residual is 0 (93%) or +4 ticks (7%,
  KCU1 one bunch late), with no drift.
- **Recipe for a new run:**
  1. Coarse offset from the RTCLOCK epochs: L = rt - 6.25 ns x ts per board, then
     L(KCU1) - L(KCU0) to ~ms. Drop the first stale `ts~0` buffer flush in each stream.
  2. Histogram ts(KCU1) - ts(KCU0) in 4-tick bins within +-20 ms of the coarse value.
  3. Take the exact mode inside the peak bin.
- **What the offset is:** a counter-epoch gap, not a latency. Matched triggers reach the host within a
  few ms of each other.
- **Physics check:** KCU1 (SFIL5) fires without KCU0 (SFIL1-4) in 0.049% of matched pairs, against
  6.55% of decoys.
- **Do not** use DAQ-slide-11 timeframe-index pairing: on 3062 it is off by 5.1 s.
- **Status:** matches the analyst's 9/29 slides exactly (same offset, same 93/7 split).

### 4.2 CAEN to KCU0 (done)

- **Timestamp:** use `GROUPTRIGGERTIME`, which the stock decoder joins into one 37-bit counter that never
  wraps over the run. **Do not use `EVENTTIME`**: it is a 32-bit counter with 2^31 backward jumps, and
  it produced an earlier false "CAEN shares no triggers" conclusion.
- **Coarse offset:** take the 1st percentile of rt - t_CAEN (not the median), because CAEN data reach
  the host up to 2 s late. That gives +12.62 s.
- **Fine offset:** a blind +-3 s scan finds one peak at -0.32 s, 35-50x above background, drifting
  +5.38 ppm.
- **Model:** a linear clock ratio 1.365340679 plus a running median (201 pairs) of matched residuals,
  which removes ~3 us of slow wander.
- **Result (2026-10-21 re-run):** 57.1% of CAEN triggers (861,778) match within 100 ns, with residual
  MAD 6.4 ns; it was 49.7% with the pre-fix decoder (fewer complete KCU0 waveforms). The +10 ms control
  gives 0.0%.
- **Comparison with the CAEN sync deck** (`context/JLab_Test_Beam_CAEN_KCU_sync.pdf`): its gap-landmark
  method gives ratio 1.36534063, which agrees with ours to 0.04 ppm, but only a +-45 us box. Ours is
  ~3000x sharper because of `GROUPTRIGGERTIME`.
- **No pulser in 3062:** the 10 Hz pulser (TIC slide 33) is not visible in 3062, so it was not used.

### 4.3 AstroPix to KCU0 (done statistically; per-event purity limited)

- **fpga_ts** (40 MHz, wraps at 2^32 every 107 s; unwrap it) is the SPI **readout** stamp, not the
  hit time. Readout runs in 1425-tick (35.6 us) slots, saturated 2/3 of the time, in ~40 ms windows
  with 70-100 ms gaps. Its drift against KCU0 is -0.30 +- 0.17 ppm.
- **On-chip `ts`** (8 bit) is a **10 MHz counter** (100 ns, wraps every 25.6 us) locked to the KCU
  clock. That is not the documented 2 or 2.5 MHz. For the true partner trigger K, in KCU ticks:
  `ts == ((K + 4) // 16 - 189) mod 256`. The evidence is a 20.3 sigma spike (19.3 pre-fix) against 2.8 sigma with
  `ts` shuffled between clusters; the divisors 32-640 are flat.
- **The phase is the same on every chip (2026-10-21, `astro_phase_split.py`, log in
  `docs/astro_phase_split_3062.txt`).** Measured from scratch over 2.35M clusters and 242M
  (cluster, trigger) pairs in the two latency windows, every chip peaks at phase 189: chips 0-8 at
  20.5 / 5.3 / 17.3 / 58.5 / 4.5 / 40.9 / 46.0 / 5.3 / 26.7 sigma (the weak ones are off the beam).
  Also 189 for row and column records separately, for both latency populations, in every 200 s bin,
  for every ToT bin and for single and multi-hit clusters. So the nine chips share one timestamp
  counter and one reset, and the matching is not tied to any particular chip. The 19 ms and 88 ms
  latency populations are present on every chip in the same proportions.
- **`fpga_ts` is not a valid trigger-level alignment on its own.** It stamps the SPI read, 19 or
  88 +- 5 ms after the hit; at 2.8 kHz a +-5 ms window holds ~30 triggers. It is still needed to
  place the cluster in run time and to select which 25.6 us wraps of `ts` are in play. The on-chip
  `ts` is the only sub-ms handle, and its 8 bits are why purity stops at ~80%: ~0.3 random
  candidates per cluster against ~0.07 true. Hardware fixes for the next run: paddle OR into a spare
  A-STEP input (stamped with `fpga_ts`), or a periodic sync pulse into both the KCU trigger card and
  the A-STEP, or a common counter reset; and log the A-STEP timestamp-clock setting.
- **Readout latency** (trigger to stamp) takes one of two values, 19 +- 5 ms or 88 +- 5 ms, measured
  relative to offset O = (envelope coarse offset, -8.828 s) - 0.43 s. The -0.43 s comes from beam on/off
  edges, with CAEN defining beam-on. Hits read more than ~30 ms before the end of a readout window fall
  in the 88 ms group; hits read in the last ~13 ms fall in the 19 ms group. The latency is the same for
  all chips and stable over the run. **Nothing tested predicts the +-5 ms jitter:** window position,
  queue depth, gaps, ToT, pixel, chip, cluster size, and a FIFO model.
- **Clustering:** consecutive hits whose `ts` is within +-2 of the first hit and which are read within
  1 ms form one cluster. There are 2.35M clusters: 1.93M singles and 397k row+col pairs.
- **Event matching (common-trigger version, 2026-10-20):**
  - **Candidates:** rows of the common-trigger table (4.4) in the right phase and latency window, so a
    partner trigger counts even when KCU0 lost its waveform.
  - **Ranking:** a likelihood ratio built from:
    - paddle pattern (KCU0, else CAEN) x (chip, pixel class);
    - SFIL2-3 beam-row energy and hottest SFIL2 cell, with an extra "no KCU0 waveform" class;
    - phase residual;
    - latency bin x read position.
  - **Training and purity:** the LR is trained on half the clusters and applied to the other half. Purity
    is calibrated on the control phase bands (60, 120).

  | log10 LR | matches | purity | est. true | previous (KCU0 only, max-pedestal paddles) |
  | --- | --- | --- | --- | --- |
  | > 2.0 | 106.4k | 70% | 74.6k | 70.8k, 68%, 48.4k |
  | > 2.5 | 50.8k | 81% | 41.0k | 34.8k, 80%, 27.9k |
  | > 3.0 | 8.6k | 88% | 7.6k | 9.6k, 88%, 8.4k |

  After the 2026-10-21 re-run (corrected decoders, 15% more complete KCU0 waveforms) the per-match
  purity selection gives: purity >= 0.5: 111.7k matches (79.1k true); >= 0.7: 53.4k (43.5k true);
  >= 0.8: 25.7k (22.0k true); 246k matches saved with log10LR > 1 (106k true); 62.8k of them on
  triggers without a KCU0 waveform (KCU1/CAEN only).

  - **Partner reach:** the matched partner has no KCU0 waveform (KCU1/CAEN only) in 7.2k of the
    log10 LR > 2.5 matches.
  - **Source of the gain** (matcher test half, true partners per cluster at 80% purity):
    | Candidates and paddle variable | True partners per cluster |
    | --- | --- |
    | old (KCU0 only, max - pedestal paddles) | 0.0103 |
    | KCU0 only, dip paddles (`TRIG=kcu0`) | 0.0135 |
    | common triggers, dip paddles | 0.0163 |
  - **Output:** `astro_kcu0_matches.npz` now holds `trig_ts`/`trig_row` and the `i_kcu0`, `i_kcu1`,
    `i_caen` and `pad` of the matched trigger, instead of `kcu0_ts`.

  - **Caveat:** the selection uses the paddles, so **do not study paddle response on the matched sample**.
    Use the unbiased background-subtracted method in `astro_geo_common.py`: phase-matched candidates
    weighted +1, control-band candidates weighted -1/2.
  - **AstroPix sees the beam:** it is ~92% beam-driven (1622 Hz with beam, 126 Hz in the
    CAEN-silent seconds).

### 4.4 Common trigger table (`trigger_table.py` -> `out/triggers_3062.npz`)

The NIM paddle OR feeds the KCU trigger card and both CAENs (TIC signal-routing slide 6). AstroPix gets
only the common clock and stays self-triggered.

- **Table:** one row per trigger seen by any of KCU0 (complete waveforms), KCU1 (complete) or CAEN.
  - All three are placed on the KCU0 clock and grouped within 100 ns.
  - No stream ever falls twice into one group.
  - Each row carries `pad` (bits P3=1, P4=2, P5=4, P6=8), taken from KCU0 or else from CAEN.

| Streams | Triggers | Share | With a KCU0 fragment |
| --- | --- | --- | --- |
| KCU0 + KCU1 + CAEN | 214k (155k pre-fix) | 4.7% | |
| KCU0 + CAEN | 648k (594k) | 14.3% | |
| KCU0 + KCU1 | 627k (457k) | 13.8% | |
| KCU0 only | 1963k (1795k) | 43.3% | |
| KCU1 + CAEN | 117k (120k) | 2.6% | 11k |
| KCU1 only | 436k (431k) | 9.6% | 55k |
| CAEN only | 530k (639k) | 11.7% | 37k |

- **Rates (2026-10-21 re-run):** the union is 4.53M triggers over 1486 s (3050 Hz). By stream: KCU0
  complete 2322 Hz (2391 Hz with fragment-only triggers), KCU1 937 Hz, CAEN 1014 Hz. Pre-fix: 4.19M,
  2819 Hz, KCU0 2018 / 2202 Hz, KCU1 782 Hz. The FPGA trigger counter says ~4150 Hz were issued.
- **No common busy:** the streams share the trigger but not a busy/veto, so each keeps a different
  subset.
- **KCU1 vs CAEN timing:** on the 120k triggers without KCU0, KCU1 and CAEN agree to MAD 6.6 ns. This
  confirms the KCU1 offset independently of KCU0.
- **Paddle pattern:** 89.7% of triggers have one (71.6% from KCU0, 18.1% from CAEN only).
- **KCU0 fragments:** a trigger whose KCU0 waveform was split by the decoder leaves fragments at
  trigger + k x 164 ticks (k = 0..19, within +-4 ticks for 96%).
  - They are flagged (`kcu0_frag`); with ROC_pool stitching the Packets they are down from 273k to
    103k flags (the remaining ones are 1-2-sample readouts and datagram-loss cuts).
  - About 80% of unmatched CAEN triggers have no KCU0 waveform at all, which points to KCU0 deadtime.

### 4.5 Geometry cross-checks (unbiased; `sync/out/astro_geo_*.png`)

Candidates are common-trigger rows. The SFIL panels need a KCU0 waveform; the paddle panels use the
common pattern.

- **Paddles:** P(paddle hit) is ~0.96-1.0 over its own chips and ~0.05 elsewhere. That is P4 on chip 0,
  P3 on chips 2-3, P6 on chips 5-6 and P5 on chip 8.
- **Along the fibres:** the SFIL1-4 beam-row (N-S)/(N+S) falls by ~0.1 from chip 0 to chip 8.
- **Vertical, SFIL2:** cell 3 holds >90% of the SFIL2 energy below AstroPix y = 11 mm. Energy moves
  to cells 1+2 at y ~ 12 mm.
- **Vertical, SFIL5:** energy moves from array B to array D as AstroPix y rises; D/(D+B) goes from 0.4
  to 0.75.

### 4.6 Event tree (`sync/out/events_3062.root (copy in data/)`, tree `events`)

**This is the file to analyse.**

- **Rows:** one entry per common trigger (4,533,542 after the 2026-10-21 re-run; 4,189,768 pre-fix), in
  time order. A stream that missed the trigger has
  `has_<stream> = 0` and zeros in its branches.
- **Branches** (full description in the `event_tree.py` docstring):
  - `trig_ts`, `t_s`
  - `has_kcu0`, `has_kcu1`, `has_caen`, `kcu0_frag`
  - `pad`, `pad_kcu0`, `pad_caen`
  - `sfil_s[4][4]` and `sfil_n[4][4]` ([cell][layer], KCU0)
  - `pad_dip[4]`
  - `sfil5_s[72]` (KCU1 ch 144-215)
  - `caen[16]` (Baby BCAL South 0-11, paddles 12-15)
  - `n_astro` and `astro_x/y/chip/nhits/log10LR/purity/latency_ms` (best AstroPix match)
- **AstroPix match counts:** 188.5k triggers carry a match with log10LR > 1, about 80.6k of them true
  (178k / 76k pre-fix). Cut on `astro_purity` (25.7k matches at >= 0.8).
- **Paddle caveat:** the AstroPix selection uses the paddles. For paddle studies, use the
  background-subtracted method (4.5), not `astro_purity` cuts.
- **Reading it:** works in plain ROOT
  (`TFile f("events_3062.root"); events->Draw("sfil_s[2][1]","has_kcu0")`) and in uproot. 392 MB.

### 4.7 Diagnostics, waveforms, event displays, shower shape (2026-10-20; `reproduce.sh diag`)

**Sync diagnostics (`sync_diag.py`, `out/sync_diag.png`).** Every link against run time, in 10 s bins.
- **KCU0-KCU1 offset:** the residual is 0 ticks (92.66%) or -4 ticks (7.34%), with no other value
  over 612k triggers.
- **CAEN-KCU0:** residual MAD 6.4 ns, flat over the run.
- **CAEN-KCU1:** on the 120k triggers without KCU0, MAD 6.9 ns.
- **Paddle pattern:** KCU0 and CAEN agree on 96.1% of triggers (21.8% for decoys), flat over the run.
- **AstroPix:** the phase-lock excess per cluster and both latency bands (19 and 88 ms) are flat
  over the run.
- **One dropout:** a KCU0 outage around t ~ 550 s shows in every KCU0-based panel.

**Waveform-level background subtraction (`wf_acc.cc`, `wf_keys.py`, `wf_bkgsub.py`,
`out/wf_bkgsub.png`).**
- **Inputs:** `wf_acc` (stock decoders) sums pedestal-subtracted waveforms for weighted trigger lists.
- **Classes:** AstroPix chip region (0, 2-3, 5-6, 8) x (signal phase band, control bands).
- **True-partner waveform:** (S_sig - 0.5 S_ctl) / (W_sig - 0.5 W_ctl).
- **Paddles:** in KCU0 and CAEN alike, only the paddle over each chip region saturates (KCU0) or pulses
  (CAEN); the others stay flat.
- **Showers:** SFIL3 beam-row and SFIL5 B/D pulses rise clearly above the random-trigger mean.
- **True partners:** 13.8k / 42.0k / 45.0k / 14.0k (KCU0) per region.

**Event displays.**
- **Static:** `evdisp.py` writes `out/evdisp/ev<row>.png` for 8 triggers seen by all three streams with
  one AstroPix match, the two most confident per paddle region. Purity is low on chip 0 (0.14), so
  those two are mostly illustrative. Each display has AstroPix + paddles, S/N cell maps in
  AstroPix-y order, the longitudinal profile, and KCU0/KCU1/CAEN raw waveforms.
- **Interactive 3D:** `evdisp3d.py` builds `out/evdisp/bic_event_viewer.html` (three.js) from `evdisp3d_template.html`.
  - **Measured:** x and y from AstroPix, paddle x ranges, cell order and SFIL5 array order.
  - **Placeholders:** z positions, thicknesses, cell heights, fibre length and the BBCAL layout.
  - **For a real GEANT4 view:** use the setup metrology or simulation geometry (TIC slide 37,
    J. Richards).

**Shower shape and relative energy (`shower.py`, `out/shower.png`, `out/shower.log`).** No per-event
beam energy exists (the August paddles were never calibrated). The paddle is used as an ordered energy
tag: in May, P4 < P3 < P6 < P5 (0.48 / 0.55 / 0.57 / 0.68 GeV), rising from North to South.
- **Energy proxies:**
  - SFIL1-4: per-cell sqrt(S x N), cells 1, 3 and 4. This cancels attenuation along the fibre.
  - SFIL5 South sum and Baby BCAL South sum, both single-ended.
- **Low-energy triggers are geometric:**
  - 39-46% of triggers have SFIL1-4 sum <= 1200 ADC, the low peak.
  - Among AstroPix true partners the low share is 3-7% at y = 3-9 mm, rising to 50-80% at y > 10 mm
    and at the bottom edge.
  - **Raw waveforms** (`k1empty_wf.py`, `out/k1empty_wf.png`) show no pulse in SFIL1-4, SFIL5 or the Baby
    BCAL. The KCU1 pedestal is normal. Only the paddle fired.
  - **Correction:** these triggers do *not* show in the BBCAL. The BBCAL peak-height sum has a noise
    floor of ~900 ADC, which earlier read as signal.
  - **Interpretation:** paddle-only triggers, from particles crossing a paddle outside the SFIL/BBCAL
    acceptance (the AstroPix low share rises at the y edges).
  - **Selection:** analyses use showering triggers (sum > 1200).
  - **The KCU0 sum is never 0:** max - pedestal of noise is 10-20 ADC per channel, plus ~300 ADC from
    the noisy SFIL1 North column, giving ~730 ADC with no signal. That is the KCU1 = 0 band of the 9/29
    ADC-sum plot.
  - **The KCU0 >= 1500 tail (3.3% of matched triggers):**
    - It is a slow baseline ramp of about +40 ADC over the 500 ns window, in every channel, not a pulse.
    - It is not caused by a recorded trigger: time since the previous trigger is the same as for showers,
      and none is within 20 us.
    - Use a pulse-shape amplitude, not max - pedestal, for such events.
- **Longitudinal (showering triggers):**
  - SFIL4/SFIL2 = 0.638 / 0.675 / 0.695 / 0.721 (P4 -> P5): the shower maximum moves deeper with
    energy.
  - SFIL2, SFIL4, SFIL5 and the BBCAL medians all rise in paddle order.
  - SFIL3 is saturated in ~63% of showers.
- **Peak width sigma/mu** (Gaussian core, uncalibrated, an upper bound on resolution: includes beam
  spread, leakage, saturation and the single-ended BBCAL). These are superseded by the
  pedestal-subtracted values in 4.8.
- **Lateral:** an erf fit to the cell-boundary step against AstroPix y gives sigma = 0.92, 0.95, 1.65
  and 2.40 mm for SFIL1-4. The shower broadens with depth. The boundary sits at y = 9.3-12.1 mm,
  depending on the layer.
- **TOT can't fix saturation:** it is enabled on only some channels, reads at most 63, and is absent on
  others. `ts_dump kcutot` dumps it for inspection.
- **Absolute energy resolution needs** an energy calibration (PS paddles or a run with known
  energies) and a saturation-robust amplitude (template fit of the waveform tail).

### 4.9 TOA / TOT timing between the North and South fibre ends (2026-10-21, `timing/toa_timing.py`)

**Question:** can the H2GCROC3 TDC (TOA 25 ps LSB, TOT) give ps timing between the two ends of an
SFIL fibre, both on KCU0?  **Answer for run 3062: no, and the reason is configuration, not the chip.**
Inputs: `ts_dump kcutime` (KCU0, 36 SFIL/paddle channels) and `ts_dump kcu1time` (KCU1 ch 144-179).
Figure `out/toa_timing.png`, log `out/toa_timing.log`.

- **The TOA response of the KCU0 SFIL channels is not uniform (observed; the cause is open).** Of
  the 32 SFIL channels, 21 carry TOA words on 15-47% of their pedestal samples, in runs that last many
  samples rather than as independent 25 ns crossings (the per-waveform count of pedestal samples with a
  word is flat or bimodal, not binomial), 10 carry no TOA word at all, not even on pulses that saturate
  the ADC (e.g. S3c3 ch52: 0 of 37k pulses above 800 ADC, and TOT on 0.7%), and one, N3c3 ch196, is
  clean (1% on pedestal samples, 85% at the pulse). Neither the ADC pedestal nor the pedestal noise of
  a channel predicts which group it is in (rank correlations -0.3 / +0.3), so this is not simply a
  common threshold meeting different baselines; a per-channel TDC mask or trim is the kind of thing
  that would do it. The paddles on ASIC1 do fire, at sample 9 (one after the ADC dip) in 91-97% of
  hits. The South partner of ch196, S3c3 ch52, has no TOA, so there is no North-South pair with a
  working TOA on both ends. On KCU1 (single SiPMs) all 16 live channels are clean (0.2% on pedestal
  samples, 99% at the pulse).
- **The TOA word is stored one sample after the sample in which the ADC crosses the pulse edge**
  (84% at +1, 16% at +0, on KCU1 with a 20% crossing), not in the same sample; the earlier verdict
  "TOA is noise" looked at the wrong sample and treated the 10-bit code as a linear number.
- **The code's fields do carry time, but not with the datasheet mapping.** Bit by bit against the ADC
  crossing time on KCU1: bits 9 and 8 carry +10.3 / +6.9 ns (the 160 MHz coarse counter works), bit 7
  +2.0 ns, bit 6 +0.9 ns, bits 5:0 < 0.2 ns. Fitting `t = 25 x sample + a(coarse) + b(coarse) x mid`
  per channel gives slopes of -0.2, 0.0, +0.26 and +0.45 ns per mid code for coarse = 0, 1, 2, 3 and
  intercept differences that are not multiples of 6.25 ns (the same on every KCU1 channel). Part of
  that is the reference: a linear interpolation of a 2-sample rising edge is itself non-linear at the
  ns level. The TDC therefore needs a real calibration (a delay scan with an LED or laser, as in the
  `PlotEventBICSingle.C` approach, or the datasheet's DLL calibration), not an ADC-based one.
- **With the ADC-referenced per-channel calibration, TOA vs TOA between two KCU1 channels in the
  same waveform (A 250-900, timewalk-corrected with 1/sqrt(A)):** MAD 0.84-1.10 ns per pair, i.e.
  600-780 ps per channel. That is an upper limit set by the reference: on these fast pulses (rise < 1
  sample) the ADC interpolation is pinned to the sampling phase (panel d: for two of the four coarse
  values it does not move at all), which distorts the fitted tables; the small ADC-vs-ADC pair
  difference (0.36-0.47 ns) is that pinning cancelling between channels, not a resolution. On KCU0 the
  summing-board pulses rise over two samples, the ADC crossing is a valid clock there and gives
  1.6-2.0 ns per end for North-South (section 4.3 and panel f).
- **The TDC is a real per-channel time measurement.** Two KCU1 channels of the same waveform land
  in the same 6.25 ns quarter (sample + coarse bits) in 85% of cases and in adjacent quarters in 13%;
  within the same quarter their mid codes agree to MAD 4.4 codes against 13 for a shuffled pairing;
  only 1% of the codes are identical, so this is not a shared per-event value. The coarse bits are
  uniform (0.23 / 0.26 / 0.26 / 0.24), as a free-running 160 MHz counter against random trigger phase
  must be.
- **Self-calibration without the ADC (`timing/toa_selfcal.py`): ~190 ps per channel.** Per channel a
  table over (coarse, mid bin of 4 codes) is fitted so that all 16 KCU1 channels of a waveform agree,
  on the even waveforms; on the odd waveforms the pair difference has MAD 0.22-0.35 ns, weighted mean
  0.27 ns over 94 pairs (58k test pairs), i.e. **192 ps per channel** after timewalk correction. That
  includes the SiPMs' own time jitter, so it is an upper limit on the TDC; the fine bits add nothing
  measurable at this level. (A coarse-only version of the same fit is trivially 0 because the pair
  shares the quarter; the mid bits are what carry the sub-quarter time.)
- **TOT** fires only on saturated pulses (Tc set): 23-51% of beam-row hits on the channels whose TOT
  threshold is sane, 0.4% on S3c3 whose thresholds are too high. It is a usable amplitude where the
  ADC clips (10-bit code, decompress 8(v-512)+4 for v >= 512), hence a timewalk handle, and nothing
  more for timing.
- **For the next run:** (1) a per-channel TOA (and TOT) threshold scan, set ~5 sigma above the
  pedestal noise; (2) the TDC DLL calibration of datasheet 1.2.2 (`CALIB_CHANNEL_DLL`, the CTDC/FTDC
  bias DACs) and an LED/laser delay scan per channel to map the 10-bit code to ns; (3) check the
  TOA-to-ADC sample alignment (`GLOBAL_LATENCY_TIME`); (4) keep the ADC waveform as the cross-check.
  With the thresholds set, the self-calibration of `toa_selfcal.py` already gives ~190 ps per channel
  on single SiPMs, so North-South timing at the 200 ps level per end is realistic for the next run even
  without a delay scan; a delay scan would make the tables absolute and show how much of the 190 ps
  is the TDC. Without the thresholds the ADC crossing is the only clock (1.6-2.0 ns per end on the
  summing boards).

### 4.10 What the KCU0 waveforms look like (2026-10-21; `decoding/kcu_hists.cc`, `kcu_hists_plot.py`, `wf_overlay.py`)

Figures in `sync/out/kcuhists/` (KCU0 segment 0001): `SFIL_South_*.png`, `SFIL_North_*.png`,
`paddles_*.png` are the panels of the LED macro `PlotEventBICSingle.C`, one panel per channel
labelled "ch n  detector", for `ADC` vs sample, `TOA` code vs sample, decompressed `TOT` vs sample,
`ADCcorr` (ADC vs time aligned by the FIRST TOA word of the waveform, the macro recipe) and `ADCpulse`
(aligned by the pulse's TOA word, one sample after the largest ADC step); `ASIC_<q>.png` are channel x
value maps per ASIC; `wf_overlay_pulses.png` / `wf_overlay_random.png` show 30 individual waveforms per
channel. The channel diagram along the beam line is `docs/kcu0_channel_map.png` (`channel_map_fig.py`).
First 40k complete waveforms for the histograms, the first 100k of the `kcutime` file for the overlays.

**Reconstruction check.** `decoder_tests/kcu_decoder_check` now compares every channel word of every
sample with the raw bytes of the line: ADC (bits 29:20) and TOA (bits 9:0) equal the raw bits for
75,549,924 of 75,549,924 words on KCU0 and 75,530,916 of 75,530,916 on KCU1 (1500 packets each), on
top of the TOT / CRC / Hamming / timestamp-order checks (`decoder_tests/README.md`). Individual
waveforms are smooth: flat pedestal in samples 0-4, rise at 5-6, peak at 7-8, no sample scrambling.

**What the data contain** (first 100k waveforms of segment 0001):

- Pulse shape, pedestal-subtracted mean of 300-600 ADC pulses, ch 52 (S3 cell 3):
  `0 0 0 0 0 9 204 452 526 513 477 438 402 370 343 322 305 293 284 279` (samples 0-19). Rise over two
  samples, peak at sample 8, then a slow tail still at 55% of the peak 275 ns later. Smaller pulses
  (ch 10, ch 29) fall to 10-30% by sample 12 and rise again towards sample 19 (slow component).
  Paddles: negative, one or two samples wide, rail at ADC 0 at sample 9, with a positive blip at
  sample 8 before the drop; many waveforms rail again at samples 15-19.
- Saturation at 1023: ch 52 36%, ch 196 34%, ch 0 22%, ch 29 20%, ch 173 12%, ch 144 3% of ALL
  waveforms (60% of beam-row waveforms carry a pulse > 100 ADC). TOT words appear only there
  (0.002-0.30 per waveform) and nowhere else.
- Noise within a waveform (samples 0-3 about their mean): 4-6 ADC on most SFIL channels, 2.5-3 on the
  paddles, 43-52 on the North SFIL1 board (ch 146, 150, 158; 98% of waveforms spread > 50 ADC in
  samples 0-3) and 20 on ch 154: that board is noisy, as the shift notes say.
- Per-event baseline shift: the event-to-event pedestal rms is 17-35 ADC on quiet SFIL channels (5-95%
  range 60 ADC wide) with NO slow drift (same rms and median in every 40k block over 90 s). It is
  shared: rho = 0.82-0.99 between quiet ASIC0 channels, 0.95 within ASIC2, 0.77 between ASIC0 and
  ASIC2, 0.01 with the paddles on ASIC1; subtracting the per-event mean of the other quiet channels
  leaves 7-11 ADC, and it does not depend on the event's own shower amplitude. So it is a common shift
  of the SFIL readout (bias / summing-board side), not the chip and not a trigger-phase effect. A
  per-event pedestal from samples 0-3 of the same channel removes it (4.8 does that); the HGCROC
  common-mode words are not exposed by the decoder (a candidate accessor, not added).
- Pile-up: 2.0-2.5% of beam-row waveforms start on the falling tail of an earlier pulse (pedestal > 150
  and falling through samples 0-4). Stored triggers are > 100 us apart, so these tails come from
  untriggered beam particles in the preceding few hundred ns, consistent with the long pulse tail.
- Why the 2D density plots look smeared: 40-60% of the beam-row waveforms carry pulses spanning a
  10x amplitude range with long tails and 20-36% saturation of all beam-row waveforms (53-66% of shower waveforms, rising P4 < P3 < P6 < P5), drawn at log z over 40k triggers, plus
  the 2% pile-up and the baseline shift. The macro's first-TOA-word alignment gives a 25 ns comb on
  KCU0 because the first TOA word is usually a pedestal-sample word (4.9); aligning by the pulse's
  own TOA word (`*_ADCpulse.png`) still leaves a comb on ch 196 and only partial structure on the
  paddles, i.e. the nominal fine x 25 ps + mid x 0.2 ns + coarse x 6.25 ns formula is not linear in
  time here (4.9, self-calibration).

Reproduction:
```
cd bic_testbeam/sync && make kcu_hists
./kcu_hists out/kcuhists/kcu0_hists.root ../data/raw/beam_ROC_0-00003062-0001.evt 40000
python3 kcu_hists_plot.py && python3 wf_overlay.py && python3 channel_map_fig.py
```

### 4.11 Anatomy of one KCU0 trigger (2026-10-21; `decoder_tests/kcu_event_anatomy.cc`, `results/anatomy_kcu0_seg{0,1}.txt`, `anatomy_kcu1_seg0.txt`)

One KCU0 trigger = 20 samples x 8 links = 160 lines of 192 bytes (30,720 bytes). Every line carries
the FPGA wrapper (64-bit timestamp, TRIGGER_IN, TRIGGER_OUT, EVENT_COUNTER) and one H2GCROC3 frame
(header with BX / EC / orbit / Hamming flags, two CM words, 36 channel words with ADC / TOA / TOT /
Tc / Tp, the calibration channel, CRC-32). So per trigger: 288 channels x 20 samples of ADC, TOA and
TOT, plus 16 CM and 8 calibration waveforms. All 8 links are assembled on the same FPGA timestamp,
so the SFIL North and South ends, the paddles and ASIC3 are synchronous by construction; the only
loss is a missing line (UDP), which blanks 36 channels of that sample (1.5% of samples in complete
waveforms miss a link).

What the 288 channels contain (first 30k complete waveforms of segments 0000 and 0001, same numbers):
- 32 SFIL channels (16 South on ASIC0, 16 North on ASIC2) and 4 paddles on ASIC1 carry signal.
  81% of triggers have an SFIL pulse > 100 ADC, 65% a beam-row (cell 3) pulse.
- ASIC3 (Baby BCAL North, ch 216-287): ch 217-251 sit at a clean pedestal (21-81 ADC) with a
  pulse > 20 ADC in at most 0.1% of triggers; ch 216 and 252-287 are railed at 1023. No Baby BCAL
  North signal in this run; its South end is in the CAEN.
- The other 180 channels never show a pulse > 100 ADC (0 of 30k triggers).
- TOA words: ~230 per trigger over all channels, mostly the pedestal-sample words of 4.9; TOT words:
  0.8 per trigger (saturated beam-row pulses only).

Timing information per trigger:
- The FPGA timestamp of sample 0 (6.25 ns ticks) is the trigger time; samples 1-19 follow at 164
  ticks (1.025 us) each, i.e. the timestamp is the readout cadence, while the ADC samples are 25 ns
  apart (BX step 1). ts19 - ts0 = 3116 ticks in every complete waveform.
- TRIGGER_OUT is the board's trigger counter: +1 per trigger, constant within a waveform, never
  reset at run start (KCU0 starts at 13,047,448, KCU1 at 222,293), 4.10 kHz on both boards. Its steps
  between consecutive complete waveforms count the triggers without a complete waveform exactly.
  Timestamp-matched KCU0/KCU1 triggers have TRIGGER_OUT(KCU1) - TRIGGER_OUT(KCU0) = -12,854,408
  for the first 12 s of the run, then -12,854,406 at 22 s: the difference slips across gaps, so it is
  a check, not a key (as noted on 2026-09-15); the timestamp residual of matched triggers is 0 or -4
  ticks.
- EVENT_COUNTER is a sample counter (+1 per sample, 20 per trigger), kept per link; the decoder
  stores the last line's value per sample, so links mix and it is not usable for trigger counting.
- TRIGGER_IN is always 0 (not wired).
- BX / EC / orbit of the ROC header are not exposed by the decoder.

Are all KCU0 triggers paddle triggers? Yes: 99.85% of complete waveforms have a paddle pulse (drop
> 40 ADC below pedestal; the 0.15% without also fail a 20 ADC cut), multiplicity 1: 87.8%, 2: 9.5%,
3: 1.5%, 4: 1.1%; per paddle P3 27.8%, P4 25.4%, P5 35.1%, P6 27.3%. Only 2 of 19.5k beam-row
pulses come without a paddle pulse.

Missing triggers: segment 0000 opens with 4 triggers, then 7.14 s with 29,228 triggers issued and
nothing recorded, then KCU1 starts (its counter epoch is 7.147 s after KCU0's, see
`bic-kcu-clock-epoch`); after that 136 of 30k triggers (0.45%) lack a complete waveform in the first
30k. Segment 0001 (UDP-loss stretches) lacks 44% (23,443 of 53,443 issued). Over the run, 3.45M
complete KCU0 waveforms against 4.10 kHz x 1479 s = 6.06M issued, i.e. 57%.

**The LED macro's "magic number" (2026-10-22).** `PlotEventBICSingle.C` aligns waveforms with
`t = 25 (s + (TOA_ns < 6.37 ? 0 : 1)) - TOA_ns` for the TOA word in sample s, the 6.37 ns being the
boundary of the TOA window relative to the ADC sampling in that (LED) configuration. Tested on the
115k KCU1 pulses of run 3062 (TOA words within one sample of the 20% ADC crossing):
- the words split by position in the window, but the split is at **20.7 +- 0.5 ns**, not 6.37: words with
  TOA_ns > 21 sit in the crossing sample, words with TOA_ns < 20 one sample later (the 16% / 84% of
  the earlier note);
- the convention that makes the arrival time continuous against the ADC crossing is the plain
  forward one, **t = 25 s + TOA_ns with no window rule** (MAD 4.3 ns, 2.5% of pulses off by a sample),
  not the macro's `-TOA_ns` form (MAD 8-11 ns, 30-50% off by a sample, with either boundary). The
  TOA counts up from the window start in this configuration, and the window edge falls 4 ns before
  the next ADC sample. The remaining residual slope against TOA_ns (about 10 ns across the window) is
  the per-channel non-linearity above (use the self-calibration tables).
- `decoding/kcu_hists` now uses the forward convention by default (`toa_mode=macro` and a wrap value
  reproduce the macro), as do `toa_time()` in `common/bic_streams.py` and `ToaTime()` in
  `common/bic_streams.C`. The 6.37 (or 20.7) ns boundary is configuration dependent (GLOBAL
  latency / phase settings): measure it per run from the word split.

### 4.12 Pulse fits and noise (2026-10-22/23; `waveforms/kcu0_sfil_plots.C`, `waveforms/kcu0_sfil_fit.C`)

- `kcu0_sfil_plots.C`: raw, pedestal-subtracted (samples 1-3), TOA and TOT versus sample for every SFIL cell and
  the paddles in one figure per quantity laid out like the setup, with trigger-level flags (complete, CRC, Hamming,
  order, overflow, paddle, single paddle), waveform-level cleaning (amplitude, pedestal spread, saturation,
  pulse-window TOA), error counters per link and a cut flow.
- `kcu0_sfil_fit.C`: per-channel pulse templates (mean normalised pulse aligned on the 20 % crossing, clean
  unsaturated pulses 100-1000 ADC) fitted with amplitude and time by a chi2 scan; saturated samples masked;
  selection amplitude > 60 ADC, peak sample 7-9, chi2/ndf < 8 with 5 ADC noise and 5 % shape error; SFIL1 North
  column masked. Full run (condor): 2.91 M triggers, 6.98 M fitted amplitudes, file
  `/sphenix/tg/tg01/jets/tmengel/eic/kcu0/fit_full/`. The fit extends the amplitude range to ~1800 ADC; above
  that the fit pins (pile-up at the edge of the spectra), so the TOT is needed for the largest pulses.
- Noise on triggers without a beam-row pulse: 3-5 ADC per sample on clean channels, max(samples 4-12) -
  pedestal below ~40 ADC at the 99th percentile on the beam row; SFIL1 North column 19-50 ADC per sample.
  The trigger-to-trigger pedestal rms of 20-45 ADC (60-75 on the beam row) is the common baseline shift plus
  pile-up, not random noise.

### 4.13 Photostatistics from the two fibre ends: data-only channel gains (2026-10-23; `sim/analysis/photostat_ns.py`)

Purpose: measure the ADC per photoelectron of the SFIL channels from run 3062 alone, to test the equal-gain assumption
behind the layer-share plots without anchoring to the simulation.

Method. For a cell read at both ends, South and North see the same scintillation light L (attenuation is fixed for a
given paddle). After equalising North to South by the median amplitude ratio r, x = S and y = N/r are two measurements
of L with independent fluctuations, so in a bin of s = (x + y)/2:
sigma_S^2 = var(x) - cov(x, y) and sigma_N^2 = r^2 [var(y) - cov(x, y)]. Photostatistics gives sigma^2 = c + a A with
a = g F (ADC per photoelectron times the SiPM excess-noise factor). Pairs: template-fit amplitudes (`kcu0_sfil_fit.C`)
with both ends accepted and unsaturated; outliers beyond 5 MAD in x - y (pile-up) dropped per bin; linear fit for
100 < A < 450 ADC (the unsaturated requirement truncates the variance above ~600 ADC; a quadratic term from the 5 %
template shape error and the position spread is kept out of the range). Output `photostat_ns.txt`, `photostat_ns.pdf`.

Results (a in ADC per photoelectron x F, statistical errors):

| pair | pairs | N/S amplitude | a_S | a_N | a_N / a_S | photoelectrons N/S |
| --- | --- | --- | --- | --- | --- | --- |
| SFIL2 cell 3 | 321,836 | 0.849 | 14.09 +- 0.16 | 11.68 +- 0.14 | 0.83 | 1.02 |
| SFIL3 cell 3 | 519,678 | 1.033 | 9.47 +- 0.06 | 9.81 +- 0.07 | 1.04 | 1.00 |
| SFIL4 cell 3 | 437,680 | 0.352 | (55.8) | (35.0) | - | - |
| SFIL2 cells 1, 2 | 6.6k, 4.2k | 0.43, 0.53 | 13.7, 10.8 | 7.2, 5.6 | 0.52 | 0.83, 1.02 |
| SFIL3 cells 1, 2 | 38k, 57k | 1.15, 0.82 | 7.8, 11.1 | 9.6, 9.2 | 1.24, 0.83 | 0.93, 0.99 |

- The two ends of SFIL2 and SFIL3 see the same number of photoelectrons (last column 1.00-1.02): the SFIL2 North/South
  amplitude ratio 0.85 is a gain difference, not light.
- The South beam-row gains differ between layers: SFIL2 / SFIL3 = 1.49 +- 0.02 (stat). Equal gains, assumed in the
  share plots, is wrong by 50 % between these two layers. Channel-to-channel within a layer the spread is +-20 %.
- SFIL4 does not give photostatistics: sigma^2 grows quadratically from the lowest amplitudes. Its North end is
  anomalous: cells 1 and 2 (channels 187, 186) carry a constant ~130 ADC "amplitude" independent of the shower (SFIL3's
  neighbours scale with the shower as expected), and North cell 3 is 0.35 of South cell 3 while the layer sums give
  0.64 only because of those two channels. SFIL1 has no North end (masked). So the South gains of SFIL1 and SFIL4 stay
  unmeasured.
- Against the simulation-implied gains (layer sums, `sfil_gain_calibration_v{1,2}.txt`): the v2 geometry implies
  SFIL2 / SFIL3 = 1.29-1.41 (4-6 GeV), the pure-Pb v1 1.09-1.20. The measured 1.49 supports the v2 shower profile.

Assumptions and caveats: the same excess-noise factor F at both ends and in both layers (F rises with overvoltage, so a
breakdown-voltage spread of the SiPMs would move a by more than g alone; a 10 % allowance is carried); the pairs are
unsaturated showers only (lower energies); the number is per fitted-amplitude ADC, so it equalises gains but gives no
MeV scale; pile-up removed by the MAD cut only. `stack_shares_pe.py` applies the SFIL2/SFIL3 ratio to the layer shares
(`stack_shares_sfil_pe_v{1,2,3}.pdf`) with SFIL1 and SFIL4 assumed at the SFIL3 gain; its `nosat` / `nosat_raw` modes keep only
triggers without any saturated channel (30 % of the showers) and apply the equivalent threshold to the model: raw peaks and fits
then agree within 0.004, the shares move to 0.037 / 0.200 / 0.391 / 0.373 (model with the same cut and without cell 4, which is
empty in the data: 0.054 / 0.205 / 0.354 / 0.387 at 4 GeV).


**Single-ended readouts (SFIL5 South, Baby BCAL): no photoelectron calibration from run 3062 (tried 2026-10-23).**
(1) Intra-array pairs on SFIL5 South: the 16 single SiPMs of an array are coupled to the fibres directly (no light
guide), so two SiPMs do not see the same light: channel-channel correlation 0.16-0.21, and var(x_i) - cov(x_i, rest)
grows as ~400 x A, thirty times photostatistics (the fibre-level shower fluctuation). (2) Dark-count single-pe spectra
in the pre-pulse samples 1-5 (`pe_kcu1_prepulse.npz`): no peak above the noise tail on any KCU1 channel; local
maxima above 3 ADC occur in 0.1 % of triggers and their 4-sample integral spectrum is the noise tail (sigma 1.7 ADC per
sample); a quiet-trigger all-sample spectrum shows a bump at ~16 ADC x samples, but it comes from shower light that
survives the quiet cut, not from dark counts. The KCU0 summed channels (16 SiPMs, noise 3.5-4 ADC) show no structure
either. (3) CAEN channels: the local-extremum spectra show irregular bumps spaced 8-16 counts in both the pre-pulse
and tail windows, consistent with DRS4 sampling structure, not with a photoelectron comb. Conclusion: SFIL5 and the
Baby BCAL stay simulation-anchored until the LED runs (3179-3190) or cosmic runs (3088-3172) are analysed.
`stack_shares_calib6.py ... pe` draws the six-element figure with the SFIL2/SFIL3 photostatistics ratio applied to
SFIL1-4 and those two points anchored (`stack_shares_calib6_pe_v2.pdf`).

**SFIL5 North's four live SiPMs as a gain reference (2026-10-24, `sim/analysis/north_photostat.py`).** KCU1 ch 54, 60, 65, 66
(the D-array block of ASIC0) are correlated pairwise at 0.74-0.83 (two bare South SiPMs: 0.10), as SiPMs behind one light guide
should be, and anti-correlated with the South B sum: they are D-cell detectors seeing mixed cell light. Photostatistics between
them gives the KCU1 single-SiPM channel gain: a = 32 / 34 / 24 / 24 ADC per photoelectron (x excess noise) for ch 54 / 60 / 65 / 66,
consistent to 2 between partners; the sigma^2 curve bends above ~450 ADC (imperfect mixing and the approach to the ceiling), so take
+-15 %. Against the KCU0 summed channels (SFIL2 14.1, SFIL3 9.5) the KCU1 single-SiPM chain delivers 2.5-3.5 times more ADC per
photoelectron. The pairing with the South D sum does not work (correlation 0.2-0.3): the bare South SiPMs cover about a third of the
fibre face, so their sum is a noisy estimate of the cell light and the method attributes that sampling variance to the South. What
the North gives is therefore the KCU1 ADC/pe scale, not a calibration of the SFIL5 South sum; with it, the raw SFIL5/SFIL4 ratio of
~28 at low amplitude splits into a chain-gain ratio of ~3 and a photoelectron ratio of ~10 per MeV between bare SiPMs and the
light-guide-plus-summing chain, a number that the LED runs should confirm before it is used.

### 4.8 Pedestal-subtracted amplitudes and plotting style (2026-10-20)

- **Why max - pedestal is biased:** on a noise-only waveform the maximum of 20 samples sits above the
  pedestal. That gives +10-20 ADC per channel and ~700 ADC per event in the KCU0 SFIL1-4 sum. It produced
  the "KCU1 empty, KCU0 non-zero" band of the 9/29 ADC-sum check.
- **Fix (`ts_dump kcuadc/kcu1adc`, branches `win`, `cm`, `w1`, `cm1`):**
  - Amplitude = mean of the peak samples (KCU0 6-9, KCU1 5-8) - pedestal (mean of samples 1-3) - common
    mode. The common mode is the same quantity for the median unconnected channel of that ASIC.
  - Noise now averages to zero. The common mode also absorbs most of the event-wide baseline ramp.
  - CAEN: `qint` = sum over samples 20-159 of (pedestal - ADC), the pulse integral.
- **Event tree:** `sfil_s_w`, `sfil_n_w`, `sfil5_s_w`, `caen_q`. The old `sfil_s`, `sfil_n`, `sfil5_s` and
  `caen` are kept; use the old ones only for saturation flags (max - ped > 900).
- **`sumcheck.py` (`out/sumcheck.png`):**
  - The circled band moves to the origin: KCU0 sum of the same KCU1-empty triggers, median 703 -> -52 ADC.
  - Matched r = 0.73 (0.84 with max - ped), decoys 0.00.
- **`shower.py` now uses these amplitudes.**
  - Results: SFIL4/SFIL2 = 0.54 -> 0.60 (P4 -> P5). sigma/mu (SFIL1-4 + BBCAL) = 0.231 -> 0.184.
    Lateral width 1.05 / 1.36 / 1.73 mm for SFIL2 / 3 / 4.
  - SFIL1 has almost no signal and is dropped from the lateral fit.
  - The earlier, smaller sigma/mu values (0.164 -> 0.131) were biased low: the BBCAL peak-height noise
    floor (~900 ADC) inflated the mean.
- **Plot style:** every figure uses `plotstyle.py`:
  - the sPHENIX style (`sPhenixStyle.C` from the sPHENIX rootmacros), with no experiment label;
  - panel descriptions in the top margin, auto-shrunk to fit;
  - legends with no box, given headroom so they never cover data;
  - no axis exponents.

### 4.14 Clipped KCU1 pulses: TOT as the measure (2026-10-24, `waveforms/kcu1_tot/`)

- Extraction: `extract_kcu1_fit.py` turns the full-run fit tree into compact arrays (`kcu1/fit_full/kcu1_fit_arrays.npz`).
- Shape test on raw waveforms (300k triggers): the single-SiPM pulse narrows with amplitude (y[pk+3]/y[pk] 0.75 -> 0.50 from
  100 to 940 ADC) and, once clipped, the sample after the clipped stretch is 720-770 ADC independent of the number of clipped
  samples (1-7). The KCU0 summed channels show neither effect below 800 ADC.
- Decoded TOT versus clipped samples: median 502 (nsat 1), 668, 1148, 1588, 2116, 2396, 2820, 3300 (nsat 2-8): 166 units per
  25 ns; ceiling crossing T_c = 419. TOT present on 63 % (nsat 1) and 89-94 % (nsat >= 2) of clipped pulses; for the rest the
  median TOT of their nsat is used. Unsaturated pulses with a TOT word (A_fit 650-800, TOT 350-650) give the toe slope 0.33.
- Recovery: A = 943 + b (TOT - T_c). Three attempts to fix b from the data failed for identifiable reasons (neighbour predictor:
  position-driven, factor 2.5 spread; SFIL5/SFIL4 linearity: the model's own <L5|L4> falls as steeply; SFIL5 width: matches the
  model at b = 0 unless the paddle energy spread is 20-30 %). b = 0.40 adopted, band 0-0.76 (`kcu1_tot_shape.{txt,png}`).
- Output `sim/analysis/sfil5_recovered.npz` (A5, A5lo, A5hi, A5old per events row); `sfil5_sfil4_gain.py rec` and
  `stack_shares_sfil5.py ... rec` use it.

### 4.15 AstroPix x as a position tag: attenuation, gain ratios, energy curve (2026-10-24, `calibration/attenuation/astro_x_attenuation.py`)

- Join: fit cache -> KCU0 header ts0 -> events row -> `astro_x`, `astro_purity`; purity >= 0.3 keeps chips 3, 5, 6, 8 (45.6k
  showers). Stack frame x = astro_x - 88.75 mm (chip 4 centre at 0, North negative; P4 chip 0 = -80, P5 chip 8 = +80).
- Per layer, beam cell, both ends unsaturated and > 50 ADC: median of ln(A_N/A_S) in 10 mm bins, straight-line fit: slope
  = -2/lambda_eff, intercept = ln(g_N/g_S). Geometric mean sqrt(A_N A_S) per bin, normalised at x = 0, is E(x)/E(0) up to the
  non-exponential part of the attenuation; the layer-sum version (fits) is used for the paddle energy ratios.
- Limits: 10 mm bins, the matched-chip coverage (-32 to +88 mm), the unsaturated selection (biases the beam-cell curve low at
  P5), SFIL4 North cell 3 on a v2 board.

### 4.16 North - South timing from the ADC crossing and the KCU0 TOA (2026-10-24, `timing/toa_position/`)

- `toa_position.py`: valid TOA = word on samples ipk-2..ipk of a pulse > 200 ADC; occupancy self-calibration (as on KCU1)
  fails on KCU0 because the code occupancy is not uniform. `toa_adc_calib.py`: per channel f(code) = median(t_cross - 25 s_TOA)
  in code bins, t_cross the linear 50 % crossing of the ADC samples; even so the TOA residual is 4-8 ns rms and the N - S TOA
  difference is double-peaked (two code populations). The ADC crossing alone gives v_eff = 157 / 159 mm/ns (SFIL2 / SFIL4
  cell 3) from the four paddle positions, rms 1.8-2.6 ns per layer.

### 4.17 Per-event fibre position: timing from the template t0, amplitude asymmetry, combination (2026-10-24, `event_x_walk.py`)

- Timing per layer: dt = t0_N - t0_S from the KCU0 template fits (walk-free by construction: constant template, 20 % crossing,
  clipped samples masked). Two-parameter fit dt = c0 + c1 x_paddle on the paddle-tagged showers gives c1 = 2/v_eff and the
  fixed N - S delay; a third term in ln(A_N/A_S) is degenerate with x (attenuation) and is reported but not applied.
- Amplitude per layer (SFIL2, SFIL3): x_amp from ln(A_N/A_S) with the measured lambda_eff and g_N/g_S, both ends unclipped.
- Combination: inverse-variance weights from each estimator's rms against the paddle position; outputs per fit-cache shower in
  `event_x_walk.npz` (x_time[4], x_amp[4], x, xpad, astro_x).
- Result: 132-175 mm rms per shower, medians within 30 mm, 145 mm against the AstroPix x; the estimators are correlated and
  do not average down.

### 4.19 SFIL2 cell 3 TOA without the ADC timing (2026-10-28, `timing/toa_clean/`; supersedes 4.18)

- Why: 4.18 built the code-to-time table from the ADC 50 % crossing (most likely crossing phase per 16 codes) and kept only
  words within 6 ns of that crossing, so the TOA was not independent of the ADC. Its widths, time scale and the
  TOA/ADC/common split are withdrawn.
- ADC role: only the pulse selection (amplitude 150-950 ADC, peak in samples 5-14, raw maximum below 1000) and the 25 ns
  window (k = first sample above half the peak). Word: on k, else on k+1.
- Time: t = 25 s + T(code), s = the word's sample. Codes unwrapped at the empty code band (North c0 640; South none, c0 0).
  T = 25 ns x the cumulative code distribution (mid-bin) of pulse words on even triggers: uniform phase, because the beam is
  not locked to the 40 MHz clock. No time window; noise words stay in and form the tails.
- Optional correction: least-squares regression of the within-paddle N - S residual (even triggers, |r| < 10 ns) on ln A_S,
  ln A_N and both per-trigger pedestal shifts; applied to odd triggers. ln(A_N/A_S) also carries position (attenuation), so
  this can remove part of the slope; codes-only numbers are quoted alongside.
- N - S per paddle: Gaussian core plus flat background on +-25 ns, EM fit; slope from the four paddle means; sigma_x =
  core sigma / slope; separation = (mean P5 - mean P4) / average core sigma. The ADC crossing is processed identically on
  the same showers as an independent estimator.
- Good-TOA cuts (`toa_good.py`): quiet (no word on samples 1..k-1), one word (not both k and k+1), pure code (16-code bin
  purity = 1 - r N f_bg / h_pulse > 0.5, with r the per-sample word rate before the pulse, N the number of pulses, f_bg the
  normalised no-pulse code shape, h_pulse the pulse-word counts; even triggers), at both ends; the amplitude/baseline regression
  is refitted per stage on even triggers.
- Event time t0 (`toa_t0.py`, `caen_padtime.py`, figure `fig_t0.png`; user's question 2026-10-28 why t_N vs t_S is
  correlated, not anticorrelated): position moves (t_S, t_N) along -45 deg, the event time along +45 deg. The KCU latches
  the trigger on the 25 ns sample clock (ts0 = 3 mod 4), so the particle phase within the sample is not recorded.
  t0 = -res + t_pad, with res = 6.25 ts0 - t_CAEN from the run-3062 CAEN match (`sync/out/caen_kcu0_match.npz`; the trigger
  phase, spread over 25 ns) and t_pad the 50 % crossing of the largest CAEN paddle dip in the DRS4 window (particle time
  relative to the trigger, +-6.6 ns; P3/P4 triggers ~8 ns later than P5/P6). The KCU0 paddles cannot serve: a fired paddle
  drops to 0 ADC in one sample and writes no TOA word. On the 22,301 odd showers with a CAEN partner: (t_N + t_S)/2 against
  t0 slope 0.98 (phase alone 0.65, paddle alone 0.40); common spread 10.1 -> 4.6 ns (half 16-84 %); corr(t_N, t_S) 0.85 ->
  0.71 after removing t0; the remainder is per-end error (~2.4 ns, the half-difference), the CAEN 8.5 ns tick (~2.5 ns) and
  shared light timing. Within a paddle the beam spot gives tenths of a ns along -45 deg, invisible under ~2.4 ns. The
  anticorrelation shows across paddles: P4 -> P5, t_N - t0 goes from below to above the paddle's common time and t_S the
  other way (median N - S -2.43, -1.31, -0.49, +0.54 ns). A per-paddle common offset remains (P6 ~1.7 ns low on both ends).
  -45 deg view (`toa_m45_fig.py`, `fig_m45.png`): per-paddle 39 % contours and medians of (t_S - t0, t_N - t0); P4 -> P5 moves
  the centroid 2.13 ns along -45 deg and 0.76 ns along +45 deg. In `fig_tNtS.png` the edges and stripes sit on the 150/175/200 ns
  sample boundaries (t = 25 s + T(code), word on sample 6 or 7 for 99 %).

### 4.18 SFIL2 cell 3 TOA as a sub-sample TDC, redone (2026-10-27, `timing/toa_sfil2/`)

- Why 4.16 failed: it took the earliest word on samples ipk-2..ipk (often a pedestal-noise word: 14 % per sample South,
  36 % North), treated the code as linear from 0 (the North code wraps near 850, which made the double peak) and took
  per-code medians across a 25 ns frame jump.
- `toa_explore.py`: pulse words sit on the ADC-crossing sample k (57 % S / 53 % N) or k+1 (23 % / 20 %); on k-1 only at the
  pedestal rate. Code against crossing phase is one monotonic band per end, spanning the 25 ns frame over ~1000 codes
  (nominal 24.4 ps per code): the TOA is a working TDC on this channel.
- `toa_cdc2.py`: word on k or k+1 picked by consistency with the crossing; circular unwrap (c0 = 0 S, 528 N); code-to-phase
  table from the mode of the crossing phase per 16 codes (beats code density, which the noise words' fixed code bands
  distort); |t_TOA - t_cross| < 6 ns keeps 68 % (S) / 51 % (N) of pulses; walk in ln A (4 ns span).
- `toa_baseline.py`: the TOA threshold acts on the unsubtracted input, so the per-trigger common baseline offset (20-28 ADC rms)
  moves the TOA by -0.05 ns/ADC (corr -0.46); a (baseline, A) table built on even triggers, tested on odd, takes the per-end
  TOA - crossing error from 2.37 / 2.60 to 1.91 / 2.12 ns.
- Result (odd triggers, both ends good, 40 % of showers): sigma(t_N - t_S) TOA 2.1-2.5 ns, ADC 50 % crossing 1.56-2.25 ns, combined
  = ADC. Three-cornered split per bin of min(A): TOA-only 1.6-1.8 ns, ADC-only 0.8-1.2 ns, common 1.4-1.9 ns. The common term
  is the light timing both estimators see (beam spread inside a paddle is ~0.1 ns) and is the floor for this channel:
  ~100 mm; ADC alone gives 114-162 mm per shower. The TOA error does not fall with amplitude (not noise jitter); remaining
  candidates are first-photon statistics at a low threshold on the slow summed pulse and the ADC-referenced table (TOA
  v_eff 226 vs ADC 146 mm/ns: the table is compressed by the crossing-phase pinning).
- Horizontal position (`toa_waveforms.py`, needs `toa_baseline.npz`): x = (dt - b)/a with each estimator's own paddle fit
  (TOA 0.0089, ADC 0.0137 ns/mm); core sigma_x per paddle TOA 225 / 248 / 270 / 294 mm, ADC 112 / 121 / 134 / 149 mm (P4/P3/P6/P5);
  per bin of min(A) TOA 238-284 mm, ADC 114-165 mm, combined = ADC. Upper limits (beam spread within a paddle included).
  Figures for the TOA deck: `slides/wf_events.png` (single showers with the TOA time), `slides/wf_persist.png` (pulse shape
  aligned on the 25 ns tick vs on the TOA time), `slides/xres.png`; other slide figures from `toa_slides.py`.
- Against AstroPix (`toa_astropix.py`, odd triggers, purity >= 0.3, x = astro_x - 88.75): 2,379 showers, P3/P6/P5 only, each a
  5-10 mm spot (AstroPix medians -20 / +28 / +86 mm). Per-paddle mean x: TOA -11 / 12 / 87 mm, ADC -13 / 20 / 83 mm (scale fixed on
  the nominal paddle positions). sigma(x - x_AstroPix) TOA 251 mm, ADC 126 mm; within-paddle correlation consistent with zero
  (AstroPix spread far below the resolution).
