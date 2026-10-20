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
- **Result:** 49.7% of CAEN triggers match within 100 ns, with residual MAD 6.4 ns. The +10 ms control
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
  `ts == ((K + 4) // 16 - 189) mod 256`. The evidence is a 19.3 sigma spike against 3.0 sigma with
  `ts` shuffled between clusters; the divisors 32-640 are flat.
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
| KCU0 + KCU1 + CAEN | 155k | 3.7% | |
| KCU0 + CAEN | 594k | 14.2% | |
| KCU0 + KCU1 | 457k | 10.9% | |
| KCU0 only | 1795k | 42.8% | |
| KCU1 + CAEN | 120k | 2.9% | 32k |
| KCU1 only | 431k | 10.3% | 113k |
| CAEN only | 639k | 15.2% | 128k |

- **Rates:** the union is 4.19M triggers over 1486 s (2819 Hz). By stream: KCU0 complete 2018 Hz
  (2202 Hz with fragment-only triggers), KCU1 782 Hz, CAEN 1014 Hz.
- **No common busy:** the streams share the trigger but not a busy/veto, so each keeps a different
  subset.
- **KCU1 vs CAEN timing:** on the 120k triggers without KCU0, KCU1 and CAEN agree to MAD 6.6 ns. This
  confirms the KCU1 offset independently of KCU0.
- **Paddle pattern:** 89.7% of triggers have one (71.6% from KCU0, 18.1% from CAEN only).
- **KCU0 fragments:** a trigger whose KCU0 waveform was split by the decoder leaves fragments at
  trigger + k x 164 ticks (k = 0..19, within +-4 ticks for 96%).
  - They are flagged (`kcu0_frag`) but not yet reassembled.
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

- **Rows:** one entry per common trigger (4,189,768), in time order. A stream that missed the trigger has
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
- **AstroPix match counts:** 178k triggers carry a match with log10LR > 1, about 76k of them true.
  Cut on `astro_purity` (20k triggers at > 0.8).
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

