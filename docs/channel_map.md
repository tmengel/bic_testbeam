# Channel mapping, run 3062


### 3.1 Authority and comparison with `roc_analysis/channel_mapping.txt`

`channel_mapping.txt` is the shift-crew table "After 13/08 (Thu entry - beam studies - Jared +
Halen), GOOD FOR RUNS 3061-3063". **It is the map used here, and it matches the data channel for
channel**, with these qualifications:

- **Paddles:** the table has two "On ASIC" columns. The first (84, 86, 81, 97) is the 8/13 hookup
  and is correct for 3062. The second (82, 86, 93, 97) is the pre-8/13 hookup, also on TIC slide 19.
  In 3062, ch 84/86/81/97 fire at 43-59%, and 82 and 93 are quiet.
- **Paddle 5's connector:** the table lists it on J5, but in the edge-board table (Kim deck p17, run log
  p60) J5 = ASIC1 ch 93 and J17 = ch 81. The channel (81) is right and the connector label is a typo.
- **KCU1 (bottom of the table):** A2 = SFIL5 South (works); A0/A1 = SFIL5 North on the board moved
  from KCU0 (dead in 3062).
- **Rows not in the table:** it gives no SFIL5 per-SiPM channels and no CAEN map. Both are filled in
  below from other sources and from data.

Precedence when sources disagree: data > 8/13 hookup table > Kim "SiPM info" deck (10/8 version) >
run-log mapping pages > TIC slides. TIC slides 13 and 19 describe the pre-8/13 setup.

### 3.2 KCU0 (ROC_0)

| ASIC | Channels | Detector | State in 3062 |
| --- | --- | --- | --- |
| A0 | 0-71 | SFIL1-4 South, passive summing boards | 16 live |
| A1 | 72-143 | trigger paddles P3-P6 | 4 live |
| A2 | 144-215 | SFIL1-4 North, passive summing boards | 16 live; SFIL1 North column noisy (~99% above 20 ADC) |
| A3 | 216-287 | Baby BCAL North (BBS01-40 on J1-J40) | no signal: ch 252-287 railed at 1023; ch 216-251 clean pedestal, nothing above 20 ADC |

A2+A3 are the board swapped in from KCU1 on 8/13.

SFIL1-4 grid (`kSouth` / `kNorth` in `maps/BicMap.h`, row = hookup cell - 1):

| Cell | SiPM (Kim deck, 10/8 update) | SFIL1 S/N | SFIL2 S/N | SFIL3 S/N | SFIL4 S/N |
| --- | --- | --- | --- | --- | --- |
| 1 | SiPM1+SiPM2 mixed (v2 boards) | 2 / 146 | 21 / 165 | 63 / 207 | 43 / 187 |
| 2 | SiPM1+SiPM2 mixed (v2 boards) | 6 / 150 | 25 / 169 | 59 / 203 | 37 / 186 |
| 3 | SiPM3, **beam row** | 10 / 154 | 29 / 173 | 52 / 196 | 0 / 144 |
| 4 | SiPM4 | 14 / 158 | 33 / 177 | 48 / 192 | 4 / 148 |

- **v2 mixing:** v2 ("new") SiPM boards mix SiPM1 and SiPM2, so cells 1 and 2 carry the same signal.
  The v2 boards are all South SFIL1-4 and North SFIL2-4. North SFIL1 is v1 (one-to-one). In data,
  r(cell1, cell2) = 0.97-0.99 on v2 and 0.20 on v1, so treat cells 1+2 as one row.
- **Kim deck correction:** the 9/29 version of the deck said SiPM3/SiPM4. The 10/8 update corrected it
  to SiPM1/SiPM2 (bench test: LED on SiPM1 or SiPM2 lights J5 and J6; SiPM3 lights J7; SiPM4 lights J8).
- **No physical up/down claims:** the SiPM boards were flipped between May and August, so don't state
  which way is up.
- **Edge-board connectors:** SFIL4 South cells 1,2 are on J63/J64; North cell 1 is on J63 and cell 2 on J13.
- **SFIL4 gain:** SFIL4 has a large North/South gain offset (beam-row (N-S)/(N+S) = -0.54).
- **SFIL5 North (8/14):** the 8/14 change put SFIL5 North on KCU0 A2 at J17-J20 = ch 153, 156, 163,
  167. It is **not** in run 3062 (those channels fire 0.1%).

| Paddle | KCU0 ch (connector) | CAEN ch | Position along AstroPix module (N to S) |
| --- | --- | --- | --- |
| P4 | 86 (J4) | 13 | chip 0 |
| P3 | 84 (J18) | 12 | chips 2-3 |
| P6 | 97 (J6) | 15 | chips 5-6 |
| P5 | 81 (J17) | 14 | chip 8 |

**Paddle signal shape on KCU0 (checked 2026-10-20, `ts_dump kcuwf`, 200k waveforms of segment 0001).**

- **Shape:** the R580 pulse is negative. A hit drives the H2GCROC ADC to **0 (saturated)** for ~9 samples,
  starting at sample 8 or 9, followed by an overshoot at samples 17-19. Paddles that did not fire show only
  a +40-60 ADC crosstalk spike at sample 8.
- **Hit flag:** `>= 3 samples <= 5 ADC in samples 8-16`. The flag carries no amplitude information.
- **The old `amp` is wrong for paddles:** max - pedestal (`kcuadc` amp[32:36]) mixes the crosstalk spike with
  the overshoot. It is a poor proxy, and the source of the 43-59% "fire fractions" quoted before.
- **Every trigger has a paddle (KCU0 and CAEN agree, so the paddle OR triggers both):**
  - A paddle is saturated on 99.8% of KCU0 triggers, with exactly one paddle on 94%.
  - Per paddle P3-P6: 0.243 / 0.226 / 0.329 / 0.259.
  - CAEN ch 12-15 (neg > 200) gives 0.248 / 0.230 / 0.335 / 0.265, 99.8% any.
  - On 41.7k matched CAEN-KCU0 pairs, the 4-paddle pattern is identical in 98.7% (P(same paddle) = 0.999),
    against 22% for decoys.
- **In the pipeline:** `kcuadc` now stores the dip for paddles: amp[32:36] = pedestal - min(samples 8-16),
  about 80 for a hit. The common-trigger table uses dip > 65, the best agreement with CAEN.
  - **Full run:** the KCU0 and CAEN patterns are identical on 96.1% of 749k shared triggers, against 21.8%
    for the neighbouring trigger. P(same paddle in CAEN | KCU0 paddle) = 0.96-0.98.
  - **Two-paddle triggers:** KCU0 has more of them than CAEN (8.6% against 5.8%).
  - The old max - pedestal files are kept in `sync/out/oldpad/`.

### 3.3 KCU1 (ROC_1): SFIL5 South, individual SiPMs (v1 board, no light guides)

| 16-ch array | Channels | In 3062 | Height (data, bottom = 1) | May-board analogy |
| --- | --- | --- | --- | --- |
| A | 144-151, 153-160 | mostly dead (only 154, 156-160 fire) | 1 (weak evidence) | SiPM1 |
| B | 162-169, 171-178 | live, strong | 2 | SiPM2 |
| D | 198-205, 207-214 | live, strong | 3 | SiPM3 |
| C | 180-187, 189-196 | live, weak (~10-20% share) | 4 | SiPM4 |

- **Gaps:** ch 152, 161, 170, 179, 188, 197, 206 and 215 are not connected (Olaiya's layout,
  Kim deck p31). 53 of 64 SiPMs fire.
- **Vertical order from data** (`sync/sfil5_map.py`, `out/sfil5_map.png`, 2026-10-20): bottom to top
  in AstroPix y, the order is **A < B < D < C**.
  - **AstroPix (background-subtracted):** B and D cross near y = 7 mm. D peaks at y ~ 12 mm and turns
    down while C keeps rising. A's share is low and roughly flat, highest at the lowest y.
  - **SFIL3-4 upper-cell fraction** (36k triggers, no AstroPix): as the shower moves up, A and B lose
    share (log top/bottom -0.32, -0.46) and C and D gain (+0.51, +0.31).
  - **Confidence:** B/D and D/C are firm; A at the bottom is weaker (6/16 live).
  - **Response-weighted centroids** (`dy`) put A and C between B and D. Do not use them for the outer
    arrays: those respond non-monotonically over the 17.5 mm AstroPix range.
  - **Calibration:** the same estimator on SFIL1-4 puts cells 1+2 at +1 to +2.5 mm above cell 3
    (-0.5 to -1.1). Cell 4 is mostly flat (no signal); SFIL4 South cell 4 is slightly below cell 3.
    So in SFIL1-4, SiPM number increases downward in AstroPix y.
- **Within an array:** all 16 channels move together. Per-channel spread is chi2/ndf ~3-5 (gains?),
  with no sub-rows.
- **Array-to-SiPM assignment:** the May-board analogy (A=SiPM1, B=SiPM2, D=SiPM3, C=SiPM4) matches the
  measured order with SiPM number increasing upward. That is opposite to SFIL1-4, so the SFIL5 South
  board is mounted flipped relative to SFIL1-4, if the analogy holds.
  - The beam sits at the B/D (SiPM2/SiPM3) boundary in SFIL5, but inside cell 3 in SFIL1-4.
  - A bench check of the pin map or orientation would still settle the SiPM numbers; the analysis only
    needs the measured order.
- **SFIL5 North (A0/A1):** dead. Ch 0-35 and 72-107 are railed at 1023. Ch 36-71 and 108-143 have
  normal pedestals, but only ch 47, 54, 60, 65 and 66 ever fire.
- **A3:** unconnected.

### 3.4 CAEN DT5742 (run log, "Changed on Aug 10")

ch 0-11 = Baby BCAL South (S A0, A1, A3, A4, A5, 6, 7, 8, 9, A14 -> 10, 11). ch 12-15 = paddles
P3-P6. Pulses are negative, 256 samples at 1 GS/s, groups 0-1. The paddle assignment was found
independently in data: r = 0.26/0.29/0.40/0.30 against the KCU0 paddles over 749k matched triggers.

### 3.5 AstroPix (9-chip v3 module, layer 1)

- **Valid hits:** only `layer==1` and `id<=8` are real; the rest are corrupt packets.
- **Orientation:** chip 0 is the North end and chip 8 the South end. Columns (`iscol=1`) run along
  the module in the chip-index direction; rows (`iscol=0`) are transverse.
- **Position:** x = chip x 20 mm + col x 0.5 mm, using nominal pitch. The beam band sits at
  rows 12-22.
- **Occupancy:** columns 0-2 are nearly empty. Chips 3, 5, 6 and 8 carry most trigger-matched hits.

