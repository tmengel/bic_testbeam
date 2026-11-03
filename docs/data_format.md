# Raw data format of the run 3062 streams: KCU (H2GCROC3), CAEN DT5742, AstroPix

Checked against the raw bytes on 2026-10-21 with the tools in `decoding/rawcheck/` (1500-3000 RCDAQ events
per stream: KCU0 segment 0001, KCU1 segment 0000, CAEN segment 0000 with 20k events, AstroPix
segment 0000). A longer version with a diagram is the Claude doc
"BIC Raw Data Format: KCU, CAEN and AstroPix Streams", https://claude.ai/artifact/HyG8TDatu9MYd2v6F78c3n.

Decoders: `/eic/data/purschke/eic/online_distribution/eventlibraries/oncsSub_idh2gcroc3_10G.cc`,
`oncsSub_idcaenv1742.cc`, `oncsSub_idastrov1.cc`, `oncsSub_idrtclock.cc`.

| Stream | Files (`data/raw/`) | Packet, hit format | RTCLOCK packet | RCDAQ events | RCDAQ event is |
| --- | --- | --- | --- | --- | --- |
| KCU0 | `beam_ROC_0-00003062-0000..0007.evt` | 12001, 302 `IDH2GCROC3_10G` | 1010 | 404,593 | a fixed-size buffer (~8.75 triggers) |
| KCU1 | `beam_ROC_1-00003062-0000..0002.evt` | 12002, 302 `IDH2GCROC3_10G` | 1011 | 162,905 | same |
| CAEN | `beam_CAEN-00003062-0000..0002.evt` | 2001, 85 `IDCAENV1742` | 1013 | 1,507,461 | one trigger |
| AstroPix | `beam_AstroPix-00003062-0000.evt` | 13001, 350 `IDASTROV1` | 1012 | 12,148 | a chunk of the SPI stream (~230 hits) |

## 1. RCDAQ layer (all streams)

- `fileEventiterator it(file, status); while (Event *e = it.getNextEvent())`; keep `DATAEVENT`.
  `e->getEvtSequence()`, `e->getPacket(id)`, `e->getPacketVector()`; delete what you get.
- `Packet`: `getIdentifier()`, `getHitFormat()`, `getLength()` (words), `getDataLength()`,
  `getPadding()`, `iValue(...)`, `lValue(...)`, `copyMe(buf, n)` (raw words, 16-byte packet header first).
- RTCLOCK (hit format 9): `lValue(0,"TOTAL")` = CLOCK_MONOTONIC_RAW ns when the RCDAQ process closed the
  event; `lValue(1,"TOTAL")` previous event; `lValue(2,"TOTAL")` CLOCK_MONOTONIC; `lValue(0,"TIMEDIFF")`.
  Shared host, but it is the buffer-close time (many triggers per KCU buffer; CAEN up to 2 s late).
- Files roll over by size (3 GB), not time.

## 2. KCU stream (H2GCROC3 via KCU105)

### Nesting

| Level | Size | Content |
| --- | --- | --- |
| RCDAQ packet | 68,304 words | 16-byte header + 200 UDP frames |
| UDP frame | 1366 bytes | 22-byte header + 7 lines |
| Line | 192 bytes = 48 big-endian words | one ASIC half (36 ch), one sample |
| Sample | 8 lines | 4 ASICs x 2 halves = 288 ch, one 25 ns crossing |
| Waveform | 20 samples = 160 lines | one trigger, 500 ns |

`ch = 72 x ASIC + 36 x half + k`, k = 0-35.

### UDP frame header (22 bytes)

| Bytes | Content |
| --- | --- |
| 0-3 | 0xAABBCCDD little-endian |
| 4-7 | frame length - 8, little-endian (1358) |
| 8-11 | frame sequence number, big-endian; +1 per frame, jumps of 27-297 at most buffer boundaries = UDP loss |
| 12-15 | constant 0xD0F9EE00 (KCU0) |
| 16-21 | 0 |

### Line, bit allocation (w0-w47, big-endian; bit 31 first on the wire)

FPGA wrapper:

| Word | Bits | Field | Observed |
| --- | --- | --- | --- |
| w0 | 31:16 | marker | 0xAA5A |
| w0 | 15:12 | board/link tag | 0 on KCU0, 1 on KCU1 |
| w0 | 11:8 | ASIC | 0-3 |
| w0 | 7:0 | half | 0x24 (ch 0-35), 0x25 (ch 36-71) |
| w1 | 31:0 | trigger_in | 0 |
| w2 | 31:0 | trigger_out | FPGA trigger counter, +1 per trigger, gaps where triggers were lost |
| w3 | 31:0 | event counter | per-link L1A counter, +1 per sample, differs between links |
| w4-w5 | 63:0 | FPGA timestamp | 6.25 ns ticks, +164 per sample |
| w6-w7 | 63:0 | spare | 0 |

ASIC frame:

| Word | Bits | Field | Observed |
| --- | --- | --- | --- |
| w8 | 31:28 | header | 0xF |
| w8 | 27:16 | BX | 1-3564, +1 per sample, own phase per ASIC (not a common clock) |
| w8 | 15:10 | event counter | 6 bit, +1 per sample |
| w8 | 9:7 | orbit counter | 3 bit, uniform |
| w8 | 6:4 | error/Hamming flags | 0 in 99.97% |
| w8 | 3:0 | trailer | 0x5 |
| w9 | 19:10, 9:0 | common mode 0 and 1, 10-bit ADC each | KCU0: 30-80 and 0/1; KCU1: 50-60 and 137/139/261... |
| w10-w27 | | channel words k = 0-17 | |
| w28 | | calibration channel, same layout as a channel word | ADC 48-77 |
| w29-w46 | | channel words k = 18-35 | |
| w47 | 31:0 | CRC-32 | not identified, not checked by any decoder |

Channel word (characterization mode: ADC, TOT, TOA every crossing):

| Bits | 31 | 30 | 29:20 | 19:10 | 9:0 |
| --- | --- | --- | --- | --- | --- |
| Field | Tc (TOT complete) | Tp (TOT in progress) | ADC | TOT (10-bit compressed code) | TOA |
| Stock decoder | dropped | dropped | `iValue(wf,s,ch)` | `(w>>10)&0x3f`: **low 6 bits only** | `iValue(wf,s,ch,"TOA")` |

### Stock decoder behaviour (`oncsSub_idh2gcroc3_10G`)

1. Walks UDP frames, visits lines starting 0xAA5A.
2. A line within 5 ticks of the previous one joins the current sample; otherwise a new sample starts.
3. `parse_timeline`: a step > 164 ticks starts a new waveform.
4. Decodes one RCDAQ event at a time: waveforms crossing a buffer boundary become two fragments, the
   last sample of each buffer is dropped, a link running a sample ahead fragments the waveform.
   ~60% of KCU0 triggers give a complete 20-sample waveform. Use `SAMPLESIZE == 20` for timing.
5. The first buffer of each file is a stale ts ~ 0 flush (KCU0 13243, KCU1 35713).

Accessors: `iValue(0,"NR_WF")`, `iValue(wf,"SAMPLESIZE")`, `iValue(0,"CHANNELS")`, `iValue(wf,s,ch)` ADC,
`iValue(wf,s,ch,"TOA"|"TOT"|"TRIGGER_OUT"|"EVENT_COUNTER"|"TRIGGER_IN")`, `lValue(wf,s,ch)` timestamp,
`iValue(n,"TIMESTAMP"|"CMN"|"CALIB")` per packet-wide sample n (last line of the sample).
Local build only (fork, 2026-10-21): `iValue(wf,s,ch,"CRCERROR")` = 1 if the line that carried `ch`
failed its CRC-32, `iValue(wf,s,ch,"ERROR")` = the header's H1 H2 H3 flags for that line,
`iValue(0,"CRCERRORS")` / `iValue(0,"HAMMINGERRORS")` = counts of such lines in the packet; `ddump`
prints both counts. Lines are not dropped. `"TOT"` is the full 10-bit code in the local build.
Also local: `iValue(0,"BADFRAMES"|"BADLINES"|"LATELINES"|"RESETS"|"DECODEERRORS")`, and samples are
assembled by timestamp (a line joins the open sample within 5 ticks; samples close 4 x 164 ticks behind
the newest timestamp), so out-of-order lines no longer split waveforms; steps 4-5 of the stock
behaviour above describe Purschke's build.

### Samples

Sampling is 25 ns (BX and the 6-bit counter step by exactly 1 per sample) although the timestamp
steps by 164 ticks. Triggers sit on a 4-tick grid. Peak at sample 6 (KCU1) or 7 (KCU0); pedestals
79-102 ADC, MAD ~1.5; railed groups read 1023 with MAD 0; SiPM pulses saturate at 1023; the R580
paddles on KCU0 ASIC1 are negative and saturate at 0.

### TOA and TOT (run 3062)

TOT: 10-bit HGCROC3 compressed code of a 12-bit value: v < 512 as is, else TOT = 8 (v - 512) + 4
(codes 528-575 never occur). Fires on 6.7% (KCU0) / 10.8% (KCU1) of hit channels, almost only when
the ADC saturated, in one sample at the peak (6-7), with Tc set. Overflow code 1023 = 13.7% of KCU0
TOT hits. The stock value equals code & 0x3F in every case checked (27k) and is a sawtooth; the
"cap at 63" is 1023 & 0x3F. Full code vs number of saturated samples: Spearman 0.74 (KCU0), 0.89
(KCU1) against 0.50/0.20 for the 6-bit value; ~55-60 ps per count above 4 saturated samples. To use
it, read bits 19:10 of the raw word (`decoding/rawcheck/kcu3.cc`).

TOA: the 10-bit TDC word, fine = bits 2:0 (FTDC, 25 ps LSB), mid = bits 7:3 (CTDC, 0.2 ns), coarse =
bits 9:8 (160 MHz counter, 6.25 ns); stored in the sample AFTER the one in which the ADC crosses the
pulse edge. In run 3062 the TOA threshold was set per channel on none of KCU0's SFIL channels (21 in
the noise, 10 above the signal, 1 right: N3c3 ch196) and on all of KCU1's. Where the threshold is
right the TDC works: the coarse bits are the 6.25 ns quarter (uniform), the mid bits the time within
it (two channels of one waveform agree to 4.4 codes), with a per-channel non-linear mapping. A
self-calibration across channels (`timing/toa_selfcal.py`) reaches ~190 ps per channel on KCU1; the
ADC crossing gives 1.6-2.0 ns per end for North-South on KCU0 (summing-board pulses) and is pinned to
the sampling phase on the fast single-SiPM pulses, so it is no reference there. Per-channel
thresholds (and ideally a delay scan) are needed before North-South ps timing exists. Details:
`method.md` 4.9, `timing/toa_timing.py`.

## 3. CAEN DT5742 (packet 2001, 1544 words, one trigger)

| Word | Bits | Field | Accessor |
| --- | --- | --- | --- |
| w0 | 31:28 / 27:0 | 0xA / event size 1544 | checked / - |
| w1 | 31:27 / 23:8 / 3:0 | board 0x1F / LVDS 0 / group mask 0x3 | `"LVDSPATTERN"`, `iValue(g,"GROUPPRESENT")` |
| w2 | 21:0 | event counter | `iValue(0,"EVNR")` |
| w3 | 31:0 | event time tag, 17.07 ns, steps back by 2^31 | `iValue(0,"EVENTTIME")` **do not use** |
| g+0 | 29:20 / 17:16 / 12 / 11:0 | index cell / freq code 2 = 1 GS/s / TR present 0 / 768 sample words | `iValue(g,"INDEXCELL")`, `"FREQUENCY"`, `"TRPRESENT"` |
| g+1..g+768 | | 256 slices x 3 words, 8 ch x 12 bit | `iValue(sample, ch)` |
| g+769 | 29:0 | group trigger time, 8.533 ns | `lValue(0,"GROUPTRIGGERTIME")` |

g = w4 (group 0), w774 (group 1). Slice packing: ch0 W0[11:0], ch1 W0[23:12], ch2 W0[31:24]+W1[3:0]<<8,
ch3 W1[15:4], ch4 W1[27:16], ch5 W1[31:28]+W2[7:0]<<4, ch6 W2[19:8], ch7 W2[31:20].
Group 1's trigger-time word carries the high rollover bits (0 early, 0x21 later); the decoder forms
`GTT0 | GTT1 << 30`, monotonic over the run. Channels 0-11 Baby BCAL South, 12-15 paddles P3-P6;
negative pulses from sample ~20; pedestal = mean of samples 0-19; `qint` = sum over 20-159 of
(ped - ADC). Samples are raw DRS4 cells (start cell available; no correction applied).

### CAEN bit map (little-endian words; observed on 3000 events)

| Word | Bits | Field | Observed |
| --- | --- | --- | --- |
| w0 | 31:28 / 27:0 | marker / event size | 0xA / 1544 |
| w1 | 31:27 / 26:24 / 23:8 / 7:4 / 3:0 | board ID / reserved / LVDS pattern / reserved / group mask | 0x1F / 0 / 0 / 0 / 0x3 |
| w2 | 31:22 / 21:0 | reserved / event counter | 0 / +1 per event |
| w3 | 31:0 | event time tag, 17.07 ns | 2^31 back-steps; do not use |
| g+0 | 31:30 / 29:20 / 19:18 / 17:16 / 15:13 / 12 / 11:0 | reserved / index cell / reserved / freq code / reserved / TR present / sample words | 0 / 0-1023 (g1 = g0 +-1) / 0 / 2 / 0 / 0 / 768 |
| g+769 | 31:30 / 29:0 | reserved / group trigger time | 0 / g0 counts, g1 = rollover count (0 early, 0x21 later) |

Slice: W0 = [k2 7:0 | k1 | k0]; W1 = [k5 3:0 | k4 | k3 | k2 11:8]; W2 = [k7 | k6 | k5 11:4] (12-bit
fields, channel 8 x group + k). Pedestals 3450-3530 (Baby BCAL) and 3720-3770 (paddles) of 4095,
pulses negative; paddle channels reach 0 on large pulses.

## 4. AstroPix (packet 13001, byte-oriented SPI chunks)

Records `[L][layer][L-5 SPI bytes][fpga_ts 4 B big-endian]`, L + 1 bytes. Only 11-byte records are
hits; 7-10 and 12-14-byte records are SPI filler (0xBA, 0xD9, 0x1A, 0x21, 0x1E ...). The decoder
keeps 11-byte records only, does no validity check, and its multiset comparator never compares
`fpga_ts`, so hits come back in arrival order.

| Byte | Bits | Field | Accessor | Valid |
| --- | --- | --- | --- | --- |
| 1 | 7:0 | layer | `"LAYER"` | 1 |
| 2 | 7:3 / 2:0 | chip id / payload | `"ID"` / `"PAYLOAD"` | 0-8 / 4 |
| 3 | 7 / 6 / 5:0 | row-or-column flag / undecoded / location | `"ISCOL"` / - / `"LOCATION"` | 0-34 (2% read 35-63) |
| 4 | 7:0 | on-chip timestamp, 10 MHz, 8 bit | `"TIMESTAMP"` | |
| 5-6 | 3:0 + 7:0 | ToT, 12 bit (byte 5 high nibble undecoded) | `"TOT"` | median 586, 21% zero |
| 7-10 | 31:0 | FPGA readout stamp, 40 MHz, wraps 107 s | `lValue(i)` | unwrap |

A pixel hit = a row record + a column record with the same on-chip timestamp. Clustering: ts within
+-2 of the first hit, read within 1 ms (2.35M clusters, 397k row+col pairs).

### AstroPix hit-record bit map (bit 7 first; observed on 693k valid hits)

| Byte | 7 | 6 | 5 | 4 | 3 | 2 | 1 | 0 | Observed |
| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |
| 0 | length - 1 | | | | | | | | 0x0A |
| 1 | layer | | | | | | | | 0x01 |
| 2 | id4 | id3 | id2 | id1 | id0 | pl2 | pl1 | pl0 | id 0-8, payload 4 |
| 3 | iscol | flag | loc5 | loc4 | loc3 | loc2 | loc1 | loc0 | bit 6 set 1.1%, then 37% loc > 34 |
| 4 | ts7 | ts6 | ts5 | ts4 | ts3 | ts2 | ts1 | ts0 | 10 MHz counter, uniform |
| 5 | x3 | x2 | x1 | x0 | tot11 | tot10 | tot9 | tot8 | nibble x non-zero 7.9%, then ToT != 0 |
| 6 | tot7 | tot6 | tot5 | tot4 | tot3 | tot2 | tot1 | tot0 | |
| 7-10 | fpga_ts 31:24, 23:16, 15:8, 7:0 | | | | | | | | 40 MHz, uniform |

Byte 3 bit 6 and byte 5 bits 7:4 are dropped by the stock decoder; corrupt records keep layer 1 and
payload 4 but have id 9-31. The meaning of byte 5's upper nibble is unsettled (check the A-STEP
firmware hit format before using ToT quantitatively).

## 5. Clocks

| Clock | Tick | Width | Relation to KCU0 (run 3062) |
| --- | --- | --- | --- |
| KCU0 timestamp | 6.25 ns | 64 bit | reference |
| KCU1 timestamp | 6.25 ns | 64 bit | -1,143,665,946 ticks, residual 0 (93%) or +4 |
| ASIC BX | 25 ns | 12 bit | own phase per ASIC; not usable across ASICs |
| CAEN GROUPTRIGGERTIME | 8.533 ns | 30 + 30 bit | 1.365340679 KCU ticks per tick + running median; MAD 6.4 ns |
| CAEN EVENTTIME | 17.07 ns | 32 bit | broken (2^31 steps) |
| AstroPix on-chip ts | 100 ns | 8 bit | ts = ((K + 4) // 16 - 189) mod 256, same phase on all 9 chips |
| AstroPix fpga_ts | 25 ns | 32 bit | readout stamp, 19 or 88 +- 5 ms after the trigger |
| RTCLOCK | 1 ns | 64 bit | coarse anchor only |

All counters run from one 40 MHz card with no common reset; measure every offset per run.

## Begin-run packets 900 / 910 / 920 (checked 2026-10-24, `ddump -p 910 -t 9 beam_ROC_0-00003062-0000.evt`)

Event type 9 of every ROC file holds three text packets (IDCSTR): **900** = the DAQ setup script (`setup_all_oneCAEN.sh`, 3174 bytes,
identical on both boards), **920** = a power-supply table (channels U500-U907 with set and measured voltages and currents; U900-U903 are
5 V at 0.75-0.77 A, U708-U715 and U712/U713 are set to 73.5-75 V and read back, the others off), **910** = the KCU/H2GCROC configuration
(782 bytes): the FPGA frame decoded as DAQ_FCMD 0x4B, GEN_PRE_FCMD 0x2D, EXT TRIG ENABLE 1, EXT TRIG DELAY 0, **EXT TRIG DEADTIME 4000**
(x 25 ns = the 100 us busy seen in the trigger spacing), **MACHINE_GUN 19** (20 samples per trigger), GEN interval 1e6 / 100 cycles,
ASIC enables 0x30 for ASIC0-7 (D1 = D0 = 1, T = 0), followed by a sample of ASIC register reads: `Top`, `CM_2`, `CM_3` and
`Channel_36`, `Channel_37`, `Channel_38` of ASIC0 only. The two boards differ only in the FPGA address (0 / 1) and in those three channel
registers (ROC_0 ch36 reg0 = 0xac: Inputdac 44, Gain_conv<1> set; ch38 reg0 = 0x64; ROC_1 ch36 reg0 = 0x00, ch38 = 0xfc). **The packet does
not contain the Global_Analog block (Cf, Rf, Gain_conv<3>) or the channel registers of the live SFIL5 South ASIC (ASIC2 of KCU1), so the
preamp gain settings of the two boards cannot be compared from the data file**; the YAML/configuration used by the DAQ is needed for that.
