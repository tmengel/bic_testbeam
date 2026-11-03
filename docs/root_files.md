# Run 3062 raw data in ROOT trees (`/sphenix/tg/tg01/jets/tmengel/eic/`)

Every trigger (or record) of every stream, one entry per trigger in stream order, with everything the
raw packet carries, organised by subsystem. Symlinks in `data/`: `kcu0_3062.root`, `kcu1_3062.root`,
`caen_3062.root`, `astropix_3062.root`.

| File | Maker | Entries | Trees |
| --- | --- | --- | --- |
| `kcu0/kcu0_3062.root` | `decoding/kcu_to_root.cc` (packet 12001), 8 segment jobs + `kcu_merge` + `hadd` | one per KCU0 trigger | `header sfil1 sfil2 sfil3 sfil4 paddles bbcal unused packets map` |
| `kcu1/kcu1_3062.root` | same, packet 12002, 3 segments | one per KCU1 trigger | `header sfil5 sfil5north unused packets map` |
| `caen/caen_3062.root` | `decoding/caen_to_root.cc` (packet 2001), 3 segments + `hadd` | one per CAEN trigger | `header bbcal paddles packets` |
| `astropix/astropix_3062.root` | `decoding/astro_to_root.cc` (packet 13001) | one per hit record | `hits records packets` |

Condor: `sync/condor/*.job` (plain `condor_submit`, no DAGMan): `kcu_to_root.job` (11 segment jobs,
both boards), `caen_to_root.job`, `astro_to_root.job`; after they finish `kcu_merge.job` (board, nseg)
and `caen_merge.job`. Logs in `output/root_logs/`. Checks: `decoding/kcu0_root_check.py` ->
`kcu<b>_3062_check.txt`, `caen_3062_check.txt` next to the files. Per-segment files stay beside the
merged ones.


## Status (2026-10-22, all jobs done)

| File | size | entries | checks |
| --- | --- | --- | --- |
| `kcu0/kcu0_3062.root` | 18.9 GB | 3,766,023 triggers (3,324,643 complete = 88.3%; 188,809 single-sample fragments), 405,101 RCDAQ events | ts0 in order except 16 flagged backward steps, no duplicate ts0, sample steps all 164 ticks, RCDAQ sequence 2-405100 without gaps, 7 boundary duplicates removed |
| `kcu1/kcu1_3062.root` | 6.1 GB | 1,533,713 triggers (82.8% complete), 162,907 RCDAQ events | 5 backward steps, 33 of the first 200k triggers with a sample step other than 164 ticks (dead-board lines), 2 boundary duplicates removed |
| `caen/caen_3062.root` | 6.8 GB | 1,508,010 triggers | event counter +1 throughout, gtt monotonic, 0 decoder mismatches |
| `astropix/astropix_3062.root` | 60 MB | 2,864,554 hit records, 735,804 other records, 12,162 RCDAQ events | 0 hit-count mismatches |

Whole-run counters (packets tree): KCU0 171,303 UDP sequence jumps (55.6M datagrams lost), 7.85M
marker-less lines, 249 CRC failures, 95,862 Hamming-flagged lines; KCU1 76,266 jumps, 17.6M CRC
failures (dead board), 51,792 Hamming flags. The per-segment files beside the merged ones can be
deleted once the merged files are accepted (they double the space).

## KCU0 and KCU1

Waveforms are assembled with the stock decoder logic (ROC_pool `h2gcroc3_10G_pool`, local fork with
the timestamp-keyed assembly, 10-bit TOT, CRC and Hamming checks); the bytes of every line are then
looked up by (timestamp, link) and decoded field by field, so the file holds the raw words, not a
re-interpretation. `TTree::GetEntry(i)` on any of the per-trigger trees gives trigger i of the run;
the trees are aligned, so `header->AddFriend("sfil3")` etc. works.

## Per-trigger trees (same entry order; `[s]` = sample 0-19, `[l]` = link 0-7 = 2 x ASIC + half)

`header`

| branch | meaning |
| --- | --- |
| `seg`, `run`, `idx` | segment file, run number, entry index inside the segment's file |
| `ts0` | FPGA timestamp of sample 0 in 6.25 ns ticks = the trigger time (160 MHz counter, free running, not reset at run start; KCU1 has its own epoch 7.147 s later) |
| `nsamp`, `nlines`, `complete` | samples in the waveform (20 for a full readout), lines found (160 full), `nsamp == 20 && nlines == 160` |
| `nsamp_overflow` | the pool offered more than 20 samples (only the first 20 are stored) |
| `dts_prev`, `backward` | ts0 minus the previous entry's ts0 (ticks), and whether it went backwards |
| `evt_first`, `evt_last`, `seg_last` | RCDAQ event sequence of the first / last line's packet, segment of the last line |
| `evt_time` | RCDAQ event time (unix s) of the first packet |
| `rt_first`, `rt_last`, `rt_mono_first` | RTCLOCK packet (1010) of the first / last packet: host time in ns (`TOTAL`), and the monotonic clock |
| `udp_first`, `udp_last` | UDP datagram sequence numbers of the first / last line |
| `ts[20]` | per-sample timestamps (step 164 ticks = 1.025 us readout cadence; the ADC samples themselves are 25 ns apart) |
| `nlink[20]`, `linkmask[20]`, `present[20][8]`, `dts[20][8]` | which lines were found per sample; `dts` = line timestamp minus `ts[s]` (0 or +-4) |
| `trigger_in`, `trigger_out`, `event_counter` `[20][8]` | FPGA wrapper words w1-w3: 0 (unused); trigger counter, +1 per trigger; per-link sample counter |
| `w0`, `w6`, `w7` `[20][8]` | wrapper words 0 (marker 0xAA5A, board tag, ASIC, half) and 6-7 (spare, 0) |
| `roc_hdr[20][8]`, `bx`, `ec`, `orbit`, `ham`, `hdr_nibbles` | H2GCROC3 frame header word and its fields: bunch counter (12 bit), event counter (6 bit), orbit (3 bit), Hamming flags H1-H3, header/trailer nibbles (0xF5 when sane) |
| `cm_word`, `cm0`, `cm1` `[20][8]` | common-mode word and its two 10-bit ADCs |
| `calib_word`, `calib_adc`, `calib_toa`, `calib_tot`, `calib_tc`, `calib_tp` `[20][8]` | calibration channel word (w28) and fields |
| `crc[20][8]`, `crc_ok[20][8]` | the frame's CRC-32 word and whether it matches (poly 0x04C11DB7, init 0, over w8..w46) |
| `udp_seq`, `pkt_evt`, `seg_line` `[20][8]` | datagram sequence, RCDAQ event and segment of every line |

`sfil1` .. `sfil4`: per layer the eight live channels, branches `Sc1..Sc4` (South end, ASIC0) and
`Nc1..Nc4` (North end, ASIC2), each with `_ch` (channel number), `_adc[20]`, `_toa[20]`, `_tot[20]`,
`_tc[20]`, `_tp[20]`. Cell 3 is the beam row.

`paddles`: `P3_`, `P4_`, `P5_`, `P6_` with the same fields (KCU0 ASIC1 ch 84, 86, 81, 97; negative pulses).

KCU1 instead: `sfil5` with the four 16-SiPM arrays of SFIL5 South, `A_ch[16]`, `A_adc[16][20]`,
`A_toa`, `A_tot`, `A_tc`, `A_tp` and the same for `B`, `C`, `D` (ASIC2; vertical order from data
A < B < D < C); `sfil5north`: `ch[144]`, `adc[144][20]` ... for the SFIL5 North board on ASIC0/ASIC1
(dead in 3062: ch 0-35 and 72-107 railed); `unused`: the 80 channels without a detector.

`bbcal`: `ch[72]` (216..287), `adc[72][20]`, `toa`, `tot`, `tc`, `tp` (Baby BCAL North on ASIC3; no
signal in run 3062, ch 216 and 252-287 railed at 1023).

`unused`: `ch[180]` and the same arrays for the channels without a detector.

Channel-word conventions (characterization mode): `adc` 10 bit; `toa` 10 bit, nominal
fine[2:0] x 25 ps + mid[7:3] x 0.2 ns + coarse[9:8] x 6.25 ns (non-linear in practice, see method.md
4.9); `tot` the 10-bit compressed code, decompress with v < 512 ? v : 8 (v - 512) + 4; `tc` = TOT
complete flag (bit 31), `tp` = TOT in progress (bit 30). The 32-bit word is exactly
tc<<31 | tp<<30 | adc<<20 | tot<<10 | toa. **-1 in any field = that line (36 channels) is missing in
that sample** (UDP loss); `header.present[s][link]` says which.

## Per-segment trees

`packets` (main files, chained in order): one entry per RCDAQ event, also the non-data events:
`seg evt evt_type run evt_time evt_len has_kcu skipped pkt_len pkt_datalen rt rt_mono rt_diff
ndatagram udp_first udp_last udp_jumps udp_lost dg_const dg_tail_nonzero nlines nomarker badlink
dup_lines crc_bad ham_flagged tc_set tp_set ts_min ts_max nwf_emitted`. `nomarker` lines are the
zero-padded lines that UDP loss leaves (dropped by the decoder, counted here); `badlink` lines have
an ASIC/half outside 0-3 / 0x24-0x25.

`map` (once): one entry per channel: `ch asic half link k subsystem end layer cell label tree branch`.

`info_seg<k>` / `info_seg<k>_head` (TNamed): the job summary with the pool's counters, the ADC
cross-check (every 7th channel of every sample compared with the decoder's own value: 0 mismatches
expected) and the number of backward timestamp steps. `merge_note_*`: what the boundary step dropped.

## Segment boundaries

Each job reads its segment plus the first 6 packets of the next one, so the trigger straddling the
end of its segment is complete; a trigger whose sample 0 lies in the next segment is left to that
job. Triggers with a line in the segment's first packet go to `kcu<b>_3062_seg<k>_head.root`; for
k >= 1 `kcu_merge` drops from it the entries whose ts0 is not after the last sample timestamp of
segment k-1's main file (those are the tail of the straddling trigger, already written complete).
`hadd -ff` then chains head_0 main_0 head_1' main_1 ... head_7' main_7, which is time order. The
RCDAQ event sequence is continuous across segment files (segment 1 starts at event 53242 right after
segment 0's 53241), there are no repeated packets.

## CAEN DT5742 (`caen_3062.root`)

One trigger per RCDAQ event, nothing straddles, so the segments are simply chained. `header`: `seg run
idx evt evt_type evt_time rt rt_mono rt_diff`, the raw words `w0..w3` and their fields `marker
event_size board_id lvds_pattern group_mask w1_reserved event_counter w2_reserved event_time_tag`
(the 17.07 ns tag that steps back by 2^31: do not use), per group `g_hdr g_index_cell g_freq
g_tr_present g_sample_words g_hdr_reserved g_tt_word g_trigger_time g_tt_reserved` and the stitched
`gtt = g_trigger_time[0] | g_trigger_time[1] << 30` in 8.533 ns ticks (the clock to use), `ok_size`,
`ok_marker`. `bbcal`: `ch[12]`, `adc[12][256]` = Baby BCAL South, CAEN ch 0-11 (S A0, A1, A3, A4, A5,
6, 7, 8, 9, A14, 10, 11), 12-bit raw DRS4 cells at 1 GS/s (start cell in `g_index_cell`, no
correction applied), negative pulses from sample ~20, pedestal 3450-3530. `paddles`: `P3_adc[256]`
... = ch 12-15, pedestal 3720-3770, rail at 0 on large pulses. The converter cross-checks every 50th
trigger against the stock decoder (`info_seg<k>`).

## AstroPix (`astropix_3062.root`)

The stream is SPI records, not triggers; `hits` has one entry per 11-byte hit record in stream order:
`idx evt run evt_time rt irec ihit raw[11] layer id payload iscol flag location ts tot tot_nibble
fpga_ts valid`. `valid` = layer 1 and id 0-8; `flag` is byte 3 bit 6 and `tot_nibble` byte 5 bits
7:4, both dropped by the stock decoder and of unsettled meaning; `ts` is the 8-bit 10 MHz on-chip
stamp (same phase on all chips: ts = ((K + 4) // 16 - 189) mod 256 against KCU0 ticks K), `fpga_ts`
the 40 MHz readout stamp (wraps every 107 s, 19 or 88 ms after the trigger). A pixel hit is a row
record plus a column record with the same `ts`. `records` holds every other record (SPI filler and
the record cut by the end of a packet payload, `avail < len`) with its raw bytes; `packets` one entry
per RCDAQ event with the record-length histogram. Position: x = chip x 20 mm + col x 0.5 mm, chip 0
is the North end.

## Reading

```
TFile *f = TFile::Open("data/kcu0_3062.root");
TTree *h = (TTree *) f->Get("header"); h->AddFriend("sfil3"); h->AddFriend("paddles");
h->Draw("Nc3_adc[7]:Sc3_adc[7]", "complete && P4_adc[9] < 20");
```
Python: `uproot.open(path)["sfil3"].arrays(["Sc3_adc", "Nc3_adc"], entry_stop=10000)`.
