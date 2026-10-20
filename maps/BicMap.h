#ifndef BICMAP_H
#define BICMAP_H

// Channel map and clock constants for the Hall D BIC beam test, run 3062.
//
// PRIMARY SOURCE: the shift-crew hookup table "After 13/08 (Thu entry - beam
// studies - Jared + Halen), GOOD FOR RUNS 3061-3063". It supersedes the TIC
// slides wherever they disagree, because the boards were re-hooked on 8/13.
// Slides (HallD-BIC-09-15-2026.pdf) are used only for geometry/context; slide
// 13's KCU1 labels and slide 19's paddle channels are from the *earlier*
// hookup and are wrong for run 3062.
//
// Stream <-> board: the raw files beam_ROC_0-* / beam_ROC_1-* are KCU0 / KCU1.
// Verified in run 3062: ROC_0 chip0 lights up on exactly the 16 South channels
// of the table, ROC_0 chip2 on exactly the 16 North ones.
//
// KCU0 (= ROC_0), one H2GCROC chip per ASIC slot, 72 channels each:
//   ASIC0  ch   0- 71  SFIL 1-4 SOUTH, passive summing boards   (16 used)
//   ASIC1  ch  72-143  trigger paddles 3-6                      ( 4 used)
//   ASIC2  ch 144-215  SFIL 1-4 NORTH, passive summing boards   (16 used)
//   ASIC3  ch 216-287  Baby BCAL NORTH, BBS01-40 on J1-J40      (40 wired)
// ASIC2+ASIC3 sit on the HGCROC board that was swapped in from KCU1 on the
// night of 8/8-8/9 ("it magically started working").
//
// KCU1 (= ROC_1):
//   ASIC2  ch 144-215  SFIL5 SOUTH, individual SiPMs on Samtec, NO light
//                      guides and no cookies. This is the only live half of
//                      KCU1 and the only thing "half of SFIL5" in the run log
//                      refers to. 53 of ~64 SiPM channels respond.
//   ASIC0  ch   0- 71  SFIL5 NORTH, individual SiPMs, LGs + cookies
//   ASIC1  ch  72-143  SFIL5 NORTH (same board)
// ASIC0/ASIC1 are the board that came off KCU0; the table already says it "is
// not working well" and run 3062 agrees (see kDeadGroups below).
//
// Packet-level and clock facts come from the data itself, see README_TIMING.md.

namespace bic
{
  // ---- clocks -------------------------------------------------------------
  // H2GCROC/KCU FPGA timestamp: 160 MHz, triggers latched on the 40 MHz grid
  constexpr double kKcuTickNs = 6.25;
  constexpr int kBxTicks = 4;                  // 25 ns = one bunch = 4 ticks
  constexpr double kSampleNs = 25.0;           // waveform samples are one bunch apart
  // CAEN V1742 raw trigger time tag: increments every 2 periods of the
  // 117.1875 MHz reference, i.e. 17.0667 ns (measured 17.0662 ns)
  constexpr double kCaenTickNs = 2000.0 / 117.1875;
  // both boards run off the common clock card, so this ratio is exact
  constexpr double kKcuPerCaenTick = kCaenTickNs / kKcuTickNs;   // 2.730667
  // AstroPix FPGA timestamp: 40 MHz (slide 36)
  constexpr double kAstroTickNs = 25.0;
  // HGCROC TOA: 10 bits across one 25 ns bunch
  constexpr double kToaNsPerCount = kSampleNs / 1024.0;

  // ts(KCU1) - ts(KCU0) for run 3062, from the tick-space correlation peak.
  // The two FPGA counters are reset independently (no common reset), so this
  // is re-measured per run; kcu_align scans for it and only uses this as a
  // starting point.
  constexpr long long kOffsetRun3062 = -1143665946LL;

  // 10 Hz pulser into both CAEN and the KCUs (slide 33), empty-event signature
  constexpr double kPulserPeriodNs = 1e8;

  // ---- KCU0 ASIC0/ASIC2: SFIL 1-4, passive summing boards, both ends ------
  // [row][SFIL layer], row r = hookup-table cell r+1 ("SFIL <layer> <cell>").
  // SiPM labels per B. Kim "SiPM info" deck, 10/8 update (corrected SiPM3,4 ->
  // SiPM1,2), ANL LED bench test on SFIL2 North (J5..J8): cells 1 and 2 both
  // read SiPM1+SiPM2 mixed on v2 ("new") SiPM boards, cell 3 = SiPM3 (the beam
  // row), cell 4 = SiPM4. v2 boards: all South SFIL1-4 and North SFIL2-4;
  // North SFIL1 is v1 (one-to-one). Run 3062: r(cell1,cell2) = 0.97-0.99 on
  // v2, 0.20 on the v1 board. So rows 0+1 are one signal except North SFIL1.
  // SFIL4 South cells 1,2 sit on edge-board J63/J64; North cell 1 on J63,
  // cell 2 on J13.
  constexpr int kNRow = 4;
  constexpr int kNLayer = 4;
  constexpr int kSouth[kNRow][kNLayer] = {{2, 21, 63, 43},
                                          {6, 25, 59, 37},
                                          {10, 29, 52, 0},
                                          {14, 33, 48, 4}};
  constexpr int kNorth[kNRow][kNLayer] = {{146, 165, 207, 187},
                                          {150, 169, 203, 186},
                                          {154, 173, 196, 144},
                                          {158, 177, 192, 148}};
  // Trigger paddles on KCU0 ASIC1, in the 8/13 hookup:
  //   Paddle 3 -> J18 -> ch84     Paddle 5 -> J5 -> ch81
  //   Paddle 4 -> J4  -> ch86     Paddle 6 -> J6 -> ch97
  // Slide 19 lists the PREVIOUS hookup (82/86/93/97). Run 3062 settles it:
  // 84/86/81/97 fire at 43-59%, 82 and 93 are quiet.
  constexpr int kPaddleCh[4] = {84, 86, 81, 97};   // paddles 3, 4, 5, 6
  constexpr int kPaddle3 = 84, kPaddle4 = 86, kPaddle5 = 81, kPaddle6 = 97;

  // ---- KCU0 ASIC3: Baby BCAL North ---------------------------------------
  // BBS01-BBS40 on J1-J40. No channel numbers in the hookup table, and it does
  // not matter for run 3062: ch252-287 sit railed at ADC 1023 and ch216-251
  // have a clean pedestal but not one channel above 20 ADC. This is slide 16's
  // "Baby BCAL - Noise only -> Added CAEN on South Side". Use the CAEN stream
  // for Baby BCAL, never KCU0 ASIC3.
  constexpr int kBBCalNorthLo = 216, kBBCalNorthHi = 288;

  // ---- KCU1 ASIC2: SFIL5 SOUTH, individual SiPMs -------------------------
  // Four 16-SiPM arrays (Olaiya layout, Kim deck slides 6/31); every 9th channel
  // (152, 161, ..., 215) has no connector. 53 of 64 fire in run 3062 (A mostly dead).
  //   A = 144-151,153-160   B = 162-169,171-178   C = 180-187,189-196   D = 198-205,207-214
  // Vertical order from data (sync/sfil5_map.py; AstroPix y and the SFIL3-4
  // upper-cell fraction agree), bottom to top in AstroPix y: A < B < D < C
  // (B/D and D/C firm, A weak). No height structure is resolved within an array.
  // With the May-board analogy A=SiPM1, B=SiPM2, D=SiPM3, C=SiPM4 this is
  // monotonic but runs opposite to SFIL1-4 (there SiPM1+2 sit above SiPM3), so the
  // SFIL5 board is mounted flipped relative to SFIL1-4 if the analogy holds.
  constexpr int kKcu1ChipLo = 144;
  constexpr int kKcu1ChipMid = 180;
  constexpr int kKcu1ChipHi = 216;

  // ---- 36-channel groups sitting railed at ADC 1023 in run 3062 -----------
  // Slide 23's "half of channels, i.e. 36-group channels were having very high
  // pedestal, or were not responding", pinned down. Pedestal (sample 0) median
  // and MAD over 4000 waveforms; a railed group has median 1023 and MAD 0.
  //   KCU0: ch252-287 (ASIC3 grp1, Baby BCAL North)
  //   KCU1: ch0-35 (ASIC0 grp0) and ch72-107 (ASIC1 grp0), both SFIL5 North;
  //         their sister groups ch36-71 and ch108-143 have a normal ~80 ADC
  //         pedestal but only 5 channels (47,54,60,65,66) ever cross 20 ADC.
  //   KCU1: all of ASIC3 (nothing connected).
  // Everything else pedestals at 79-102 ADC with MAD ~1.

  constexpr int kChannels = 288;
}  // namespace bic

#endif  // BICMAP_H
