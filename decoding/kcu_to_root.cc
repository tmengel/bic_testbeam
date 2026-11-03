// kcu_to_root: every trigger of one KCU segment into ROOT trees, one entry per trigger, in stream
// order, organised by subsystem. Nothing of the raw line is left out: the FPGA wrapper words, the
// H2GCROC3 frame header, common-mode and calibration words, CRC word and check, and every channel's
// ADC / TOA / TOT / Tc / Tp for every sample. Waveforms are assembled by the stock decoder logic
// (ROC_pool, h2gcroc3_10G_pool); the raw bytes of every line are looked up by (timestamp, link).
//
//   kcu_to_root out.root segment.evt seg [next_segment.evt] [n_head_packets=6] [packetid=12001|12002] [max_packets=0] [run=0]
//
// run selects the hookup epoch of the run-log sheet (ASIC_Mapping): run < 3085 (or 0) = the 13/08 table (runs 3061-3069);
// run >= 3085 = "After 14/08, changed on Aug 14 9:55 PM CT": the SFIL5 North summed boards are on KCU0 ASIC2 (J17-J20 =
// ch 153, 156, 163, 167 = SFIL5 cells 1-4 North), written to an extra KCU0 tree sfil5n (Nc1..Nc4) instead of unused.
//
// Packet 12001 = KCU0 (SFIL1-4 both ends, paddles, Baby BCAL North), 12002 = KCU1 (SFIL5 South arrays
// A-D on ASIC2, the dead SFIL5-North board on ASIC0/1, nothing on ASIC3).
//
// Output: out.root with the trees of the triggers whose lines all come from packets after the first
// one, and out_head.root with the triggers that have a line in the segment's FIRST packet (the
// continuation of the trigger that straddles the segment boundary lands there; kcu_merge drops it
// when the previous segment already holds that trigger complete). The next segment's first
// n_head_packets packets are read too, so the trigger straddling the end of this segment is complete;
// triggers whose sample 0 lies in the next segment are left to that segment's job.
//
// Trees (all one entry per trigger, same order; AddFriend them or GetEntry(i) on each):
//   header   stream / RCDAQ / FPGA / frame-header information, per sample and per line
//   KCU0: sfil1..4 the 8 channels of each SFIL layer: {S,N}c{1..4}_{adc,toa,tot,tc,tp}[20], _ch
//         paddles  P3 P4 P5 P6;  bbcal  ASIC3 (Baby BCAL North), 72 channels: adc[72][20] ...
//         unused   the 180 channels without a detector: adc[180][20] ...
//   KCU1: sfil5 the four 16-SiPM arrays of SFIL5 South: {A,B,C,D}_ch[16], {A,B,C,D}_{adc,toa,tot,tc,tp}[16][20]
//         sfil5north the SFIL5 North board on ASIC0/1 (dead in 3062): ch[144], adc[144][20] ...
//         unused   the 80 channels without a detector (ASIC2 gaps, ASIC3)
//   packets  one entry per RCDAQ event of the segment (also the non-data events)
//   map      one entry per channel: where it is and which tree/branch holds it (segment 0 only)
// A missing line (UDP loss) leaves its 36 channels at -1 in that sample; header.present[s][link] says so.
#include <Event/Event.h>
#include <Event/Eventiterator.h>
#include <Event/fileEventiterator.h>
#include <Event/packet.h>
#include <Event/EventTypes.h>
#include <Event/oncsSubConstants.h>
#include <ROC_pool/h2gcroc3_10G_pool.h>
#include <TFile.h>
#include <TTree.h>
#include <TNamed.h>
#include <TString.h>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <cstdint>
#include <deque>
#include <map>
#include <vector>
#include <string>
#include <algorithm>

static uint32_t be32(const uint8_t *p) { return (uint32_t(p[0]) << 24) | (uint32_t(p[1]) << 16) | (uint32_t(p[2]) << 8) | p[3]; }
static uint64_t be64(const uint8_t *p) { return (uint64_t(be32(p)) << 32) | be32(p + 4); }
static uint32_t le32(const uint8_t *p) { return (uint32_t(p[3]) << 24) | (uint32_t(p[2]) << 16) | (uint32_t(p[1]) << 8) | p[0]; }
// CRC-32, polynomial 0x04C11DB7, init 0, no reflection, no final xor, over w8..w46 (datasheet 1.3.1)
static uint32_t crc32_frame(const uint8_t *d, int n)
{
  uint32_t crc = 0;
  for (int i = 0; i < n; i++) { crc ^= uint32_t(d[i]) << 24; for (int b = 0; b < 8; b++) crc = (crc & 0x80000000u) ? (crc << 1) ^ 0x04C11DB7u : crc << 1; }
  return crc;
}

static const int NCH = 288, NS = 20, NL = 8;
static const int SFIL_S[4][4] = {{2, 6, 10, 14}, {21, 25, 29, 33}, {63, 59, 52, 48}, {43, 37, 0, 4}};          // [layer][cell] South end, ASIC0
static const int SFIL_N[4][4] = {{146, 150, 154, 158}, {165, 169, 173, 177}, {207, 203, 196, 192}, {187, 186, 144, 148}};   // North end, ASIC2
static const int PADDLE[4] = {84, 86, 81, 97}; static const char *PADDLE_NAME[4] = {"P3", "P4", "P5", "P6"};
static const char *PADDLE_CAEN[4] = {"CAEN ch 12", "CAEN ch 13", "CAEN ch 14", "CAEN ch 15"};
// KCU1: SFIL5 South, 16-SiPM arrays on ASIC2 (gaps 152, 161, 170, 179, 188, 197, 206, 215 not connected)
static const int SFIL5[4][16] = {{144, 145, 146, 147, 148, 149, 150, 151, 153, 154, 155, 156, 157, 158, 159, 160},
                                 {162, 163, 164, 165, 166, 167, 168, 169, 171, 172, 173, 174, 175, 176, 177, 178},
                                 {180, 181, 182, 183, 184, 185, 186, 187, 189, 190, 191, 192, 193, 194, 195, 196},
                                 {198, 199, 200, 201, 202, 203, 204, 205, 207, 208, 209, 210, 211, 212, 213, 214}};
static const char *SFIL5_NAME[4] = {"A", "B", "C", "D"};
static const int SFIL5N_KCU0[4] = {153, 156, 163, 167};   // run >= 3085: SFIL5 North cells 1-4 on KCU0 ASIC2
static bool sfil5n_on_kcu0 = false;
static const char *SFIL5_NOTE[4] = {"bottom, mostly dead (154, 156-160 fire); May-board analogy SiPM1", "second from bottom, live; analogy SiPM2", "top, live but weak; analogy SiPM4", "third from bottom, live; analogy SiPM3"};
static int board = 0;   // 0 = KCU0 (packet 12001), 1 = KCU1 (12002)
static const int NUNUSED0 = 180, NUNUSED1 = 80;

struct Line
{
  uint8_t b[192];
  uint64_t ts; int link; uint32_t udp; int evt; int seg; bool crc_ok; uint8_t ham;
};
struct PktMeta { int seg, evt; long long evt_time, rt, rt_mono, rt_diff; };

struct Out   // one set of event trees (main or head)
{
  TFile *f = nullptr; TTree *header = nullptr, *sfil[4] = {nullptr, nullptr, nullptr, nullptr}, *paddles = nullptr, *bbcal = nullptr, *unused = nullptr, *sfil5n = nullptr;
  long long n = 0;
  void fill() { header->Fill(); for (int L = 0; L < 4; L++) if (sfil[L]) sfil[L]->Fill(); if (paddles) paddles->Fill(); if (bbcal) bbcal->Fill(); if (sfil5n) sfil5n->Fill(); unused->Fill(); n++; }
  void write() { f->cd(); header->Write(); for (int L = 0; L < 4; L++) if (sfil[L]) sfil[L]->Write(); if (paddles) paddles->Write(); if (bbcal) bbcal->Write(); if (sfil5n) sfil5n->Write(); unused->Write(); }
};

// ---- the per-trigger record (branch addresses)
struct Rec
{
  Int_t sfil_ch[4][8], sfil5n_ch[4], paddle_ch[4], bbcal_ch[72], unused_ch[180], sfil5_ch[4][16], north_ch[144];   // constants, before the per-trigger block that emit() clears
  Int_t seg, run, nsamp, nlines, evt_first, evt_last, seg_last; Long64_t idx, ts0, evt_time, rt_first, rt_last, rt_mono_first, dts_prev;
  Bool_t complete, backward, nsamp_overflow; UInt_t udp_first, udp_last;
  Long64_t ts[NS]; UChar_t nlink[NS], linkmask[NS];
  Bool_t present[NS][NL]; Char_t dts[NS][NL];
  UInt_t trigger_in[NS][NL], trigger_out[NS][NL], event_counter[NS][NL], w0[NS][NL], w6[NS][NL], w7[NS][NL], roc_hdr[NS][NL], cm_word[NS][NL], calib_word[NS][NL], crc[NS][NL], udp_seq[NS][NL], pkt_evt[NS][NL];
  UShort_t bx[NS][NL], cm0[NS][NL], cm1[NS][NL], calib_adc[NS][NL], calib_toa[NS][NL], calib_tot[NS][NL];
  UChar_t ec[NS][NL], orbit[NS][NL], ham[NS][NL], hdr_nibbles[NS][NL], calib_tc[NS][NL], calib_tp[NS][NL], seg_line[NS][NL]; Bool_t crc_ok[NS][NL];
  Short_t adc[NCH][NS], toa[NCH][NS], tot[NCH][NS]; Char_t tc[NCH][NS], tp[NCH][NS];
  Short_t uadc[180][NS], utoa[180][NS], utot[180][NS]; Char_t utc[180][NS], utp[180][NS];
  Short_t s5adc[4][16][NS], s5toa[4][16][NS], s5tot[4][16][NS]; Char_t s5tc[4][16][NS], s5tp[4][16][NS];
};
static Rec R; static int unused_list[180]; static int nunused = 180;

static void book(Out &o, const char *fname, const char *title)
{
  o.f = new TFile(fname, "RECREATE", title, ROOT::RCompressionSetting::EDefaults::kUseGeneralPurpose);
  o.f->SetCompressionAlgorithm(ROOT::RCompressionSetting::EAlgorithm::kZSTD); o.f->SetCompressionLevel(5);
  TTree *h = o.header = new TTree("header", "per trigger: stream, RCDAQ, FPGA wrapper and H2GCROC3 frame-header information; [sample][link], link = 2*ASIC + half");
  h->Branch("seg", &R.seg, "seg/I"); h->Branch("run", &R.run, "run/I"); h->Branch("idx", &R.idx, "idx/L");
  h->Branch("ts0", &R.ts0, "ts0/L"); h->Branch("nsamp", &R.nsamp, "nsamp/I"); h->Branch("nlines", &R.nlines, "nlines/I"); h->Branch("complete", &R.complete, "complete/O");
  h->Branch("nsamp_overflow", &R.nsamp_overflow, "nsamp_overflow/O");
  h->Branch("dts_prev", &R.dts_prev, "dts_prev/L"); h->Branch("backward", &R.backward, "backward/O");
  h->Branch("evt_first", &R.evt_first, "evt_first/I"); h->Branch("evt_last", &R.evt_last, "evt_last/I"); h->Branch("seg_last", &R.seg_last, "seg_last/I");
  h->Branch("evt_time", &R.evt_time, "evt_time/L"); h->Branch("rt_first", &R.rt_first, "rt_first/L"); h->Branch("rt_last", &R.rt_last, "rt_last/L"); h->Branch("rt_mono_first", &R.rt_mono_first, "rt_mono_first/L");
  h->Branch("udp_first", &R.udp_first, "udp_first/i"); h->Branch("udp_last", &R.udp_last, "udp_last/i");
  h->Branch("ts", R.ts, "ts[20]/L"); h->Branch("nlink", R.nlink, "nlink[20]/b"); h->Branch("linkmask", R.linkmask, "linkmask[20]/b");
  h->Branch("present", R.present, "present[20][8]/O"); h->Branch("dts", R.dts, "dts[20][8]/B");
  h->Branch("trigger_in", R.trigger_in, "trigger_in[20][8]/i"); h->Branch("trigger_out", R.trigger_out, "trigger_out[20][8]/i"); h->Branch("event_counter", R.event_counter, "event_counter[20][8]/i");
  h->Branch("w0", R.w0, "w0[20][8]/i"); h->Branch("w6", R.w6, "w6[20][8]/i"); h->Branch("w7", R.w7, "w7[20][8]/i");
  h->Branch("roc_hdr", R.roc_hdr, "roc_hdr[20][8]/i"); h->Branch("bx", R.bx, "bx[20][8]/s"); h->Branch("ec", R.ec, "ec[20][8]/b"); h->Branch("orbit", R.orbit, "orbit[20][8]/b"); h->Branch("ham", R.ham, "ham[20][8]/b"); h->Branch("hdr_nibbles", R.hdr_nibbles, "hdr_nibbles[20][8]/b");
  h->Branch("cm_word", R.cm_word, "cm_word[20][8]/i"); h->Branch("cm0", R.cm0, "cm0[20][8]/s"); h->Branch("cm1", R.cm1, "cm1[20][8]/s");
  h->Branch("calib_word", R.calib_word, "calib_word[20][8]/i"); h->Branch("calib_adc", R.calib_adc, "calib_adc[20][8]/s"); h->Branch("calib_toa", R.calib_toa, "calib_toa[20][8]/s"); h->Branch("calib_tot", R.calib_tot, "calib_tot[20][8]/s"); h->Branch("calib_tc", R.calib_tc, "calib_tc[20][8]/b"); h->Branch("calib_tp", R.calib_tp, "calib_tp[20][8]/b");
  h->Branch("crc", R.crc, "crc[20][8]/i"); h->Branch("crc_ok", R.crc_ok, "crc_ok[20][8]/O");
  h->Branch("udp_seq", R.udp_seq, "udp_seq[20][8]/i"); h->Branch("pkt_evt", R.pkt_evt, "pkt_evt[20][8]/i"); h->Branch("seg_line", R.seg_line, "seg_line[20][8]/b");
  if (board == 1)
    {
      TTree *t = o.sfil[0] = new TTree("sfil5", "SFIL5 South, individual SiPMs: arrays A-D of 16 channels on KCU1 ASIC2 (vertical order from data A < B < D < C); -1 = line missing");
      for (int a = 0; a < 4; a++)
        {
          const char *nm = SFIL5_NAME[a]; for (int i = 0; i < 16; i++) R.sfil5_ch[a][i] = SFIL5[a][i];
          t->Branch(Form("%s_ch", nm), R.sfil5_ch[a], Form("%s_ch[16]/I", nm));
          t->Branch(Form("%s_adc", nm), R.s5adc[a], Form("%s_adc[16][20]/S", nm)); t->Branch(Form("%s_toa", nm), R.s5toa[a], Form("%s_toa[16][20]/S", nm)); t->Branch(Form("%s_tot", nm), R.s5tot[a], Form("%s_tot[16][20]/S", nm));
          t->Branch(Form("%s_tc", nm), (void *) R.s5tc[a], Form("%s_tc[16][20]/B", nm)); t->Branch(Form("%s_tp", nm), (void *) R.s5tp[a], Form("%s_tp[16][20]/B", nm));
        }
      TTree *n = o.sfil[1] = new TTree("sfil5north", "SFIL5 North board on KCU1 ASIC0/ASIC1 (ch 0-143; dead in run 3062: 0-35 and 72-107 railed, the rest quiet); -1 = line missing");
      for (int i = 0; i < 144; i++) R.north_ch[i] = i;
      n->Branch("ch", R.north_ch, "ch[144]/I"); n->Branch("adc", R.adc[0], "adc[144][20]/S"); n->Branch("toa", R.toa[0], "toa[144][20]/S"); n->Branch("tot", R.tot[0], "tot[144][20]/S"); n->Branch("tc", (void *) R.tc[0], "tc[144][20]/B"); n->Branch("tp", (void *) R.tp[0], "tp[144][20]/B");
      TTree *u = o.unused = new TTree("unused", "the 80 KCU1 channels without a detector (ASIC2 gaps and ASIC3); ch[i] gives the channel number; -1 = line missing");
      u->Branch("ch", R.unused_ch, "ch[80]/I"); u->Branch("adc", R.uadc, "adc[80][20]/S"); u->Branch("toa", R.utoa, "toa[80][20]/S"); u->Branch("tot", R.utot, "tot[80][20]/S"); u->Branch("tc", (void *) R.utc, "tc[80][20]/B"); u->Branch("tp", (void *) R.utp, "tp[80][20]/B");
      return;
    }
  for (int L = 0; L < 4; L++)
    {
      TTree *t = o.sfil[L] = new TTree(Form("sfil%d", L + 1), Form("SFIL%d: South end (KCU0 ASIC0) and North end (ASIC2) of cells 1-4, 20 samples of 25 ns; -1 = line missing", L + 1));
      for (int e = 0; e < 2; e++) for (int c = 0; c < 4; c++)
        {
          const int ch = e == 0 ? SFIL_S[L][c] : SFIL_N[L][c]; const char *nm = Form("%cc%d", e == 0 ? 'S' : 'N', c + 1);
          R.sfil_ch[L][4 * e + c] = ch;
          t->Branch(Form("%s_ch", nm), &R.sfil_ch[L][4 * e + c], Form("%s_ch/I", nm));
          t->Branch(Form("%s_adc", nm), R.adc[ch], Form("%s_adc[20]/S", nm)); t->Branch(Form("%s_toa", nm), R.toa[ch], Form("%s_toa[20]/S", nm)); t->Branch(Form("%s_tot", nm), R.tot[ch], Form("%s_tot[20]/S", nm));
          t->Branch(Form("%s_tc", nm), (void *) R.tc[ch], Form("%s_tc[20]/B", nm)); t->Branch(Form("%s_tp", nm), (void *) R.tp[ch], Form("%s_tp[20]/B", nm));
        }
    }
  TTree *p = o.paddles = new TTree("paddles", "trigger paddles P3-P6 on KCU0 ASIC1 (negative pulses); -1 = line missing");
  for (int k = 0; k < 4; k++)
    {
      const int ch = PADDLE[k]; const char *nm = PADDLE_NAME[k]; R.paddle_ch[k] = ch;
      p->Branch(Form("%s_ch", nm), &R.paddle_ch[k], Form("%s_ch/I", nm));
      p->Branch(Form("%s_adc", nm), R.adc[ch], Form("%s_adc[20]/S", nm)); p->Branch(Form("%s_toa", nm), R.toa[ch], Form("%s_toa[20]/S", nm)); p->Branch(Form("%s_tot", nm), R.tot[ch], Form("%s_tot[20]/S", nm));
      p->Branch(Form("%s_tc", nm), (void *) R.tc[ch], Form("%s_tc[20]/B", nm)); p->Branch(Form("%s_tp", nm), (void *) R.tp[ch], Form("%s_tp[20]/B", nm));
    }
  if (sfil5n_on_kcu0)
    {
      TTree *t = o.sfil5n = new TTree("sfil5n", "SFIL5 North summed boards on KCU0 ASIC2 (runs >= 3085): cells 1-4 on ch 153, 156, 163, 167; -1 = line missing");
      for (int c = 0; c < 4; c++)
        {
          const int ch = SFIL5N_KCU0[c]; const char *nm = Form("Nc%d", c + 1); R.sfil5n_ch[c] = ch;
          t->Branch(Form("%s_ch", nm), &R.sfil5n_ch[c], Form("%s_ch/I", nm));
          t->Branch(Form("%s_adc", nm), R.adc[ch], Form("%s_adc[20]/S", nm)); t->Branch(Form("%s_toa", nm), R.toa[ch], Form("%s_toa[20]/S", nm)); t->Branch(Form("%s_tot", nm), R.tot[ch], Form("%s_tot[20]/S", nm));
          t->Branch(Form("%s_tc", nm), (void *) R.tc[ch], Form("%s_tc[20]/B", nm)); t->Branch(Form("%s_tp", nm), (void *) R.tp[ch], Form("%s_tp[20]/B", nm));
        }
    }
  TTree *b = o.bbcal = new TTree("bbcal", "Baby BCAL North on KCU0 ASIC3: index i = channel 216+i; -1 = line missing");
  for (int i = 0; i < 72; i++) R.bbcal_ch[i] = 216 + i;
  b->Branch("ch", R.bbcal_ch, "ch[72]/I"); b->Branch("adc", R.adc[216], "adc[72][20]/S"); b->Branch("toa", R.toa[216], "toa[72][20]/S"); b->Branch("tot", R.tot[216], "tot[72][20]/S"); b->Branch("tc", (void *) R.tc[216], "tc[72][20]/B"); b->Branch("tp", (void *) R.tp[216], "tp[72][20]/B");
  TTree *u = o.unused = new TTree("unused", Form("the %d KCU0 channels without a detector (ASIC0/1/2); ch[i] gives the channel number; -1 = line missing", nunused));
  u->Branch("ch", R.unused_ch, Form("ch[%d]/I", nunused)); u->Branch("adc", R.uadc, Form("adc[%d][20]/S", nunused)); u->Branch("toa", R.utoa, Form("toa[%d][20]/S", nunused)); u->Branch("tot", R.utot, Form("tot[%d][20]/S", nunused)); u->Branch("tc", (void *) R.utc, Form("tc[%d][20]/B", nunused)); u->Branch("tp", (void *) R.utp, Form("tp[%d][20]/B", nunused));
}

int main(int argc, char **argv)
{
  if (argc < 4) { fprintf(stderr, "usage: kcu_to_root out.root segment.evt seg [next_segment.evt|-] [n_head_packets=6] [packetid=12001] [max_packets=0] [run=0]\n"); return 1; }
  const char *outname = argv[1]; const char *fname = argv[2]; const int seg = atoi(argv[3]);
  const char *nextname = argc > 4 && strcmp(argv[4], "-") ? argv[4] : nullptr;
  const int nhead = argc > 5 ? atoi(argv[5]) : 6; const int pid = argc > 6 ? atoi(argv[6]) : 12001; const long maxpk = argc > 7 ? atol(argv[7]) : 0;
  const int run_hookup = argc > 8 ? atoi(argv[8]) : 0;

  // channel bookkeeping
  board = (pid == 12002) ? 1 : 0;
  sfil5n_on_kcu0 = board == 0 && run_hookup >= 3085;
  bool used[NCH] = {false};
  if (board == 0)
    {
      for (int L = 0; L < 4; L++) for (int c = 0; c < 4; c++) { used[SFIL_S[L][c]] = used[SFIL_N[L][c]] = true; }
      for (int k = 0; k < 4; k++) used[PADDLE[k]] = true;
      for (int c = 216; c < NCH; c++) used[c] = true;
      if (sfil5n_on_kcu0) for (int c = 0; c < 4; c++) used[SFIL5N_KCU0[c]] = true;
    }
  else
    {
      for (int a = 0; a < 4; a++) for (int i = 0; i < 16; i++) used[SFIL5[a][i]] = true;
      for (int c = 0; c < 144; c++) used[c] = true;
    }
  nunused = board == 0 ? NUNUSED0 - (sfil5n_on_kcu0 ? 4 : 0) : NUNUSED1;
  { int n = 0; for (int c = 0; c < NCH; c++) if (!used[c]) unused_list[n++] = c; if (n != nunused) { fprintf(stderr, "unused channel count %d != %d\n", n, nunused); return 1; } for (int i = 0; i < nunused; i++) R.unused_ch[i] = unused_list[i]; }

  Out main_out, head_out;
  book(main_out, outname, Form("KCU packet %d, run segment %d: triggers in stream order", pid, seg));
  std::string headname = outname; if (headname.size() > 5 && headname.substr(headname.size() - 5) == ".root") headname = headname.substr(0, headname.size() - 5); headname += "_head.root";
  book(head_out, headname.c_str(), Form("KCU packet %d, segment %d: triggers with a line in the segment's first packet", pid, seg));

  // packets tree (main file)
  main_out.f->cd();
  TTree *tp = new TTree("packets", "one entry per RCDAQ event of this segment, in file order");
  Int_t p_seg, p_evt, p_type, p_run, p_evtlen, p_pktlen, p_datalen, p_ndg, p_jumps, p_nlines, p_nomarker, p_badlink, p_dup, p_crcbad, p_ham, p_tc, p_tp, p_tailnz, p_nwf; Long64_t p_time, p_rt, p_rtmono, p_rtdiff, p_lost, p_tsmin, p_tsmax; UInt_t p_udpfirst, p_udplast, p_const; Bool_t p_haskcu, p_skipped;
  tp->Branch("seg", &p_seg, "seg/I"); tp->Branch("evt", &p_evt, "evt/I"); tp->Branch("evt_type", &p_type, "evt_type/I"); tp->Branch("run", &p_run, "run/I"); tp->Branch("evt_time", &p_time, "evt_time/L"); tp->Branch("evt_len", &p_evtlen, "evt_len/I");
  tp->Branch("has_kcu", &p_haskcu, "has_kcu/O"); tp->Branch("skipped", &p_skipped, "skipped/O"); tp->Branch("pkt_len", &p_pktlen, "pkt_len/I"); tp->Branch("pkt_datalen", &p_datalen, "pkt_datalen/I");
  tp->Branch("rt", &p_rt, "rt/L"); tp->Branch("rt_mono", &p_rtmono, "rt_mono/L"); tp->Branch("rt_diff", &p_rtdiff, "rt_diff/L");
  tp->Branch("ndatagram", &p_ndg, "ndatagram/I"); tp->Branch("udp_first", &p_udpfirst, "udp_first/i"); tp->Branch("udp_last", &p_udplast, "udp_last/i"); tp->Branch("udp_jumps", &p_jumps, "udp_jumps/I"); tp->Branch("udp_lost", &p_lost, "udp_lost/L");
  tp->Branch("dg_const", &p_const, "dg_const/i"); tp->Branch("dg_tail_nonzero", &p_tailnz, "dg_tail_nonzero/I");
  tp->Branch("nlines", &p_nlines, "nlines/I"); tp->Branch("nomarker", &p_nomarker, "nomarker/I"); tp->Branch("badlink", &p_badlink, "badlink/I"); tp->Branch("dup_lines", &p_dup, "dup_lines/I");
  tp->Branch("crc_bad", &p_crcbad, "crc_bad/I"); tp->Branch("ham_flagged", &p_ham, "ham_flagged/I"); tp->Branch("tc_set", &p_tc, "tc_set/I"); tp->Branch("tp_set", &p_tp, "tp_set/I");
  tp->Branch("ts_min", &p_tsmin, "ts_min/L"); tp->Branch("ts_max", &p_tsmax, "ts_max/L"); tp->Branch("nwf_emitted", &p_nwf, "nwf_emitted/I");

  // map tree (segment 0 only)
  if (seg == 0)
    {
      TTree *tm = new TTree("map", Form("one entry per KCU%d channel: location and where to find it", board));
      Int_t m_ch, m_asic, m_half, m_link, m_k, m_layer, m_cell; Char_t m_sub[16], m_end[2], m_label[48], m_tree[16], m_branch[32];
      tm->Branch("ch", &m_ch, "ch/I"); tm->Branch("asic", &m_asic, "asic/I"); tm->Branch("half", &m_half, "half/I"); tm->Branch("link", &m_link, "link/I"); tm->Branch("k", &m_k, "k/I");
      tm->Branch("subsystem", m_sub, "subsystem/C"); tm->Branch("end", m_end, "end/C"); tm->Branch("layer", &m_layer, "layer/I"); tm->Branch("cell", &m_cell, "cell/I"); tm->Branch("label", m_label, "label/C"); tm->Branch("tree", m_tree, "tree/C"); tm->Branch("branch", m_branch, "branch/C");
      for (int ch = 0; ch < NCH; ch++)
        {
          m_ch = ch; m_asic = ch / 72; m_half = (ch % 72) / 36; m_link = 2 * m_asic + m_half; m_k = ch % 36; m_layer = m_cell = 0; strcpy(m_end, ""); strcpy(m_sub, "unused"); strcpy(m_tree, "unused"); strcpy(m_label, "no detector"); int ui = -1;
          for (int i = 0; i < nunused; i++) if (unused_list[i] == ch) ui = i;
          if (ui >= 0) sprintf(m_branch, "adc[%d]", ui);
          if (board == 1)
            {
              for (int a = 0; a < 4; a++) for (int i = 0; i < 16; i++) if (SFIL5[a][i] == ch) { strcpy(m_sub, "SFIL5"); strcpy(m_tree, "sfil5"); strcpy(m_end, "S"); m_layer = 5; m_cell = a + 1; sprintf(m_branch, "%s_adc[%d]", SFIL5_NAME[a], i); sprintf(m_label, "SFIL5 South array %s SiPM %d: %s", SFIL5_NAME[a], i, SFIL5_NOTE[a]); }
              if (ch < 144) { strcpy(m_sub, "SFIL5N"); strcpy(m_tree, "sfil5north"); strcpy(m_end, "N"); m_layer = 5; sprintf(m_branch, "adc[%d]", ch); strcpy(m_label, (ch < 36 || (ch >= 72 && ch < 108)) ? "SFIL5 North board, railed at 1023 in run 3062" : "SFIL5 North board, quiet in run 3062"); }
              tm->Fill(); continue;
            }
          for (int L = 0; L < 4; L++) for (int c = 0; c < 4; c++) for (int e = 0; e < 2; e++) if ((e == 0 ? SFIL_S[L][c] : SFIL_N[L][c]) == ch)
            { sprintf(m_sub, "SFIL%d", L + 1); sprintf(m_tree, "sfil%d", L + 1); strcpy(m_end, e == 0 ? "S" : "N"); m_layer = L + 1; m_cell = c + 1; sprintf(m_branch, "%cc%d_adc", e == 0 ? 'S' : 'N', c + 1);
              sprintf(m_label, "SFIL%d cell %d %s end%s", L + 1, c + 1, e == 0 ? "South" : "North", c == 2 ? " (beam row, SiPM3)" : c == 3 ? " (SiPM4)" : (L == 0 && e == 1) ? (c == 0 ? " (v1 board, SiPM1)" : " (v1 board, SiPM2)") : " (SiPM1+2 sum, v2 board)"); }
          for (int k = 0; k < 4; k++) if (PADDLE[k] == ch) { strcpy(m_sub, "paddles"); strcpy(m_tree, "paddles"); sprintf(m_branch, "%s_adc", PADDLE_NAME[k]); sprintf(m_label, "trigger paddle %s (R580 PMT, negative; also %s)", PADDLE_NAME[k], PADDLE_CAEN[k]); }
          for (int c = 0; c < 4; c++) if (sfil5n_on_kcu0 && SFIL5N_KCU0[c] == ch) { strcpy(m_sub, "SFIL5"); strcpy(m_tree, "sfil5n"); strcpy(m_end, "N"); m_layer = 5; m_cell = c + 1; sprintf(m_branch, "Nc%d_adc", c + 1); sprintf(m_label, "SFIL5 cell %d North end (summed board)", c + 1); }
          if (ch >= 216) { strcpy(m_sub, "bbcal"); strcpy(m_tree, "bbcal"); sprintf(m_branch, "adc[%d]", ch - 216); strcpy(m_label, (ch == 216 || ch >= 252) ? "Baby BCAL North, railed at 1023 in run 3062" : "Baby BCAL North (BBS on J1-J40), no signal in run 3062"); }
          tm->Fill();
        }
      tm->Write();
    }

  // raw line store: the last KEEP packets, indexed by timestamp
  const size_t KEEP = 12;
  std::deque<std::vector<Line>> store; std::multimap<uint64_t, const Line *> index; std::map<std::pair<int, int>, PktMeta> meta;
  auto prune = [&]()
    {
      while (store.size() > KEEP)
        {
          for (const Line &l : store.front()) { auto r = index.equal_range(l.ts); for (auto it = r.first; it != r.second; ++it) if (it->second == &l) { index.erase(it); break; } }
          const Line &l0 = store.front().front(); meta.erase(std::make_pair(l0.seg, l0.evt));
          store.pop_front();
        }
    };
  auto find_line = [&](uint64_t ts, int link) -> const Line *
    {
      for (auto it = index.lower_bound(ts > 5 ? ts - 5 : 0); it != index.end() && it->first <= ts + 5; ++it) if (it->second->link == link) return it->second;
      return nullptr;
    };

  h2gcroc3_10G_pool *pool = nullptr;
  long long nkept = 0, nhead_kept = 0, ndropped_next = 0, nwf_total = 0, adc_mismatch = 0, nbackward = 0, nolines = 0;
  int first_evt = -1; long long prev_ts0 = -1;
  int cur_seg = seg;

  // ---- emit the waveforms on offer
  auto emit = [&](int &nwf_counter)
    {
      const int n = pool->iValue(0, "NR_WF");
      for (int w = 0; w < n; w++)
        {
          nwf_total++; nwf_counter++;
          memset(&R.ts, 0, sizeof(Rec) - offsetof(Rec, ts));   // everything after the scalars
          const int ns_real = pool->iValue(w, "SAMPLESIZE"); const int ns = std::min(ns_real, NS);
          R.seg = seg; R.run = 0; R.nsamp = ns_real; R.nsamp_overflow = ns_real > NS; R.nlines = 0;
          for (int c = 0; c < NCH; c++) for (int s = 0; s < NS; s++) { R.adc[c][s] = R.toa[c][s] = R.tot[c][s] = -1; R.tc[c][s] = R.tp[c][s] = -1; }
          int owner = -1; bool inhead = false; long long ev_first = -1, ev_last = -1; int sg_first = -1, sg_last = -1; uint32_t u_first = 0, u_last = 0; bool have_first = false;
          for (int s = 0; s < ns; s++)
            {
              const uint64_t ts = (uint64_t) pool->lValue(w, s, 0); R.ts[s] = (Long64_t) ts;
              for (int l = 0; l < NL; l++)
                {
                  const Line *L = find_line(ts, l);
                  if (!L) continue;
                  const uint8_t *d = L->b;
                  R.present[s][l] = true; R.nlink[s]++; R.linkmask[s] |= (1 << l); R.nlines++;
                  R.dts[s][l] = (Char_t) ((long long) L->ts - (long long) ts);
                  R.w0[s][l] = be32(d); R.trigger_in[s][l] = be32(d + 4); R.trigger_out[s][l] = be32(d + 8); R.event_counter[s][l] = be32(d + 12); R.w6[s][l] = be32(d + 24); R.w7[s][l] = be32(d + 28);
                  const uint32_t hdr = be32(d + 32); R.roc_hdr[s][l] = hdr; R.bx[s][l] = (hdr >> 16) & 0xfff; R.ec[s][l] = (hdr >> 10) & 0x3f; R.orbit[s][l] = (hdr >> 7) & 7; R.ham[s][l] = (hdr >> 4) & 7; R.hdr_nibbles[s][l] = ((hdr >> 28) << 4) | (hdr & 0xf);
                  const uint32_t cm = be32(d + 36); R.cm_word[s][l] = cm; R.cm0[s][l] = (cm >> 10) & 0x3ff; R.cm1[s][l] = cm & 0x3ff;
                  const uint32_t cw = be32(d + 112); R.calib_word[s][l] = cw; R.calib_adc[s][l] = (cw >> 20) & 0x3ff; R.calib_tot[s][l] = (cw >> 10) & 0x3ff; R.calib_toa[s][l] = cw & 0x3ff; R.calib_tc[s][l] = cw >> 31; R.calib_tp[s][l] = (cw >> 30) & 1;
                  R.crc[s][l] = be32(d + 188); R.crc_ok[s][l] = L->crc_ok; R.udp_seq[s][l] = L->udp; R.pkt_evt[s][l] = L->evt; R.seg_line[s][l] = L->seg;
                  for (int k = 0; k < 36; k++)
                    {
                      const uint32_t wd = be32(d + (k < 18 ? 40 + 4 * k : 116 + 4 * (k - 18))); const int ch = 36 * l + k;
                      R.adc[ch][s] = (wd >> 20) & 0x3ff; R.tot[ch][s] = (wd >> 10) & 0x3ff; R.toa[ch][s] = wd & 0x3ff; R.tc[ch][s] = wd >> 31; R.tp[ch][s] = (wd >> 30) & 1;
                    }
                  if (s == 0 && owner < 0) owner = L->seg;
                  if (L->seg == seg && L->evt == first_evt) inhead = true;
                  if (!have_first) { ev_first = L->evt; sg_first = L->seg; u_first = L->udp; have_first = true; }
                  ev_last = L->evt; sg_last = L->seg; u_last = L->udp;
                }
            }
          // cross-check with the decoder's own ADC values (both come from the same bytes)
          for (int s = 0; s < ns; s++) for (int c = 0; c < NCH; c += 7) if (R.present[s][2 * (c / 72) + (c % 72) / 36] && R.adc[c][s] != pool->iValue(w, s, c)) adc_mismatch++;
          if (owner < 0) { nolines++; owner = cur_seg; }
          if (owner != seg) { ndropped_next++; continue; }   // sample 0 lies in the next segment: that job writes it
          R.ts0 = R.ts[0]; R.complete = (ns_real == NS && R.nlines == NS * NL);
          R.evt_first = (Int_t) ev_first; R.evt_last = (Int_t) ev_last; R.seg_last = sg_last; R.udp_first = u_first; R.udp_last = u_last;
          const PktMeta *mf = nullptr, *ml = nullptr; auto itf = meta.find(std::make_pair(sg_first, (int) ev_first)); if (itf != meta.end()) mf = &itf->second; auto itl = meta.find(std::make_pair(sg_last, (int) ev_last)); if (itl != meta.end()) ml = &itl->second;
          R.evt_time = mf ? mf->evt_time : 0; R.rt_first = mf ? mf->rt : 0; R.rt_mono_first = mf ? mf->rt_mono : 0; R.rt_last = ml ? ml->rt : 0;
          R.dts_prev = prev_ts0 < 0 ? 0 : R.ts0 - prev_ts0; R.backward = prev_ts0 >= 0 && R.ts0 < prev_ts0; if (R.backward) nbackward++; prev_ts0 = R.ts0;
          for (int i = 0; i < nunused; i++) { const int c = unused_list[i]; memcpy(R.uadc[i], R.adc[c], sizeof(R.adc[c])); memcpy(R.utoa[i], R.toa[c], sizeof(R.toa[c])); memcpy(R.utot[i], R.tot[c], sizeof(R.tot[c])); memcpy(R.utc[i], R.tc[c], sizeof(R.tc[c])); memcpy(R.utp[i], R.tp[c], sizeof(R.tp[c])); }
          if (board == 1) for (int a = 0; a < 4; a++) for (int i = 0; i < 16; i++) { const int c = SFIL5[a][i]; memcpy(R.s5adc[a][i], R.adc[c], sizeof(R.adc[c])); memcpy(R.s5toa[a][i], R.toa[c], sizeof(R.toa[c])); memcpy(R.s5tot[a][i], R.tot[c], sizeof(R.tot[c])); memcpy(R.s5tc[a][i], R.tc[c], sizeof(R.tc[c])); memcpy(R.s5tp[a][i], R.tp[c], sizeof(R.tp[c])); }
          if (inhead) { R.idx = head_out.n; head_out.fill(); nhead_kept++; } else { R.idx = main_out.n; main_out.fill(); nkept++; }
        }
    };

  // ---- one RCDAQ event: raw parse (store + packets tree), feed the pool, emit
  auto process_event = [&](Event *e, int this_seg, bool record_packet)
    {
      p_seg = this_seg; p_evt = e->getEvtSequence(); p_type = e->getEvtType(); p_run = e->getRunNumber(); p_time = (Long64_t) e->getTime(); p_evtlen = e->getEvtLength();
      p_haskcu = false; p_skipped = false; p_pktlen = p_datalen = p_ndg = p_jumps = p_nlines = p_nomarker = p_badlink = p_dup = p_crcbad = p_ham = p_tc = p_tp = p_tailnz = p_nwf = 0; p_lost = 0; p_tsmin = p_tsmax = 0; p_udpfirst = p_udplast = p_const = 0; p_rt = p_rtmono = p_rtdiff = 0;
      for (int id = 1010; id <= 1013; id++) if (Packet *q = e->getPacket(id)) { p_rt = q->lValue(0, "TOTAL"); p_rtmono = q->lValue(2, "TOTAL"); p_rtdiff = q->lValue(0, "TIMEDIFF"); delete q; break; }
      if (R.run == 0) R.run = p_run;
      Packet *p = (p_type == DATAEVENT) ? e->getPacket(pid) : nullptr;
      if (p)
        {
          p_haskcu = true; p_pktlen = p->getLength(); p_datalen = p->getDataLength();
          if (first_evt < 0 && this_seg == seg) first_evt = p_evt;
          meta[std::make_pair(this_seg, p_evt)] = PktMeta{this_seg, p_evt, p_time, p_rt, p_rtmono, p_rtdiff};
          std::vector<int> raw(p->getLength()); p->copyMe(raw.data(), raw.size());
          const uint8_t *c = reinterpret_cast<const uint8_t *>(raw.data()) + 16; const long nbytes = 4L * p->getDataLength();
          std::vector<Line> lines; lines.reserve(1400);
          long pos = 0; uint32_t prevseq = 0; bool haveseq = false; uint64_t tmin = ~0ULL, tmax = 0;
          while (pos + 22 <= nbytes)
            {
              const long plen = long(le32(c + pos + 4)) + 8;
              if (le32(c + pos) != 0xaabbccdd || plen < 22 || pos + plen > nbytes) break;
              const uint32_t sq = be32(c + pos + 8);
              if (p_ndg == 0) { p_udpfirst = sq; p_const = le32(c + pos + 12); }
              p_udplast = sq;
              if (haveseq && sq != prevseq + 1) { p_jumps++; p_lost += (long long) sq - prevseq - 1; }
              prevseq = sq; haveseq = true;
              bool tailnz = false; for (int i = 16; i < 22; i++) if (c[pos + i]) tailnz = true; if (tailnz) p_tailnz++;
              p_ndg++;
              for (long l = pos + 22; l + 192 <= pos + plen; l += 192)
                {
                  const uint8_t *d = c + l; p_nlines++;
                  if (d[0] != 0xaa || d[1] != 0x5a) { p_nomarker++; continue; }
                  const int link = 2 * (d[2] & 0xf) + (d[3] - 0x24);
                  if (link < 0 || link > 7) { p_badlink++; continue; }
                  Line L; memcpy(L.b, d, 192); L.ts = be64(d + 16); L.link = link; L.udp = sq; L.evt = p_evt; L.seg = this_seg;
                  L.crc_ok = crc32_frame(d + 32, 156) == be32(d + 188); L.ham = (be32(d + 32) >> 4) & 7;
                  if (!L.crc_ok) p_crcbad++;
                  if (L.ham) p_ham++;
                  for (int k = 0; k < 36; k++) { const uint32_t wd = be32(d + (k < 18 ? 40 + 4 * k : 116 + 4 * (k - 18))); if (wd >> 31) p_tc++; if ((wd >> 30) & 1) p_tp++; }
                  if (L.ts < tmin) tmin = L.ts;
                  if (L.ts > tmax) tmax = L.ts;
                  if (find_line(L.ts, link)) p_dup++;
                  lines.push_back(L);
                }
              pos += plen;
            }
          p_tsmin = (Long64_t) tmin; p_tsmax = (Long64_t) tmax;
          store.push_back(std::move(lines));
          for (const Line &L : store.back()) index.emplace(L.ts, &L);
          if (!pool) pool = new h2gcroc3_10G_pool();
          pool->addPacket(p); delete p;
          while (pool->depth_ok()) { emit(p_nwf); pool->next(); }
          prune();
        }
      if (record_packet) tp->Fill();
    };

  int status = 0;
  {
    fileEventiterator it(fname, status);
    if (status) { fprintf(stderr, "cannot open %s\n", fname); return 1; }
    long npk = 0;
    while (Event *e = it.getNextEvent())
      {
        process_event(e, seg, true); delete e;
        if (maxpk && ++npk >= maxpk) break;
      }
  }
  if (nextname)
    {
      cur_seg = seg + 1;
      fileEventiterator it(nextname, status);
      if (status) { fprintf(stderr, "cannot open %s\n", nextname); return 1; }
      int n = 0;
      while (Event *e = it.getNextEvent())
        {
          const bool data = e->getEvtType() == DATAEVENT;
          process_event(e, seg + 1, false); delete e;
          if (data && ++n >= nhead) break;
        }
    }
  if (pool) { int dummy = 0; pool->drain(); while (pool->depth_ok()) { emit(dummy); pool->next(); } }

  main_out.f->cd(); tp->Write();
  TString info = Form("kcu_to_root: %s (segment %d, packet %d), next segment head: %s (%d packets). Triggers kept %lld (+%lld in the head file), left to the next segment %lld, waveforms seen %lld; "
                      "pool: CRC-flagged lines %llu, hamming-flagged %llu, bad lines %llu, late lines %llu, resets %llu; ADC cross-check mismatches %lld; backward ts0 steps %lld; waveforms without any line found %lld. "
                      "TOT is the 10-bit compressed code (v<512: v, else 8*(v-512)+4 counts); TOA code = fine[2:0] x 25 ps + mid[7:3] x 0.2 ns + coarse[9:8] x 6.25 ns nominal; ts in 6.25 ns ticks; -1 = line missing.",
                      fname, seg, pid, nextname ? nextname : "-", nextname ? nhead : 0, nkept, nhead_kept, ndropped_next, nwf_total,
                      pool ? pool->crc_errors() : 0ULL, pool ? pool->hamming_errors() : 0ULL, pool ? pool->bad_lines() : 0ULL, pool ? pool->late_lines() : 0ULL, pool ? pool->resets() : 0ULL, adc_mismatch, nbackward, nolines);
  TNamed(Form("info_seg%d", seg), info.Data()).Write();
  main_out.write(); main_out.f->Close();
  head_out.f->cd(); TNamed(Form("info_seg%d_head", seg), info.Data()).Write(); head_out.write(); head_out.f->Close();
  printf("%s\n", info.Data());
  return 0;
}
