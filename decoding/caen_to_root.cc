// caen_to_root: every CAEN DT5742 trigger of one segment (packet 2001, one trigger per RCDAQ event)
// into ROOT trees, one entry per trigger in file order, every word of the packet decoded.
//
//   caen_to_root out.root segment.evt seg [max_events=0]
//
// Trees (aligned, one entry per trigger):
//   header   RCDAQ event, RTCLOCK, the four CAEN header words and both group headers / trigger times
//   bbcal    Baby BCAL South, CAEN ch 0-11: adc[12][256] (12-bit raw DRS4 samples, 1 GS/s, negative pulses)
//   paddles  P3 P4 P5 P6 = CAEN ch 12-15: P3_adc[256] ...
//   packets  one entry per RCDAQ event (also the non-data events)
// Words are cross-checked against the stock decoder (iValue(sample, ch)): see the info object.
#include <Event/Event.h>
#include <Event/Eventiterator.h>
#include <Event/fileEventiterator.h>
#include <Event/packet.h>
#include <Event/EventTypes.h>
#include <TFile.h>
#include <TTree.h>
#include <TNamed.h>
#include <TString.h>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <cstdint>
#include <vector>

int main(int argc, char **argv)
{
  if (argc < 4) { fprintf(stderr, "usage: caen_to_root out.root segment.evt seg [max_events]\n"); return 1; }
  const int seg = atoi(argv[3]); const long maxev = argc > 4 ? atol(argv[4]) : 0; const int pid = 2001;
  TFile *f = new TFile(argv[1], "RECREATE", Form("CAEN DT5742 packet 2001, run segment %d", seg)); f->SetCompressionAlgorithm(ROOT::RCompressionSetting::EAlgorithm::kZSTD); f->SetCompressionLevel(5);
  // header
  Int_t seg_b, run, evt, evt_type; Long64_t idx, evt_time, rt, rt_mono, rt_diff;
  UInt_t w0, w1, w2, w3, marker, evsize, board_id, lvds, gmask, w1_res, w2_res, evnr, evttime;
  UInt_t g_hdr[2], g_idxcell[2], g_freq[2], g_tr[2], g_nwords[2], g_res[2], g_ttword[2], g_tt[2], g_ttres[2]; Long64_t gtt;
  Bool_t ok_size, ok_marker;
  TTree *h = new TTree("header", "per trigger: RCDAQ event, RTCLOCK, CAEN DT5742 event header words w0-w3 and the two group headers / trigger-time words");
  h->Branch("seg", &seg_b, "seg/I"); h->Branch("run", &run, "run/I"); h->Branch("idx", &idx, "idx/L"); h->Branch("evt", &evt, "evt/I"); h->Branch("evt_type", &evt_type, "evt_type/I"); h->Branch("evt_time", &evt_time, "evt_time/L");
  h->Branch("rt", &rt, "rt/L"); h->Branch("rt_mono", &rt_mono, "rt_mono/L"); h->Branch("rt_diff", &rt_diff, "rt_diff/L");
  h->Branch("w0", &w0, "w0/i"); h->Branch("w1", &w1, "w1/i"); h->Branch("w2", &w2, "w2/i"); h->Branch("w3", &w3, "w3/i");
  h->Branch("marker", &marker, "marker/i"); h->Branch("event_size", &evsize, "event_size/i"); h->Branch("board_id", &board_id, "board_id/i"); h->Branch("lvds_pattern", &lvds, "lvds_pattern/i"); h->Branch("group_mask", &gmask, "group_mask/i"); h->Branch("w1_reserved", &w1_res, "w1_reserved/i");
  h->Branch("event_counter", &evnr, "event_counter/i"); h->Branch("w2_reserved", &w2_res, "w2_reserved/i"); h->Branch("event_time_tag", &evttime, "event_time_tag/i");
  h->Branch("g_hdr", g_hdr, "g_hdr[2]/i"); h->Branch("g_index_cell", g_idxcell, "g_index_cell[2]/i"); h->Branch("g_freq", g_freq, "g_freq[2]/i"); h->Branch("g_tr_present", g_tr, "g_tr_present[2]/i"); h->Branch("g_sample_words", g_nwords, "g_sample_words[2]/i"); h->Branch("g_hdr_reserved", g_res, "g_hdr_reserved[2]/i");
  h->Branch("g_tt_word", g_ttword, "g_tt_word[2]/i"); h->Branch("g_trigger_time", g_tt, "g_trigger_time[2]/i"); h->Branch("g_tt_reserved", g_ttres, "g_tt_reserved[2]/i"); h->Branch("gtt", &gtt, "gtt/L");
  h->Branch("ok_size", &ok_size, "ok_size/O"); h->Branch("ok_marker", &ok_marker, "ok_marker/O");
  // channel trees
  Short_t adc[16][256]; Int_t bb_ch[12], pad_ch[4];
  TTree *b = new TTree("bbcal", "Baby BCAL South, CAEN ch 0-11 (S A0, A1, A3, A4, A5, 6, 7, 8, 9, A14, 10, 11): adc[i][256], 12-bit DRS4 samples at 1 GS/s, negative pulses, pedestal ~3450-3530");
  for (int i = 0; i < 12; i++) bb_ch[i] = i; b->Branch("ch", bb_ch, "ch[12]/I"); b->Branch("adc", adc[0], "adc[12][256]/S");
  TTree *p = new TTree("paddles", "trigger paddles P3-P6 = CAEN ch 12-15 (negative pulses, pedestal ~3720-3770, rail at 0 on large pulses)");
  const char *PN[4] = {"P3", "P4", "P5", "P6"};
  for (int k = 0; k < 4; k++) { pad_ch[k] = 12 + k; p->Branch(Form("%s_ch", PN[k]), &pad_ch[k], Form("%s_ch/I", PN[k])); p->Branch(Form("%s_adc", PN[k]), adc[12 + k], Form("%s_adc[256]/S", PN[k])); }
  // packets
  Int_t p_seg, p_evt, p_type, p_run, p_evtlen, p_pktlen, p_datalen; Long64_t p_time, p_rt; Bool_t p_has;
  TTree *tp = new TTree("packets", "one entry per RCDAQ event of this segment, in file order");
  tp->Branch("seg", &p_seg, "seg/I"); tp->Branch("evt", &p_evt, "evt/I"); tp->Branch("evt_type", &p_type, "evt_type/I"); tp->Branch("run", &p_run, "run/I"); tp->Branch("evt_time", &p_time, "evt_time/L"); tp->Branch("evt_len", &p_evtlen, "evt_len/I"); tp->Branch("has_caen", &p_has, "has_caen/O"); tp->Branch("pkt_len", &p_pktlen, "pkt_len/I"); tp->Branch("pkt_datalen", &p_datalen, "pkt_datalen/I"); tp->Branch("rt", &p_rt, "rt/L");

  int status = 0; fileEventiterator it(argv[2], status); if (status) { fprintf(stderr, "cannot open %s\n", argv[2]); return 1; }
  long nev = 0, ntrig = 0, badsize = 0, badmarker = 0, mismatch = 0, nsamp_cmp = 0;
  while (Event *e = it.getNextEvent())
    {
      nev++;
      p_seg = seg; p_evt = e->getEvtSequence(); p_type = e->getEvtType(); p_run = e->getRunNumber(); p_time = (Long64_t) e->getTime(); p_evtlen = e->getEvtLength(); p_has = false; p_pktlen = p_datalen = 0; p_rt = 0;
      Long64_t mono = 0, diff = 0;
      for (int id = 1010; id <= 1013; id++) if (Packet *q = e->getPacket(id)) { p_rt = q->lValue(0, "TOTAL"); mono = q->lValue(2, "TOTAL"); diff = q->lValue(0, "TIMEDIFF"); delete q; break; }
      Packet *pk = (p_type == DATAEVENT) ? e->getPacket(pid) : nullptr;
      if (pk)
        {
          p_has = true; p_pktlen = pk->getLength(); p_datalen = pk->getDataLength();
          std::vector<int> raw(pk->getLength()); pk->copyMe(raw.data(), raw.size());
          const uint32_t *w = reinterpret_cast<const uint32_t *>(raw.data()) + 4; const int nw = pk->getDataLength();
          seg_b = seg; run = p_run; idx = ntrig; evt = p_evt; evt_type = p_type; evt_time = p_time; rt = p_rt; rt_mono = mono; rt_diff = diff;
          w0 = nw > 0 ? w[0] : 0; w1 = nw > 1 ? w[1] : 0; w2 = nw > 2 ? w[2] : 0; w3 = nw > 3 ? w[3] : 0;
          marker = w0 >> 28; evsize = w0 & 0x0fffffff; board_id = w1 >> 27; w1_res = ((w1 >> 24) & 7) | (((w1 >> 4) & 0xf) << 3); lvds = (w1 >> 8) & 0xffff; gmask = w1 & 0xf; evnr = w2 & 0x3fffff; w2_res = w2 >> 22; evttime = w3;
          ok_marker = marker == 0xA; ok_size = (int) evsize == nw && nw == 1544; if (!ok_marker) badmarker++; if (!ok_size) badsize++;
          memset(adc, 0, sizeof(adc)); for (int g = 0; g < 2; g++) { g_hdr[g] = g_idxcell[g] = g_freq[g] = g_tr[g] = g_nwords[g] = g_res[g] = g_ttword[g] = g_tt[g] = g_ttres[g] = 0; }
          int pos = 4;
          for (int g = 0; g < 2 && pos + 770 <= nw; g++)
            {
              const uint32_t gh = w[pos]; g_hdr[g] = gh; g_idxcell[g] = (gh >> 20) & 0x3ff; g_freq[g] = (gh >> 16) & 3; g_tr[g] = (gh >> 12) & 1; g_nwords[g] = gh & 0xfff; g_res[g] = (gh >> 30) | (((gh >> 18) & 3) << 2) | (((gh >> 13) & 7) << 4);
              const int nsw = g_nwords[g] == 768 ? 768 : std::min<int>(g_nwords[g], nw - pos - 2);
              for (int s = 0; s < nsw / 3; s++)
                {
                  const uint32_t W0 = w[pos + 1 + 3 * s], W1 = w[pos + 2 + 3 * s], W2 = w[pos + 3 + 3 * s];
                  Short_t v[8] = {Short_t(W0 & 0xfff), Short_t((W0 >> 12) & 0xfff), Short_t((W0 >> 24) | ((W1 & 0xf) << 8)), Short_t((W1 >> 4) & 0xfff), Short_t((W1 >> 16) & 0xfff), Short_t((W1 >> 28) | ((W2 & 0xff) << 4)), Short_t((W2 >> 8) & 0xfff), Short_t(W2 >> 20)};
                  if (s < 256) for (int k = 0; k < 8; k++) adc[8 * g + k][s] = v[k];
                }
              const uint32_t tw = w[pos + 1 + nsw]; g_ttword[g] = tw; g_tt[g] = tw & 0x3fffffff; g_ttres[g] = tw >> 30;
              pos += 2 + nsw;
            }
          gtt = (Long64_t) g_tt[0] | ((Long64_t) g_tt[1] << 30);
          // cross-check with the stock decoder
          if (ntrig % 50 == 0) { for (int ch = 0; ch < 16; ch += 5) for (int s = 0; s < 256; s += 7) { nsamp_cmp++; if (pk->iValue(s, ch) != adc[ch][s]) mismatch++; } if (pk->lValue(0, "GROUPTRIGGERTIME") != gtt) mismatch++; }
          delete pk;
          h->Fill(); b->Fill(); p->Fill(); ntrig++;
        }
      tp->Fill(); delete e;
      if (maxev && nev >= maxev) break;
    }
  TString info = Form("caen_to_root: %s (segment %d): %ld RCDAQ events, %ld triggers; bad size %ld, bad marker %ld; stock-decoder cross-check: %ld sample comparisons, %ld mismatches (incl. GROUPTRIGGERTIME). gtt = g_trigger_time[0] | g_trigger_time[1] << 30 (8.533 ns ticks); event_time_tag is the broken 17.07 ns tag.", argv[2], seg, nev, ntrig, badsize, badmarker, nsamp_cmp, mismatch);
  f->cd(); TNamed(Form("info_seg%d", seg), info.Data()).Write(); h->Write(); b->Write(); p->Write(); tp->Write(); f->Close();
  printf("%s\n", info.Data());
  return 0;
}
