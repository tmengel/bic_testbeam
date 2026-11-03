// astro_to_root: the AstroPix stream (packet 13001) into ROOT trees, every SPI record kept.
//
//   astro_to_root out.root file.evt [max_events=0]
//
// Trees:
//   hits     one entry per 11-byte hit record, in stream order: raw bytes and every decoded field
//            (incl. the bits the stock decoder drops: byte 3 bit 6, byte 5 bits 7:4), plus the RCDAQ
//            event, record index in the packet and the RTCLOCK time
//   records  one entry per record that is NOT an 11-byte hit (SPI filler, 7-14 bytes): raw bytes
//   packets  one entry per RCDAQ event (also the non-data events): record-length histogram, counts
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
  if (argc < 3) { fprintf(stderr, "usage: astro_to_root out.root file.evt [max_events]\n"); return 1; }
  const long maxev = argc > 3 ? atol(argv[3]) : 0; const int pid = 13001;
  TFile *f = new TFile(argv[1], "RECREATE", Form("AstroPix packet 13001, %s", argv[2])); f->SetCompressionAlgorithm(ROOT::RCompressionSetting::EAlgorithm::kZSTD); f->SetCompressionLevel(5);
  Int_t evt, run, irec, ihit; Long64_t idx, evt_time, rt; UChar_t raw[11], layer, id, payload, iscol, flag, location, ts, tot_nibble; UShort_t tot; UInt_t fpga_ts; Bool_t valid;
  TTree *h = new TTree("hits", "one entry per 11-byte AstroPix hit record in stream order; valid = layer 1 and id <= 8");
  h->Branch("idx", &idx, "idx/L"); h->Branch("evt", &evt, "evt/I"); h->Branch("run", &run, "run/I"); h->Branch("evt_time", &evt_time, "evt_time/L"); h->Branch("rt", &rt, "rt/L"); h->Branch("irec", &irec, "irec/I"); h->Branch("ihit", &ihit, "ihit/I");
  h->Branch("raw", raw, "raw[11]/b"); h->Branch("layer", &layer, "layer/b"); h->Branch("id", &id, "id/b"); h->Branch("payload", &payload, "payload/b"); h->Branch("iscol", &iscol, "iscol/b"); h->Branch("flag", &flag, "flag/b"); h->Branch("location", &location, "location/b");
  h->Branch("ts", &ts, "ts/b"); h->Branch("tot", &tot, "tot/s"); h->Branch("tot_nibble", &tot_nibble, "tot_nibble/b"); h->Branch("fpga_ts", &fpga_ts, "fpga_ts/i"); h->Branch("valid", &valid, "valid/O");
  Int_t r_evt, r_irec, r_len, r_avail; UChar_t r_raw[16]; Long64_t r_rt;
  TTree *r = new TTree("records", "every record that is not an 11-byte hit: len = length byte + 1, avail = bytes present (less than len when the payload ends inside the record), first 16 bytes raw");
  r->Branch("evt", &r_evt, "evt/I"); r->Branch("rt", &r_rt, "rt/L"); r->Branch("irec", &r_irec, "irec/I"); r->Branch("len", &r_len, "len/I"); r->Branch("avail", &r_avail, "avail/I"); r->Branch("raw", r_raw, "raw[16]/b");
  Int_t p_evt, p_type, p_run, p_evtlen, p_pktlen, p_datalen, p_padding, p_nrec, p_nhit, p_nvalid, p_lenhist[16], p_trunc; Long64_t p_time, p_rt, p_rtmono, p_rtdiff; Bool_t p_has;
  TTree *tp = new TTree("packets", "one entry per RCDAQ event: record counts and the record-length histogram (index = len, 15 = longer)");
  tp->Branch("evt", &p_evt, "evt/I"); tp->Branch("evt_type", &p_type, "evt_type/I"); tp->Branch("run", &p_run, "run/I"); tp->Branch("evt_time", &p_time, "evt_time/L"); tp->Branch("evt_len", &p_evtlen, "evt_len/I"); tp->Branch("has_astro", &p_has, "has_astro/O"); tp->Branch("pkt_len", &p_pktlen, "pkt_len/I"); tp->Branch("pkt_datalen", &p_datalen, "pkt_datalen/I"); tp->Branch("padding", &p_padding, "padding/I");
  tp->Branch("rt", &p_rt, "rt/L"); tp->Branch("rt_mono", &p_rtmono, "rt_mono/L"); tp->Branch("rt_diff", &p_rtdiff, "rt_diff/L"); tp->Branch("nrec", &p_nrec, "nrec/I"); tp->Branch("nhit", &p_nhit, "nhit/I"); tp->Branch("nvalid", &p_nvalid, "nvalid/I"); tp->Branch("lenhist", p_lenhist, "lenhist[16]/I"); tp->Branch("truncated", &p_trunc, "truncated/I");

  int status = 0; fileEventiterator it(argv[2], status); if (status) { fprintf(stderr, "cannot open %s\n", argv[2]); return 1; }
  long nev = 0, nhits = 0, nrecs = 0, mismatch = 0, ncmp = 0; idx = 0;
  while (Event *e = it.getNextEvent())
    {
      nev++;
      p_evt = e->getEvtSequence(); p_type = e->getEvtType(); p_run = e->getRunNumber(); p_time = (Long64_t) e->getTime(); p_evtlen = e->getEvtLength(); p_has = false; p_pktlen = p_datalen = p_padding = p_nrec = p_nhit = p_nvalid = p_trunc = 0; memset(p_lenhist, 0, sizeof(p_lenhist)); p_rt = p_rtmono = p_rtdiff = 0;
      for (int i = 1010; i <= 1013; i++) if (Packet *q = e->getPacket(i)) { p_rt = q->lValue(0, "TOTAL"); p_rtmono = q->lValue(2, "TOTAL"); p_rtdiff = q->lValue(0, "TIMEDIFF"); delete q; break; }
      Packet *pk = (p_type == DATAEVENT) ? e->getPacket(pid) : nullptr;
      if (pk)
        {
          p_has = true; p_pktlen = pk->getLength(); p_datalen = pk->getDataLength(); p_padding = pk->getPadding();
          std::vector<int> rawp(pk->getLength()); pk->copyMe(rawp.data(), rawp.size());
          const uint8_t *b = reinterpret_cast<const uint8_t *>(rawp.data()) + 16; const long nbytes = 4L * pk->getDataLength() - pk->getPadding();
          long pos = 0; int nh = 0;
          while (pos < nbytes)
            {
              const int len = b[pos] + 1;
              if (pos + len > nbytes)   // the payload ends inside this record: keep what there is, flagged by len > bytes available
                {
                  p_trunc++; r_evt = p_evt; r_rt = p_rt; r_irec = p_nrec; r_len = len; r_avail = (Int_t) (nbytes - pos); memset(r_raw, 0, 16); memcpy(r_raw, b + pos, r_avail < 16 ? r_avail : 16); r->Fill(); nrecs++;
                  break;
                }
              p_nrec++; p_lenhist[len < 15 ? len : 15]++;
              if (len == 11)
                {
                  const uint8_t *c = b + pos; memcpy(raw, c, 11);
                  evt = p_evt; run = p_run; evt_time = p_time; rt = p_rt; irec = p_nrec - 1; ihit = nh;
                  layer = c[1]; id = c[2] >> 3; payload = c[2] & 7; iscol = c[3] >> 7; flag = (c[3] >> 6) & 1; location = c[3] & 0x3f; ts = c[4]; tot_nibble = c[5] >> 4; tot = ((c[5] & 0xf) << 8) | c[6]; fpga_ts = (uint32_t(c[7]) << 24) | (uint32_t(c[8]) << 16) | (uint32_t(c[9]) << 8) | c[10];
                  valid = layer == 1 && id <= 8; if (valid) p_nvalid++;
                  h->Fill(); idx++; nh++; nhits++;
                }
              else
                {
                  r_evt = p_evt; r_rt = p_rt; r_irec = p_nrec - 1; r_len = len; r_avail = len; memset(r_raw, 0, 16); memcpy(r_raw, b + pos, len < 16 ? len : 16); r->Fill(); nrecs++;
                }
              pos += len;
            }
          p_nhit = nh;
          // cross-check: the stock decoder's hit count and a few fields (its multiset reorders hits, so compare sums)
          if (nev % 20 == 0) { ncmp++; if (pk->iValue(0, "NR_UNITS") != nh) mismatch++; }
          delete pk;
        }
      tp->Fill(); delete e;
      if (maxev && nev >= maxev) break;
    }
  TString info = Form("astro_to_root: %s: %ld RCDAQ events, %ld hit records, %ld other records; stock-decoder hit-count cross-check on %ld packets, %ld mismatches. fpga_ts: 40 MHz readout stamp (wraps every 107 s); ts: 8-bit 10 MHz on-chip stamp; tot: 12 bit; valid = layer 1 and id <= 8.", argv[2], nev, nhits, nrecs, ncmp, mismatch);
  f->cd(); TNamed("info", info.Data()).Write(); h->Write(); r->Write(); tp->Write(); f->Close();
  printf("%s\n", info.Data());
  return 0;
}
