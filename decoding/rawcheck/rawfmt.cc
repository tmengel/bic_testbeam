// rawfmt: inspect raw words of the BIC streams to document the format.
//   rawfmt kcu  file.evt nevt
//   rawfmt caen file.evt nevt
//   rawfmt astro file.evt nevt
#include <Event/Event.h>
#include <Event/Eventiterator.h>
#include <Event/fileEventiterator.h>
#include <Event/packet.h>
#include <Event/EventTypes.h>

#include <cstdio>
#include <cstdint>
#include <map>
#include <string>
#include <vector>

static uint32_t le32(const uint8_t *v) { return (uint32_t(v[3]) << 24) | (uint32_t(v[2]) << 16) | (uint32_t(v[1]) << 8) | v[0]; }
static uint32_t be32(const uint8_t *v) { return (uint32_t(v[0]) << 24) | (uint32_t(v[1]) << 16) | (uint32_t(v[2]) << 8) | v[3]; }
static uint64_t be64(const uint8_t *v) { uint64_t r = 0; for (int i = 0; i < 8; i++) r = (r << 8) | v[i]; return r; }

static void hexdump(const uint8_t *c, int n)
{
  for (int i = 0; i < n; i++) { printf("%02x", c[i]); if (i % 4 == 3) printf(" "); if (i % 32 == 31) printf("\n"); }
  printf("\n");
}

template <class M> static void show(const char *name, const M &m, int maxn = 40)
{
  printf("%s (%zu distinct):", name, m.size());
  int n = 0;
  for (auto &kv : m) { if (n++ >= maxn) { printf(" ..."); break; } printf(" %llx:%ld", (unsigned long long) kv.first, kv.second); }
  printf("\n");
}

int main(int argc, char **argv)
{
  std::string mode = argv[1];
  int status = 0;
  fileEventiterator it(argv[2], status);
  long nmax = atol(argv[3]);
  long nev = 0;

  // kcu
  std::map<uint64_t, long> h_d2, h_d3, h_udpmark, h_udphdr8_21, h_rochdr, h_cmn_hi, h_calib_hi, h_spare, h_tctp, h_totfield_hi, h_lines_per_udp, h_udplen;
  long nwords = 0, n_toa = 0, n_tot = 0, n_tot_hi = 0, n_adc_any = 0, n_bits1916 = 0;
  std::vector<long> toa_hist(1024, 0), tot10_hist(1024, 0);
  long n_both = 0;
  // caen
  std::map<uint64_t, long> h_w0top, h_mask, h_board, h_freq, h_size, h_trpres, h_gdiff;
  // astro
  std::map<uint64_t, long> h_len, h_layer, h_id, h_payload, h_b3bit6, h_iscol;
  std::vector<long> atot(4096, 0), ats(256, 0);
  int printed = 0;

  while (Event *e = it.getNextEvent())
  {
    if (e->getEvtType() != DATAEVENT) { delete e; continue; }
    if (nev++ >= nmax) { delete e; break; }
    for (Packet *p : e->getPacketVector())
    {
      const int hf = p->getHitFormat();
      std::vector<int> raw(p->getLength());
      p->copyMe(raw.data(), raw.size());
      const uint8_t *c = reinterpret_cast<const uint8_t *>(raw.data()) + 16;
      const long nbytes = 4L * p->getDataLength();
      if (mode == "kcu" && hf == 302)
      {
        if (printed < 1) { printf("packet id %d len %d words, hdr: ", p->getIdentifier(), p->getLength()); hexdump(reinterpret_cast<const uint8_t *>(raw.data()), 16); printf("first 22+192*2 payload bytes:\n"); hexdump(c, 22 + 384); printed++; }
        long pos = 0;
        while (pos + 22 <= nbytes)
        {
          const long plen = long(le32(c + pos + 4)) + 8;
          h_udpmark[le32(c + pos)]++;
          if (le32(c + pos) != 0xaabbccdd || plen < 22 || pos + plen > nbytes) break;
          h_udplen[plen]++;
          h_lines_per_udp[(plen - 22) / 192 * 1000 + (plen - 22) % 192]++;
          for (long l = pos + 22; l + 192 <= pos + plen; l += 192)
          {
            const uint8_t *d = c + l;
            if (d[0] != 0xaa || d[1] != 0x5a) continue;
            h_d2[d[2]]++; h_d3[d[3]]++;
            h_spare[be64(d + 24)]++;
            h_rochdr[be32(d + 32) & 0xf0000fff]++;
            h_cmn_hi[be32(d + 36) >> 20 << 20 ? (be32(d + 36) & 0xf0000000u) | 0x1 : 0]++;
            h_calib_hi[be32(d + 112) & 0xf0000000u]++;
            for (int k = 0; k < 36; k++)
            {
              const uint32_t w = be32(d + (k < 18 ? 40 + 4 * k : 116 + 4 * (k - 18)));
              nwords++;
              const int tctp = w >> 30, adc = (w >> 20) & 0x3ff, tot10 = (w >> 10) & 0x3ff, toa = w & 0x3ff;
              h_tctp[tctp]++;
              if (adc) n_adc_any++;
              if (toa) { n_toa++; toa_hist[toa]++; }
              if (tot10) { n_tot++; tot10_hist[tot10]++; }
              if (tot10 >> 6) n_bits1916++;
              if (toa && tot10) n_both++;
            }
          }
          pos += plen;
        }
      }
      if (mode == "caen" && hf != 9)
      {
        const uint32_t *w = reinterpret_cast<const uint32_t *>(c);
        if (printed < 2) { printf("packet id %d hitformat %d len %d\n", p->getIdentifier(), hf, p->getLength()); for (int i = 0; i < 6; i++) printf(" w%d=%08x", i, w[i]); printf("\n"); printed++; }
        h_w0top[w[0] >> 28]++; h_mask[w[1] & 0xf]++; h_board[w[1] >> 27]++;
        // group 0 header at w[4]
        h_freq[(w[4] >> 16) & 3]++; h_size[w[4] & 0xfff]++; h_trpres[(w[4] >> 12) & 1]++;
        long g0 = (long) (p->lValue(0, "GROUPTRIGGERTIME") & 0x3fffffff);
        // raw group1 tag: find via size
        int sz = w[4] & 0xfff;
        uint32_t t0 = w[4 + 1 + sz], t1 = w[4 + 1 + sz + 1 + sz + 1];
        h_gdiff[(long) (t1 & 0x3fffffff) - (long) (t0 & 0x3fffffff)]++;
        if (printed < 6) { printf(" evtTT=%u g0tag=%08x g1hdr=%08x g1tag=%08x  GTT=%lld  stitched_g0=%ld\n", w[3], t0, w[4 + 1 + sz + 1], t1, p->lValue(0, "GROUPTRIGGERTIME"), g0); printed++; }
      }
      if (mode == "astro" && hf != 9)
      {
        const uint8_t *b = c;
        const long n = nbytes - p->getPadding();
        if (printed < 1) { printf("packet id %d hitformat %d datalen %d pad %d\n", p->getIdentifier(), hf, p->getDataLength(), p->getPadding()); hexdump(b, 66); printed++; }
        long i = 0;
        while (i < n)
        {
          const int l = b[i] + 1;
          h_len[l]++;
          if (l == 11 && i + l <= n)
          {
            h_layer[b[i + 1]]++;
            if (b[i + 1] == 1) {
              h_id[b[i + 2] >> 3]++; h_payload[b[i + 2] & 7]++; h_b3bit6[(b[i + 3] >> 6) & 1]++; h_iscol[b[i + 3] >> 7]++;
              ats[b[i + 4]]++; atot[((b[i + 5] & 0xf) << 8) | b[i + 6]]++;
              if (b[i + 5] >> 4) h_b3bit6[0x100 | (b[i + 5] >> 4)]++;
            }
          }
          i += l;
        }
      }
      delete p;
    }
    delete e;
  }
  printf("events %ld\n", nev);
  if (mode == "kcu")
  {
    show("udp marker", h_udpmark); show("udp len", h_udplen); show("lines*1000+rem", h_lines_per_udp);
    show("d[2]", h_d2); show("d[3]", h_d3); show("spare d[24..31]", h_spare, 10); show("rochdr&0xf0000fff", h_rochdr, 20);
    show("cmn top", h_cmn_hi); show("calib top", h_calib_hi);
    show("Tc/Tp", h_tctp);
    printf("channel words %ld, adc!=0 %ld, toa!=0 %ld (%.4f%%), tot10!=0 %ld (%.4f%%), tot bits19:16!=0 %ld, toa&tot %ld\n", nwords, n_adc_any, n_toa, 100.0 * n_toa / nwords, n_tot, 100.0 * n_tot / nwords, n_bits1916, n_both);
    // coarse histos
    printf("TOA hist (64 bins of 16):"); for (int b = 0; b < 64; b++) { long s = 0; for (int k = 0; k < 16; k++) s += toa_hist[b * 16 + k]; printf(" %ld", s); } printf("\n");
    printf("TOT10 hist (64 bins of 16):"); for (int b = 0; b < 64; b++) { long s = 0; for (int k = 0; k < 16; k++) s += tot10_hist[b * 16 + k]; printf(" %ld", s); } printf("\n");
    printf("TOT10 low values:"); for (int k = 0; k < 80; k++) printf(" %d:%ld", k, tot10_hist[k]); printf("\n");
  }
  if (mode == "caen") { show("w0>>28", h_w0top); show("groupmask", h_mask); show("w1>>27 (board)", h_board); show("freq", h_freq); show("g0 size", h_size); show("trpresent", h_trpres); show("g1tag-g0tag", h_gdiff); }
  if (mode == "astro")
  {
    show("rec len", h_len); show("layer", h_layer); show("id(layer1)", h_id); show("payload", h_payload); show("b3bit6 / b5hi", h_b3bit6); show("iscol", h_iscol);
    printf("ts hist (32 bins of 8):"); for (int b = 0; b < 32; b++) { long s = 0; for (int k = 0; k < 8; k++) s += ats[b * 8 + k]; printf(" %ld", s); } printf("\n");
    long tot = 0, nz = 0, mx = 0; for (int k = 0; k < 4096; k++) { tot += atot[k]; if (atot[k]) mx = k; }
    printf("ToT: max %ld; hist (64 bins of 64):", mx); for (int b = 0; b < 64; b++) { long s = 0; for (int k = 0; k < 64; k++) s += atot[b * 64 + k]; printf(" %ld", s); } printf("\n");
    (void) nz;
  }
  return 0;
}
