// kcu2: deeper look at H2GCROC3 line words, and TOA/TOT inside stock-decoded waveforms.
//   kcu2 file.evt nevt
#include <Event/Event.h>
#include <Event/Eventiterator.h>
#include <Event/fileEventiterator.h>
#include <Event/packet.h>
#include <Event/EventTypes.h>

#include <cstdio>
#include <cstdint>
#include <map>
#include <vector>

static uint32_t le32(const uint8_t *v) { return (uint32_t(v[3]) << 24) | (uint32_t(v[2]) << 16) | (uint32_t(v[1]) << 8) | v[0]; }
static uint32_t be32(const uint8_t *v) { return (uint32_t(v[0]) << 24) | (uint32_t(v[1]) << 16) | (uint32_t(v[2]) << 8) | v[3]; }
static uint64_t be64(const uint8_t *v) { uint64_t r = 0; for (int i = 0; i < 8; i++) r = (r << 8) | v[i]; return r; }

template <class M> static void show(const char *name, const M &m, int maxn = 40)
{
  printf("%s (%zu distinct):", name, m.size());
  int n = 0;
  for (auto &kv : m) { if (n++ >= maxn) { printf(" ..."); break; } printf(" %lld:%ld", (long long) kv.first, kv.second); }
  printf("\n");
}

int main(int argc, char **argv)
{
  int status = 0;
  fileEventiterator it(argv[1], status);
  long nmax = atol(argv[2]);
  long nev = 0;
  std::map<long long, long> h_hdr_top, h_hdr_low, h_hdr_err, h_hdr_ev6, h_bx_minus, h_cmn_a, h_cmn_b, h_cmn_top, h_cal_adc, h_cal_tctp,
      h_tctp_vs_tot, h_udp8, h_udp12, h_udp16, h_udp20, h_evctr_minus_trigout, h_trigin, h_ev6_minus_evctr;
  // waveform level
  std::map<long long, long> h_toa_sample, h_tot_sample, h_peak_minus_toa, h_ntoa_per_ch, h_tot6_vals, h_adc_at_tot, h_maxadc_tot, h_toa_per_chip;
  long n_hitch = 0, n_hitch_toa = 0, n_hitch_tot = 0, n_ped_toa = 0, n_ped_samples = 0;
  std::map<int, long> ch_toa, ch_tot;
  uint32_t prev_udp8 = 0;
  bool have_prev = false;

  while (Event *e = it.getNextEvent())
  {
    if (e->getEvtType() != DATAEVENT) { delete e; continue; }
    if (nev++ >= nmax) { delete e; break; }
    Packet *p = e->getPacket(12001);
    if (!p) p = e->getPacket(12002);
    if (!p) { delete e; continue; }
    std::vector<int> raw(p->getLength());
    p->copyMe(raw.data(), raw.size());
    const uint8_t *c = reinterpret_cast<const uint8_t *>(raw.data()) + 16;
    const long nbytes = 4L * p->getDataLength();
    long pos = 0;
    while (pos + 22 <= nbytes)
    {
      const long plen = long(le32(c + pos + 4)) + 8;
      if (le32(c + pos) != 0xaabbccdd || pos + plen > nbytes) break;
      const uint32_t u8 = le32(c + pos + 8);
      if (have_prev) h_udp8[(long long) u8 - prev_udp8]++;
      prev_udp8 = u8; have_prev = true;
      h_udp12[le32(c + pos + 12) >> 24]++;
      h_udp16[le32(c + pos + 16)]++;
      h_udp20[c[pos + 20] * 256 + c[pos + 21]]++;
      for (long l = pos + 22; l + 192 <= pos + plen; l += 192)
      {
        const uint8_t *d = c + l;
        if (d[0] != 0xaa || d[1] != 0x5a) continue;
        const uint32_t trin = be32(d + 4), trout = be32(d + 8), evctr = be32(d + 12);
        const uint64_t ts = be64(d + 16);
        if (ts == 0) continue;
        h_trigin[trin]++;
        h_evctr_minus_trigout[(long long) evctr - trout]++;
        const uint32_t hdr = be32(d + 32), cmn = be32(d + 36), cal = be32(d + 112);
        h_hdr_top[hdr >> 28]++; h_hdr_low[hdr & 0xf]++; h_hdr_err[(hdr >> 4) & 0x7]++;
        const int bx = (hdr >> 16) & 0xfff, ev6 = (hdr >> 10) & 0x3f;
        h_bx_minus[((long long) (ts / 4) - bx) % 3564]++;
        h_ev6_minus_evctr[(ev6 - (long long) evctr) & 0x3f]++;
        h_hdr_ev6[ev6 >= 0 ? 0 : 1]++;
        h_cmn_top[cmn >> 20]++; h_cmn_a[(cmn >> 10) & 0x3ff]++; h_cmn_b[cmn & 0x3ff]++;
        h_cal_adc[(cal >> 20) & 0x3ff]++; h_cal_tctp[cal >> 30]++;
        for (int k = 0; k < 36; k++)
        {
          const uint32_t w = be32(d + (k < 18 ? 40 + 4 * k : 116 + 4 * (k - 18)));
          const int tot = (w >> 10) & 0x3ff;
          h_tctp_vs_tot[(w >> 30) * 10 + (tot ? 1 : 0)]++;
        }
      }
      pos += plen;
    }

    // waveform level, stock decoder, complete 20-sample waveforms
    const int nwf = p->iValue(0, "NR_WF");
    const int nch = p->iValue(0, "CHANNELS");
    for (int w = 0; w < nwf; w++)
    {
      if (p->iValue(w, "SAMPLESIZE") != 20) continue;
      for (int ch = 0; ch < nch; ch++)
      {
        int a[20], toa[20], tot[20];
        int mx = -1, smx = 0;
        for (int s = 0; s < 20; s++) { a[s] = p->iValue(w, s, ch); toa[s] = p->iValue(w, s, ch, "TOA"); tot[s] = p->iValue(w, s, ch, "TOT"); if (a[s] > mx) { mx = a[s]; smx = s; } }
        const float ped = (a[1] + a[2] + a[3]) / 3.f;
        if (ped > 1000) continue;
        int ntoa = 0, stoa = -1, ntot = 0, stot = -1;
        for (int s = 0; s < 20; s++)
        {
          if (toa[s]) { ntoa++; if (stoa < 0) stoa = s; }
          if (tot[s]) { ntot++; if (stot < 0) stot = s; h_tot6_vals[tot[s]]++; }
        }
        if (ntoa) ch_toa[ch]++;
        if (ntot) ch_tot[ch]++;
        const bool hit = mx - ped > 50;
        if (!hit) { for (int s = 0; s < 4; s++) { n_ped_samples++; if (toa[s]) n_ped_toa++; } continue; }
        n_hitch++;
        if (ntoa) { n_hitch_toa++; h_toa_sample[stoa]++; h_peak_minus_toa[smx - stoa]++; }
        h_ntoa_per_ch[ntoa]++;
        if (ntot) { n_hitch_tot++; h_tot_sample[stot]++; h_maxadc_tot[mx / 100 * 100]++; }
      }
    }
    delete p;
    delete e;
  }
  printf("events %ld\n", nev);
  show("UDP bytes 8-11 delta", h_udp8, 10); show("UDP byte 15", h_udp12, 10); show("UDP 16-19", h_udp16, 10); show("UDP 20-21", h_udp20, 10);
  show("trig_in", h_trigin, 5); show("evctr - trig_out", h_evctr_minus_trigout, 10);
  show("ROC hdr [31:28]", h_hdr_top); show("ROC hdr [3:0]", h_hdr_low); show("ROC hdr [6:4]", h_hdr_err);
  show("(ts/4 - BX) mod 3564", h_bx_minus, 10);
  show("(ev6 - evctr) mod 64", h_ev6_minus_evctr, 10);
  show("CMN [31:20]", h_cmn_top, 10); show("CMN [19:10]", h_cmn_a, 60); show("CMN [9:0]", h_cmn_b, 10);
  show("CALIB ADC", h_cal_adc, 30); show("CALIB TcTp", h_cal_tctp);
  show("TcTp*10 + (tot!=0)", h_tctp_vs_tot);
  printf("hit channels (max-ped>50): %ld, with TOA %ld (%.1f%%), with TOT %ld (%.2f%%)\n", n_hitch, n_hitch_toa, 100.0 * n_hitch_toa / n_hitch, n_hitch_tot, 100.0 * n_hitch_tot / n_hitch);
  printf("pedestal samples 0-3 of non-hit channels: %ld, TOA!=0 %ld (%.3f%%)\n", n_ped_samples, n_ped_toa, 100.0 * n_ped_toa / n_ped_samples);
  show("first TOA sample (hit ch)", h_toa_sample); show("peak sample - first TOA sample", h_peak_minus_toa); show("# samples with TOA (hit ch)", h_ntoa_per_ch);
  show("first TOT sample", h_tot_sample); show("max ADC (100s) when TOT", h_maxadc_tot); show("stock TOT (6-bit) values", h_tot6_vals, 64);
  printf("channels with TOA in any wf (top): "); for (auto &kv : ch_toa) if (kv.second > 1000) printf(" %d:%ld", kv.first, kv.second); printf("\n");
  printf("channels with TOT: "); for (auto &kv : ch_tot) printf(" %d:%ld", kv.first, kv.second); printf("\n");
  return 0;
}
