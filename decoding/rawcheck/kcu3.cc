// kcu3: full 10-bit TOT from raw lines vs the stock 6-bit TOT, against the number of saturated ADC samples;
// also dumps (ts, BX, ev6, evctr, trig_out) of raw lines for the header study.
//   kcu3 file.evt nevt tot.txt hdr.txt
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

int main(int argc, char **argv)
{
  int status = 0;
  fileEventiterator it(argv[1], status);
  long nmax = atol(argv[2]);
  FILE *ft = fopen(argv[3], "w"), *fh = fopen(argv[4], "w");
  long nev = 0, nhdr = 0;
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
    // (ts, chip*2+group) -> 36 full words
    std::multimap<uint64_t, std::pair<int, std::vector<uint32_t>>> lines;
    long pos = 0;
    while (pos + 22 <= nbytes)
    {
      const long plen = long(le32(c + pos + 4)) + 8;
      if (le32(c + pos) != 0xaabbccdd || pos + plen > nbytes) break;
      for (long l = pos + 22; l + 192 <= pos + plen; l += 192)
      {
        const uint8_t *d = c + l;
        if (d[0] != 0xaa || d[1] != 0x5a) continue;
        const uint64_t ts = be64(d + 16);
        if (ts == 0) continue;
        const int link = (d[2] & 0xf) * 2 + (d[3] - 36);
        const uint32_t hdr = be32(d + 32);
        if (nhdr < 20000) { fprintf(fh, "%llu %d %u %u %u %u\n", (unsigned long long) ts, link, (hdr >> 16) & 0xfff, (hdr >> 10) & 0x3f, be32(d + 12), be32(d + 8)); nhdr++; }
        std::vector<uint32_t> w(36);
        bool anytot = false;
        for (int k = 0; k < 36; k++) { w[k] = be32(d + (k < 18 ? 40 + 4 * k : 116 + 4 * (k - 18))); if ((w[k] >> 10) & 0x3ff) anytot = true; }
        if (anytot) lines.emplace(ts, std::make_pair(link, w));
      }
      pos += plen;
    }
    const int nwf = p->iValue(0, "NR_WF");
    for (int wf = 0; wf < nwf; wf++)
    {
      if (p->iValue(wf, "SAMPLESIZE") != 20) continue;
      for (int ch = 0; ch < p->iValue(0, "CHANNELS"); ch++)
      {
        int nsat = 0, tot6 = 0, s6 = -1, mx = 0;
        for (int s = 0; s < 20; s++)
        {
          const int a = p->iValue(wf, s, ch);
          if (a >= 1023) nsat++;
          if (a > mx) mx = a;
          const int t = p->iValue(wf, s, ch, "TOT");
          if (t && s6 < 0) { tot6 = t; s6 = s; }
        }
        if (s6 < 0) continue;
        const uint64_t ts = p->lValue(wf, s6, 0);
        const int link = (ch / 72) * 2 + (ch % 72) / 36;
        int tot10 = -1, flags = -1;
        for (auto r = lines.lower_bound(ts > 20 ? ts - 20 : 0); r != lines.end() && r->first <= ts + 20; ++r)
          if (r->second.first == link) { const uint32_t w = r->second.second[ch % 36]; tot10 = (w >> 10) & 0x3ff; flags = w >> 30; break; }
        fprintf(ft, "%d %d %d %d %d %d %d\n", ch, s6, tot6, tot10, flags, nsat, mx);
      }
    }
    delete p;
    delete e;
  }
  fclose(ft); fclose(fh);
  return 0;
}
