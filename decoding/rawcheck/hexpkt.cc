// hexpkt file.evt pktid evt_skip nshow first_word last_word(-1=end) : print raw payload words (LE as stored)
#include <Event/Event.h>
#include <Event/fileEventiterator.h>
#include <Event/packet.h>
#include <Event/EventTypes.h>
#include <cstdio>
#include <vector>
int main(int argc, char **argv)
{
  int status = 0;
  fileEventiterator it(argv[1], status);
  int id = atoi(argv[2]); long skip = atol(argv[3]); int nshow = atoi(argv[4]); int w0 = atoi(argv[5]); int w1 = atoi(argv[6]);
  long n = 0; int shown = 0;
  while (Event *e = it.getNextEvent())
  {
    if (e->getEvtType() != DATAEVENT) { delete e; continue; }
    if (n++ < skip) { delete e; continue; }
    Packet *p = e->getPacket(id);
    if (p)
    {
      std::vector<int> raw(p->getLength());
      p->copyMe(raw.data(), raw.size());
      const unsigned char *c = (const unsigned char *) raw.data() + 16;
      int nb = 4 * p->getDataLength() - p->getPadding();
      int a = w0, b = (w1 < 0 ? nb / 4 : w1);
      printf("evt %ld pkt %d datalen %d pad %d", n - 1, id, p->getDataLength(), p->getPadding());
      if (id == 2001) printf(" GTT0=%lld GTT1=%lld EVTT=%u", p->lValue(0, "GROUPTRIGGERTIME"), p->lValue(1, "GROUPTRIGGERTIME"), (unsigned) p->iValue(0, "EVENTTIME"));
      printf("\n");
      if (argc > 7) { for (int i = a * 4; i < b * 4 && i < nb; i++) { printf("%02x ", c[i]); if (i % 22 == 21) printf("\n"); } printf("\n"); }
      else { const unsigned *w = (const unsigned *) c; for (int i = a; i < b; i++) printf("w%d=%08x ", i, w[i]); printf("\n"); }
      delete p;
      if (++shown >= nshow) { delete e; break; }
    }
    delete e;
  }
}
