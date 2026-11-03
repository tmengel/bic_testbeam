#include <Event/Event.h>
#include <Event/fileEventiterator.h>
#include <Event/packet.h>
#include <Event/EventTypes.h>
#include <cstdio>
#include <vector>
int main(int argc, char **argv)
{
  int status = 0; fileEventiterator it(argv[1], status); long nmax = atol(argv[2]);
  FILE *f = fopen(argv[3], "wb"); long n = 0;
  while (Event *e = it.getNextEvent())
  {
    if (e->getEvtType() != DATAEVENT) { delete e; continue; }
    if (n++ >= nmax) { delete e; break; }
    if (Packet *p = e->getPacket(atoi(argv[4])))
    {
      std::vector<int> raw(p->getLength()); p->copyMe(raw.data(), raw.size());
      unsigned nb = 4 * p->getDataLength() - p->getPadding();
      fwrite(&nb, 4, 1, f); fwrite((char *) raw.data() + 16, 1, nb, f);
      delete p;
    }
    delete e;
  }
  fclose(f);
}
