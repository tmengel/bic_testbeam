// wf_acc: weighted sums of pedestal-subtracted waveforms for listed triggers, stock decoders only.
//
//   wf_acc <kcu0|kcu1|caen> keys.txt out.root in1.evt [in2.evt ...]
//
// keys.txt: one line per (trigger, class): "key class weight". The key is the waveform
// timestamp lValue(wf,0,0) for kcu0/kcu1 (complete 20-sample waveforms only) and the row
// number of packet-2001 events (as in ts_dump caen) for caen. A trigger may appear in
// several classes. For every class c the output holds
//   TH2D s<c>  (channel x sample): sum of weight * (ADC - pedestal)
//   TH2D q<c>  same with weight^2 * (ADC - pedestal)^2  (for errors)
//   TH1D w<c>  bin 1: sum of weights, bin 2: sum of weights^2, bin 3: entries found
// Channels: kcu0 = the 36 kcuadc channels (SFIL1-4 South [cell][layer], North, paddles 3-6),
//           kcu1 = ASIC2 ch 144..215, caen = ch 0..15.
// Pedestal: mean of samples 1-3 (KCU), mean of samples 0-19 (CAEN); WFRAW=1 in the environment
// keeps raw ADC values instead.
#include <Event/Event.h>
#include <Event/Eventiterator.h>
#include <Event/fileEventiterator.h>
#include <Event/packet.h>
#include <Event/EventTypes.h>

#include <TFile.h>
#include <TH1D.h>
#include <TH2D.h>

#include <cstdlib>
#include <fstream>
#include <iostream>
#include <map>
#include <string>
#include <unordered_map>
#include <vector>

struct Acc { std::vector<double> s, q; double w = 0, w2 = 0, n = 0; };

int main(int argc, char **argv)
{
  if (argc < 5) { std::cerr << "usage: wf_acc <kcu0|kcu1|caen> keys.txt out.root in.evt ...\n"; return 1; }
  const std::string mode = argv[1];
  static const int kAmpCh[36] = {2, 21, 63, 43, 6, 25, 59, 37, 10, 29, 52, 0, 14, 33, 48, 4,
                                 146, 165, 207, 187, 150, 169, 203, 186, 154, 173, 196, 144, 158, 177, 192, 148,
                                 84, 86, 81, 97};
  std::vector<int> chans;
  if (mode == "kcu0") chans.assign(kAmpCh, kAmpCh + 36);
  else if (mode == "kcu1") for (int c = 144; c < 216; c++) chans.push_back(c);
  else if (mode == "caen") for (int c = 0; c < 16; c++) chans.push_back(c);
  else { std::cerr << "unknown mode " << mode << "\n"; return 1; }
  const int nch = chans.size(), nsmp = mode == "caen" ? 1024 : 20;
  const bool raw = getenv("WFRAW") != nullptr;          // WFRAW=1: no pedestal subtraction

  std::unordered_map<long long, std::vector<std::pair<int, double>>> keys;
  std::map<int, Acc> acc;
  {
    std::ifstream in(argv[2]); long long k; int c; double w;
    while (in >> k >> c >> w)
      {
        keys[k].push_back({c, w});
        if (!acc.count(c)) { acc[c].s.assign(nch * nsmp, 0); acc[c].q.assign(nch * nsmp, 0); }
      }
  }
  std::cerr << keys.size() << " keys, " << acc.size() << " classes\n";
  TFile fout(argv[3], "RECREATE");
  std::vector<double> v(nsmp);
  long long caenrow = 0;
  for (int f = 4; f < argc; f++)
    {
      int status = 0;
      fileEventiterator it(argv[f], status);
      if (status) { std::cerr << "cannot open " << argv[f] << "\n"; continue; }
      Event *e;
      while ((e = it.getNextEvent()))
        {
          if (e->getEvtType() != DATAEVENT) { delete e; continue; }
          Packet *p = e->getPacket(mode == "caen" ? 2001 : (mode == "kcu0" ? 12001 : 12002));
          if (!p) { delete e; continue; }
          // one (key, waveform index) list per event
          std::vector<std::pair<long long, int>> todo;
          if (mode == "caen") todo.push_back({caenrow++, 0});
          else
            for (int n = 0; n < p->iValue(0, "NR_WF"); n++)
              if (p->iValue(n, "SAMPLESIZE") == 20) todo.push_back({p->lValue(n, 0, 0), n});
          for (auto &[key, n] : todo)
            {
              auto hit = keys.find(key);
              if (hit == keys.end()) continue;
              const int ns = mode == "caen" ? std::min(nsmp, p->iValue(0, "SAMPLES")) : 20;
              for (int i = 0; i < nch; i++)
                {
                  const int ch = chans[i];
                  for (int s = 0; s < ns; s++) v[s] = mode == "caen" ? p->iValue(s, ch) : p->iValue(n, s, ch);
                  double ped = 0;
                  if (raw) ped = 0;
                  else if (mode == "caen") { for (int s = 0; s < 20; s++) ped += v[s]; ped /= 20; }
                  else ped = (v[1] + v[2] + v[3]) / 3.;
                  for (auto &[c, w] : hit->second)
                    {
                      Acc &a = acc[c];
                      for (int s = 0; s < ns; s++)
                        {
                          const double d = v[s] - ped;
                          a.s[i * nsmp + s] += w * d; a.q[i * nsmp + s] += w * w * d * d;
                        }
                    }
                }
              for (auto &[c, w] : hit->second) { acc[c].w += w; acc[c].w2 += w * w; acc[c].n += 1; }
            }
          delete p; delete e;
        }
      std::cerr << argv[f] << " done\n";
    }
  fout.cd();
  for (auto &[c, a] : acc)
    {
      TH2D hs(Form("s%d", c), "", nch, -0.5, nch - 0.5, nsmp, -0.5, nsmp - 0.5), hq(Form("q%d", c), "", nch, -0.5, nch - 0.5, nsmp, -0.5, nsmp - 0.5);
      TH1D hw(Form("w%d", c), "", 3, 0.5, 3.5);
      for (int i = 0; i < nch; i++)
        for (int s = 0; s < nsmp; s++) { hs.SetBinContent(i + 1, s + 1, a.s[i * nsmp + s]); hq.SetBinContent(i + 1, s + 1, a.q[i * nsmp + s]); }
      hw.SetBinContent(1, a.w); hw.SetBinContent(2, a.w2); hw.SetBinContent(3, a.n);
      hs.Write(); hq.Write(); hw.Write();
    }
  fout.Close();
  return 0;
}
