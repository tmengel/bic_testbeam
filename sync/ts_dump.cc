// ts_dump: one timestamp per trigger (or per hit) from any BIC stream, using
// only the stock decoders of the online_distribution Event library.
//
//   ts_dump <kcu|kcuadc|kcu1adc|caen|caenadc|astro> out.root in1.evt [in2.evt ...]
//
// Every row carries the RTCLOCK host time (ns) of its RCDAQ event, so all
// streams can be put on one coarse timeline before matching in device ticks.
//   kcu   packet 12001/12002  tree "kcu":   evt rt ts nsamp
//         ts = lValue(wf,0,0), FPGA timestamp of sample 0 (6.25 ns)
//   kcuwf   first 200k complete KCU0 waveforms: wf[8][20] raw samples of paddles 3-6 (ch 84,86,81,97)
//           then SFIL1-4 South beam row (ch 10,29,52,0)
//   kcutot  first 100k complete KCU0 waveforms: wf[36][20] ADC and tot[36][20] TOT for the kcuadc channels
//   kcuadc  as kcu, complete 20-sample waveforms only, plus amp[36]: max - pedestal
//         (pedestal = mean of samples 1-3, as roc1) for SFIL1-4 South [cell][layer],
//         SFIL1-4 North [cell][layer], then paddles 3,4,5,6 (8/13 hookup)
//         paddles are pedestal - min of samples 8-16 instead (a hit saturates at ADC 0, so ~80)
//         win[36]: mean of samples 6-9 (the pulse peak) - pedestal: zero mean on noise, unlike max - ped
//         cm[4]:   the same quantity, median over the unconnected channels of each ASIC A0..A3
//                  (event-wide baseline shift; subtract the channel's ASIC entry: A0 for South, A2 North)
//   kcu1adc  KCU1 (packet 12002), complete waveforms: ts, a1[72] = max - pedestal for
//         ASIC2 ch 144..215 (SFIL5 South individual SiPMs)
//   caenadc  as caen plus, for ch 0..15 (groups 0,1): neg[16] = ped - min, pos[16] = max - ped,
//         ped = mean of samples 0..19; qint[16] = sum over samples 20..159 of (ped - ADC): the
//         pedestal-subtracted pulse integral (negative pulses), zero mean on noise
//   caen  packet 2001         tree "caen":  evt rt evnr evttime gtt
//         gtt = lValue(0,"GROUPTRIGGERTIME") (8.533 ns, extended by decoder)
//   astro packet 13001        tree "astro": evt rt fpga_ts ts tot id loc iscol layer
//         fpga_ts = lValue(hit), 40 MHz
#include <Event/Event.h>
#include <Event/Eventiterator.h>
#include <Event/fileEventiterator.h>
#include <Event/packet.h>
#include <Event/EventTypes.h>

#include <TFile.h>
#include <TTree.h>

#include <algorithm>
#include <cstring>
#include <iostream>
#include <string>
#include <vector>

// median over the unconnected, non-railed channels of one ASIC of (mean of 4 samples from s0) - mean(samples 1-3)
static float common_mode(Packet *p, int n, int ch0, const std::vector<bool> &used, int s0)
{
  std::vector<float> v;
  for (int ch = ch0; ch < ch0 + 72; ch++)
    {
      if (used[ch]) continue;
      const float ped = (p->iValue(n, 1, ch) + p->iValue(n, 2, ch) + p->iValue(n, 3, ch)) / 3.f;
      if (ped > 1000 || ped < 5) continue;                       // railed or dead
      v.push_back((p->iValue(n, s0, ch) + p->iValue(n, s0 + 1, ch) + p->iValue(n, s0 + 2, ch) + p->iValue(n, s0 + 3, ch)) / 4.f - ped);
    }
  if (v.empty()) return 0;
  std::nth_element(v.begin(), v.begin() + v.size() / 2, v.end());
  return v[v.size() / 2];
}

int main(int argc, char **argv)
{
  if (argc < 4)
    {
      std::cerr << "usage: ts_dump <kcu|caen|astro> out.root in.evt [...]\n";
      return 1;
    }
  const std::string mode = argv[1];
  TFile fout(argv[2], "RECREATE");
  TTree t(mode.c_str(), mode.c_str());

  int evt = 0, nsamp = 0, evnr = 0, tot = 0, id = 0, loc = 0, iscol = 0, layer = 0, ats = 0;
  Long64_t rt = 0, ts = 0, gtt = 0;
  UInt_t evttime = 0;
  static const int kAmpCh[36] = {2, 21, 63, 43, 6, 25, 59, 37, 10, 29, 52, 0, 14, 33, 48, 4,
                                 146, 165, 207, 187, 150, 169, 203, 186, 154, 173, 196, 144, 158, 177, 192, 148,
                                 84, 86, 81, 97};
  Float_t amp[36], win[36], cm[4], w1[72], cm1, qint[16];
  std::vector<bool> used(288, false), used1(288, false);      // connected channels, excluded from the common mode
  for (int k = 0; k < 36; k++) used[kAmpCh[k]] = true;
  for (int ch = 144; ch < 216; ch++) used1[ch] = ((ch - 144) % 9 != 8);   // KCU1 ASIC2: every 9th channel is unconnected
  Short_t a1[72], neg[16], pos[16], wf[8][20], wf36[36][20], tot36[36][20];
  static const int kWfCh[8] = {84, 86, 81, 97, 10, 29, 52, 0};
  Long64_t nfill = 0;
  t.Branch("evt", &evt);
  t.Branch("rt", &rt);
  if (mode == "kcu") { t.Branch("ts", &ts); t.Branch("nsamp", &nsamp); }
  else if (mode == "kcuadc") { t.Branch("ts", &ts); t.Branch("amp", amp, "amp[36]/F"); t.Branch("win", win, "win[36]/F"); t.Branch("cm", cm, "cm[4]/F"); }
  else if (mode == "kcuwf") { t.Branch("ts", &ts); t.Branch("wf", wf, "wf[8][20]/S"); }
  else if (mode == "kcutot") { t.Branch("ts", &ts); t.Branch("wf", wf36, "wf[36][20]/S"); t.Branch("tot", tot36, "tot[36][20]/S"); }
  else if (mode == "kcu1adc") { t.Branch("ts", &ts); t.Branch("a1", a1, "a1[72]/S"); t.Branch("w1", w1, "w1[72]/F"); t.Branch("cm1", &cm1); }
  else if (mode == "caen") { t.Branch("evnr", &evnr); t.Branch("evttime", &evttime); t.Branch("gtt", &gtt); }
  else if (mode == "caenadc") { t.Branch("evnr", &evnr); t.Branch("gtt", &gtt); t.Branch("neg", neg, "neg[16]/S"); t.Branch("pos", pos, "pos[16]/S"); t.Branch("qint", qint, "qint[16]/F"); }
  else if (mode == "astro")
    {
      t.Branch("fpga_ts", &ts); t.Branch("ts", &ats); t.Branch("tot", &tot); t.Branch("id", &id);
      t.Branch("loc", &loc); t.Branch("iscol", &iscol); t.Branch("layer", &layer);
    }
  else { std::cerr << "unknown mode " << mode << "\n"; return 1; }

  for (int f = 3; f < argc; f++)
    {
      int status = 0;
      fileEventiterator it(argv[f], status);
      if (status) { std::cerr << "cannot open " << argv[f] << "\n"; continue; }
      Event *e;
      while ((e = it.getNextEvent()))
        {
          if (e->getEvtType() != DATAEVENT) { delete e; continue; }
          evt = e->getEvtSequence();
          rt = 0;
          for (int id_rt = 1010; id_rt <= 1013; id_rt++)
            if (Packet *pr = e->getPacket(id_rt)) { rt = pr->lValue(0, "TOTAL"); delete pr; break; }

          if (mode == "kcutot")
            {
              if (nfill >= 100000) { delete e; break; }
              if (Packet *p = e->getPacket(12001))
                {
                  const int nwf = p->iValue(0, "NR_WF");
                  for (int n = 0; n < nwf; n++)
                    {
                      if (p->iValue(n, "SAMPLESIZE") != 20) continue;
                      ts = p->lValue(n, 0, 0);
                      for (int k = 0; k < 36; k++)
                        for (int smp = 0; smp < 20; smp++)
                          {
                            wf36[k][smp] = (Short_t) p->iValue(n, smp, kAmpCh[k]);
                            tot36[k][smp] = (Short_t) p->iValue(n, smp, kAmpCh[k], "TOT");
                          }
                      t.Fill(); nfill++;
                    }
                  delete p;
                }
            }
          else if (mode == "kcuwf")
            {
              if (nfill >= 200000) { delete e; break; }
              if (Packet *p = e->getPacket(12001))
                {
                  const int nwf = p->iValue(0, "NR_WF");
                  for (int n = 0; n < nwf; n++)
                    {
                      if (p->iValue(n, "SAMPLESIZE") != 20) continue;
                      ts = p->lValue(n, 0, 0);
                      for (int k = 0; k < 8; k++)
                        for (int smp = 0; smp < 20; smp++) wf[k][smp] = (Short_t) p->iValue(n, smp, kWfCh[k]);
                      t.Fill(); nfill++;
                    }
                  delete p;
                }
            }
          else if (mode == "kcuadc")
            {
              Packet *p = e->getPacket(12001);
              if (!p) p = e->getPacket(12002);
              if (p)
                {
                  const int nwf = p->iValue(0, "NR_WF");
                  for (int n = 0; n < nwf; n++)
                    {
                      if (p->iValue(n, "SAMPLESIZE") != 20) continue;
                      ts = p->lValue(n, 0, 0);
                      for (int k = 0; k < 36; k++)
                        {
                          const int ch = kAmpCh[k];
                          const float ped = (p->iValue(n, 1, ch) + p->iValue(n, 2, ch) + p->iValue(n, 3, ch)) / 3.f;
                          if (k >= 32)    // paddles: negative PMT pulse saturates the ADC at 0 -> dip depth
                            {
                              int mn = 1 << 30;
                              for (int smp = 8; smp <= 16; smp++) mn = std::min(mn, p->iValue(n, smp, ch));
                              amp[k] = ped - mn;
                              continue;
                            }
                          int mx = 0;
                          for (int smp = 0; smp < 20; smp++) mx = std::max(mx, p->iValue(n, smp, ch));
                          amp[k] = mx - ped;
                          win[k] = (p->iValue(n, 6, ch) + p->iValue(n, 7, ch) + p->iValue(n, 8, ch) + p->iValue(n, 9, ch)) / 4.f - ped;
                        }
                      for (int a = 0; a < 4; a++) cm[a] = common_mode(p, n, a * 72, used, 6);
                      t.Fill();
                    }
                  delete p;
                }
            }
          else if (mode == "kcu1adc")
            {
              if (Packet *p = e->getPacket(12002))
                {
                  const int nwf = p->iValue(0, "NR_WF");
                  for (int n = 0; n < nwf; n++)
                    {
                      if (p->iValue(n, "SAMPLESIZE") != 20) continue;
                      ts = p->lValue(n, 0, 0);
                      for (int k = 0; k < 72; k++)
                        {
                          const int ch = 144 + k;
                          const float ped = (p->iValue(n, 1, ch) + p->iValue(n, 2, ch) + p->iValue(n, 3, ch)) / 3.f;
                          int mx = 0;
                          for (int smp = 0; smp < 20; smp++) mx = std::max(mx, p->iValue(n, smp, ch));
                          a1[k] = (Short_t) (mx - ped);
                          w1[k] = (p->iValue(n, 5, ch) + p->iValue(n, 6, ch) + p->iValue(n, 7, ch) + p->iValue(n, 8, ch)) / 4.f - ped;
                        }
                      cm1 = common_mode(p, n, 144, used1, 5);
                      t.Fill();
                    }
                  delete p;
                }
            }
          else if (mode == "caenadc")
            {
              if (Packet *p = e->getPacket(2001))
                {
                  evnr = p->iValue(0, "EVNR");
                  gtt = p->lValue(0, "GROUPTRIGGERTIME");
                  const int ns = p->iValue(0, "SAMPLES");
                  for (int ch = 0; ch < 16; ch++)
                    {
                      double ped = 0;
                      for (int smp = 0; smp < 20; smp++) ped += p->iValue(smp, ch);
                      ped /= 20.;
                      int mn = 1 << 30, mx = -(1 << 30);
                      for (int smp = 0; smp < ns; smp++) { const int v = p->iValue(smp, ch); mn = std::min(mn, v); mx = std::max(mx, v); }
                      neg[ch] = (Short_t) (ped - mn); pos[ch] = (Short_t) (mx - ped);
                      qint[ch] = 0;
                      for (int smp = 20; smp < std::min(ns, 160); smp++) qint[ch] += ped - p->iValue(smp, ch);
                    }
                  t.Fill();
                  delete p;
                }
            }
          else if (mode == "kcu")
            {
              Packet *p = e->getPacket(12001);
              if (!p) p = e->getPacket(12002);
              if (p)
                {
                  const int nwf = p->iValue(0, "NR_WF");
                  for (int n = 0; n < nwf; n++)
                    {
                      nsamp = p->iValue(n, "SAMPLESIZE");
                      ts = p->lValue(n, 0, 0);
                      t.Fill();
                    }
                  delete p;
                }
            }
          else if (mode == "caen")
            {
              if (Packet *p = e->getPacket(2001))
                {
                  evnr = p->iValue(0, "EVNR");
                  evttime = (UInt_t) p->iValue(0, "EVENTTIME");
                  gtt = p->lValue(0, "GROUPTRIGGERTIME");
                  t.Fill();
                  delete p;
                }
            }
          else
            {
              if (Packet *p = e->getPacket(13001))
                {
                  const int nu = p->iValue(0, "NR_UNITS");
                  for (int u = 0; u < nu; u++)
                    {
                      ts = p->lValue(u);
                      ats = p->iValue(u, "TIMESTAMP"); tot = p->iValue(u, "TOT");
                      id = p->iValue(u, "ID"); loc = p->iValue(u, "LOCATION");
                      iscol = p->iValue(u, "ISCOL"); layer = p->iValue(u, "LAYER");
                      t.Fill();
                    }
                  delete p;
                }
            }
          delete e;
        }
      std::cout << argv[f] << " done, " << t.GetEntries() << " rows so far\n";
    }
  t.Write();
  return 0;
}
