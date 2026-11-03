// kcu_hists: the histograms of the LED macro PlotEventBICSingle.C, filled from a KCU stream with the
// stock decoder + ROC_pool, for all 288 channels.
//
//   kcu_hists out.root in.evt [nwaveforms=30000] [packetid=12001] [toa_mode=forward|macro] [toa_wrap_ns=6.37]
//
// TOA time of a word found in sample s:
//   forward (default, what run 3062 supports, method.md 4.9): t = 25 s + TOA_ns, no window rule;
//   macro: the LED macro's rule, t = 25 (s + (TOA_ns < wrap ? 0 : 1)) - TOA_ns, with its "magic number" wrap
//   (6.37 ns in the LED configuration; the word split in run 3062 sits at ~20.7 ns and the sign is opposite).
//
// Per channel ch (0..287), complete 20-sample waveforms:
//   hADC_ch   sample x ADC (10 bit)            hTOA_ch   sample x TOA code (10 bit)
//   hTOT_ch   sample x TOT (decompressed 12 bit, 8 counts per bin)
//   hADCcorr_ch  (as in the macro) time-aligned waveform: t = 25 x (sample + 1) - TOA_ns vs ADC, where
//                TOA_ns = fine x 0.025 + mid x 0.2 + coarse x 6.25 of the FIRST non-zero TOA word of
//                the waveform (the macro's decode_toa_value_ns); kept literal so the result can be
//                compared with the LED plots, see method.md 4.9 for what the fields really do here
//   hADCpulse_ch the same alignment but by the PULSE's TOA word: the word one sample after the largest
//                ADC step of the waveform (method.md 4.9); only waveforms with a step > 50 ADC that have
//                a TOA word there
// Per ASIC a (0..3): hADCASIC_a, hTOAASIC_a (72 channels x value), hTOTASIC_a (72 x decompressed TOT)
// Also hNwf (number of waveforms), hTc_ch (sample x Tc flag is not available through the decoder: not filled)
#include <Event/Event.h>
#include <Event/Eventiterator.h>
#include <Event/fileEventiterator.h>
#include <Event/packet.h>
#include <Event/EventTypes.h>
#include <Event/oncsSubConstants.h>
#include <ROC_pool/h2gcroc3_10G_pool.h>

#include <TFile.h>
#include <TH1D.h>
#include <TH2F.h>

#include <cstdio>
#include <cstdlib>
#include <string>
#include <vector>

static int tot_dec(int v) { return v < 512 ? v : 8 * (v - 512) + 4; }
static double toa_ns(int v) { return (v & 7) * 0.025 + ((v >> 3) & 31) * 0.2 + ((v >> 8) & 3) * 6.25; }
static bool toa_macro = false; static double toa_wrap = 6.37;
// time (ns from the start of sample 0) of the TOA word found in sample s
static double toa_time(int s, int code) { const double t = toa_ns(code); return toa_macro ? (t < toa_wrap ? s : s + 1) * 25.0 - t : 25.0 * s + t; }

int main(int argc, char **argv)
{
  if (argc < 3) { fprintf(stderr, "usage: kcu_hists out.root in.evt [nwaveforms] [packetid]\n"); return 1; }
  const long nmax = argc > 3 ? atol(argv[3]) : 30000;
  const int pid = argc > 4 ? atoi(argv[4]) : 12001;
  if (argc > 5) toa_macro = std::string(argv[5]) == "macro";
  if (argc > 6) toa_wrap = atof(argv[6]);
  const int NCH = 288;
  TFile fout(argv[1], "RECREATE");
  std::vector<TH2F *> hADC(NCH), hTOA(NCH), hTOT(NCH), hCorr(NCH), hPulse(NCH);
  for (int c = 0; c < NCH; c++)
    {
      hPulse[c] = new TH2F(Form("hADCpulse_%d", c), Form("ch %d, time-aligned by the pulse TOA word;time [ns];ADC", c), 400, -100, 900, 1024, -0.5, 1023.5);
      hADC[c] = new TH2F(Form("hADC_%d", c), Form("ch %d ADC;sample;ADC", c), 20, -0.5, 19.5, 1024, -0.5, 1023.5);
      hTOA[c] = new TH2F(Form("hTOA_%d", c), Form("ch %d TOA;sample;TOA code", c), 20, -0.5, 19.5, 1024, -0.5, 1023.5);
      hTOT[c] = new TH2F(Form("hTOT_%d", c), Form("ch %d TOT;sample;TOT (decompressed)", c), 20, -0.5, 19.5, 512, -0.5, 4095.5);
      hCorr[c] = new TH2F(Form("hADCcorr_%d", c), Form("ch %d, time-aligned by the first TOA word;time [ns];ADC", c), 400, -100, 900, 1024, -0.5, 1023.5);
    }
  TH2F *hADCA[4], *hTOAA[4], *hTOTA[4];
  for (int a = 0; a < 4; a++)
    {
      hADCA[a] = new TH2F(Form("hADCASIC_%d", a), Form("ASIC %d ADC;channel;ADC", a), 72, -0.5, 71.5, 1024, -0.5, 1023.5);
      hTOAA[a] = new TH2F(Form("hTOAASIC_%d", a), Form("ASIC %d TOA;channel;TOA code", a), 72, -0.5, 71.5, 1024, -0.5, 1023.5);
      hTOTA[a] = new TH2F(Form("hTOTASIC_%d", a), Form("ASIC %d TOT;channel;TOT (decompressed)", a), 72, -0.5, 71.5, 512, -0.5, 4095.5);
    }
  TH1D hNwf("hNwf", "waveforms", 1, 0, 1);

  h2gcroc3_10G_pool *pool = nullptr;
  long nwf = 0;
  auto process = [&](h2gcroc3_10G_pool *p)
    {
      const int n = p->iValue(0, "NR_WF");
      for (int w = 0; w < n && nwf < nmax; w++)
        {
          if (p->iValue(w, "SAMPLESIZE") != 20) continue;
          nwf++;
          for (int c = 0; c < NCH; c++)
            {
              const int a = c / 72;
              int adc[20], toa[20], ftoa = -1, fs = -1;   // ftoa: first TOA word of the waveform (as the macro)
              for (int s = 0; s < 20; s++)
                {
                  adc[s] = p->iValue(w, s, c); toa[s] = p->iValue(w, s, c, "TOA");
                  const int tot = tot_dec(p->iValue(w, s, c, "TOT"));
                  hADC[c]->Fill(s, adc[s]); hADCA[a]->Fill(c % 72, adc[s]);
                  if (toa[s]) { hTOA[c]->Fill(s, toa[s]); hTOAA[a]->Fill(c % 72, toa[s]); if (ftoa < 0) { ftoa = toa[s]; fs = s; } }
                  if (tot) { hTOT[c]->Fill(s, tot); hTOTA[a]->Fill(c % 72, tot); }
                }
              if (ftoa > 0)
                {
                  const double t0 = toa_time(fs, ftoa) - 25.0 * (fs + 1);   // macro: samples relative to the TOA time, with the wrap rule
                  for (int s = 0; s < 20; s++) hCorr[c]->Fill(25.0 * (s + 1) - 25.0 * (fs + 1) - t0, adc[s]);
                }
              // the pulse's TOA word: one sample after the largest ADC step (positive or negative pulses)
              int r = 0;
              for (int s = 1; s < 19; s++) if (abs(adc[s + 1] - adc[s]) > abs(adc[r + 1] - adc[r])) r = s;
              const int sw = toa[r + 1] > 0 ? r + 1 : (toa[r] > 0 ? r : -1);   // the pulse's TOA word: one sample after the largest step, or in it (late in the window)
              if (abs(adc[r + 1] - adc[r]) > 50 && sw >= 0)
                {
                  const double t0 = toa_time(sw, toa[sw]);   // ns from the start of sample 0
                  for (int s = 0; s < 20; s++) hPulse[c]->Fill(25.0 * s - t0 + 200.0, adc[s]);   // +200 ns so the pulse sits mid-range as before
                }
            }
        }
    };

  int status = 0;
  fileEventiterator it(argv[2], status);
  if (status) { fprintf(stderr, "cannot open %s\n", argv[2]); return 1; }
  while (Event *e = it.getNextEvent())
    {
      if (nwf >= nmax) { delete e; break; }
      if (e->getEvtType() != DATAEVENT) { delete e; continue; }
      if (Packet *p = e->getPacket(pid))
        {
          if (!pool)
            {
              if (p->getHitFormat() != IDH2GCROC3_10G) { fprintf(stderr, "packet %d is not IDH2GCROC3_10G\n", pid); delete p; delete e; return 1; }
              pool = new h2gcroc3_10G_pool();
            }
          pool->addPacket(p); delete p;
          while (nwf < nmax && pool->depth_ok()) { process(pool); pool->next(); }
        }
      delete e;
    }
  if (pool && nwf < nmax) { pool->drain(); while (nwf < nmax && pool->depth_ok()) { process(pool); pool->next(); } }
  hNwf.SetBinContent(1, nwf);
  fout.Write();
  printf("%ld complete waveforms, wrote %s\n", nwf, argv[1]);
  return 0;
}
