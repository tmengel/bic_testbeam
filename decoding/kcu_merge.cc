// kcu_merge: the head file of segment k minus the triggers that the previous segment's main file
// already holds (its boundary-straddling trigger, complete there), written as a new file that hadd can
// chain in order: head_0 main_0 head_1' main_1 ... head_7' main_7.
//
//   kcu_merge out_head.root head_k.root main_{k-1}.root
//
// A trigger of head_k is dropped when its ts0 is not after the last sample timestamp of the previous
// main file (then its lines were all read by the previous job and are part of its last trigger).
#include <TFile.h>
#include <TTree.h>
#include <TNamed.h>
#include <TList.h>
#include <TKey.h>
#include <cstdio>
#include <cstdlib>
#include <algorithm>
int main(int argc, char **argv)
{
  if (argc < 4) { fprintf(stderr, "usage: kcu_merge out_head.root head_k.root main_prev.root\n"); return 1; }
  TFile *prev = TFile::Open(argv[3]); if (!prev || prev->IsZombie()) { fprintf(stderr, "cannot open %s\n", argv[3]); return 1; }
  TTree *ph = (TTree *) prev->Get("header"); Long64_t ts[20]; Int_t nsamp; ph->SetBranchAddress("ts", ts); ph->SetBranchAddress("nsamp", &nsamp);
  Long64_t T = 0; const Long64_t n = ph->GetEntries();
  for (Long64_t i = std::max(0LL, n - 200); i < n; i++) { ph->GetEntry(i); for (int s = 0; s < std::min(nsamp, 20); s++) T = std::max(T, ts[s]); }
  TFile *in = TFile::Open(argv[2]); if (!in || in->IsZombie()) { fprintf(stderr, "cannot open %s\n", argv[2]); return 1; }
  TTree *hh = (TTree *) in->Get("header"); Long64_t ts0; hh->SetBranchAddress("ts0", &ts0);
  const Long64_t nh = hh->GetEntries(); Long64_t nskip = 0;
  for (Long64_t i = 0; i < nh; i++) { hh->GetEntry(i); if (ts0 <= T) nskip = i + 1; else break; }   // the duplicates are the leading entries
  TFile *out = TFile::Open(argv[1], "RECREATE", "", in->GetCompressionSettings());
  const char *trees[] = {"header", "sfil1", "sfil2", "sfil3", "sfil4", "paddles", "bbcal", "sfil5n", "unused", "sfil5", "sfil5north"};
  for (const char *tn : trees)
    {
      TTree *t = (TTree *) in->Get(tn); if (!t) continue;
      out->cd(); TTree *o = t->CloneTree(0);
      for (Long64_t i = nskip; i < t->GetEntries(); i++) { t->GetEntry(i); o->Fill(); }
      o->Write();
    }
  TIter nextkey(in->GetListOfKeys());
  while (TKey *k = (TKey *) nextkey()) if (TString(k->GetClassName()) == "TNamed") { TNamed *nm = (TNamed *) k->ReadObj(); out->cd(); nm->Write(); }
  out->cd(); TNamed(Form("merge_note_%s", in->GetName()), Form("head file %s: %lld of %lld triggers dropped as already held by %s (last sample timestamp there %lld)", argv[2], nskip, nh, argv[3], T)).Write();
  out->Close();
  printf("%s: kept %lld of %lld head triggers (dropped %lld with ts0 <= %lld, the end of %s)\n", argv[1], nh - nskip, nh, nskip, T, argv[3]);
  return 0;
}
