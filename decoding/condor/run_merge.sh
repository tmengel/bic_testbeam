#!/bin/bash
# condor job, after every segment of the run is converted: one file per stream and run,
#   /sphenix/tg/tg01/jets/tmengel/eic/<stream>/run<RUN>/<stream>_<RUN>.root, symlinked as bic_testbeam/data/<stream>_<RUN>.root
#   run_merge.sh <stream kcu0|kcu1|caen|astropix> <run> [part]
# With a part number p, only segments 10p..10p+9 are merged, into <stream>_<RUN>_partPP.root (run 3064: one file would be
# ~350 GB, past ROOT's 100 GB tree limit); a part starts on a de-duplicated head file, so parts never share a trigger.
# KCU: drop from each head file the triggers the previous segment already holds (kcu_merge), then chain head_0 main_0 head_1' ...
# CAEN: one trigger per RCDAQ event, nothing straddles: chain the segments. AstroPix: one segment, renamed.
EIC=/sphenix/user/tmengel/eic
source $EIC/setup.sh
STREAM=$1; RUN=$2; PART=$3; OUT=/sphenix/tg/tg01/jets/tmengel/eic/$STREAM/run$RUN; B=$OUT/${STREAM}_${RUN}
cd $EIC/bic_testbeam/decoding
NSEG=$(ls $OUT/${STREAM}_${RUN}_seg*.root 2>/dev/null | grep -v -E "_head|_dedup" | wc -l)
echo "merge $STREAM run $RUN: $NSEG segments, start $(date) on $(hostname)"
[ $NSEG -gt 0 ] || { echo "no segments"; exit 1; }
K0=0; K1=$NSEG; T=$B
if [ -n "$PART" ]; then
  K0=$((10 * PART)); K1=$((K0 + 10 < NSEG ? K0 + 10 : NSEG)); T=${B}_part$(printf %02d $PART)
  [ $K0 -lt $NSEG ] || { echo "part $PART beyond $NSEG segments"; exit 1; }
  echo "part $PART: segments $K0..$((K1 - 1))"
fi
case $STREAM in
  kcu0|kcu1)
    LIST=""
    for ((k = K0; k < K1; k++)); do
      if [ $k -eq 0 ]; then LIST="${B}_seg0_head.root ${B}_seg0.root"; continue; fi
      ./kcu_merge ${B}_seg${k}_head_dedup.root ${B}_seg${k}_head.root ${B}_seg$((k - 1)).root || exit 1
      LIST="$LIST ${B}_seg${k}_head_dedup.root ${B}_seg${k}.root"
    done
    hadd -ff -f $T.root $LIST || exit 1
    python3 kcu0_root_check.py $T.root > ${T}_check.txt 2>&1; cat ${T}_check.txt ;;
  caen)
    LIST=""; for ((k = K0; k < K1; k++)); do LIST="$LIST ${B}_seg$k.root"; done
    hadd -ff -f $T.root $LIST || exit 1 ;;
  astropix) cp ${B}_seg0.root $T.root || exit 1 ;;
esac
ln -sfn $T.root $EIC/bic_testbeam/data/$(basename $T).root
python3 -c "
import uproot; f = uproot.open('$T.root')
print({k: f[k].num_entries for k in sorted({n.split(';')[0] for n, c in f.classnames().items() if c == 'TTree'})})"
echo "merge end $(date)"
