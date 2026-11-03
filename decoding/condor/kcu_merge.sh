#!/bin/bash
# condor job: de-duplicate the segment boundaries and chain the segments of one KCU in order
#   kcu_merge.sh <board 0|1> <number of segments>
EIC=/sphenix/user/tmengel/eic
source $EIC/setup.sh
BOARD=$1; NSEG=$2; OUT=/sphenix/tg/tg01/jets/tmengel/eic/kcu$BOARD; B=$OUT/kcu${BOARD}_3062
cd $EIC/bic_testbeam/decoding
echo "merge KCU$BOARD start $(date) on $(hostname)"
LIST="${B}_seg0_head.root ${B}_seg0.root"
for ((k = 1; k < NSEG; k++)); do
  ./kcu_merge ${B}_seg${k}_head_dedup.root ${B}_seg${k}_head.root ${B}_seg$((k - 1)).root || exit 1
  LIST="$LIST ${B}_seg${k}_head_dedup.root ${B}_seg${k}.root"
done
echo "hadd order: $LIST"
hadd -ff -f $B.root $LIST || exit 1
ln -sfn $B.root $EIC/bic_testbeam/data/kcu${BOARD}_3062.root
python3 kcu0_root_check.py $B.root > ${B}_check.txt 2>&1; cat ${B}_check.txt
echo "merge end $(date)"
