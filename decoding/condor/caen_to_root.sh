#!/bin/bash
# condor job: one CAEN segment -> /sphenix/tg/tg01/jets/tmengel/eic/caen/caen_3062_seg<k>.root;  caen_to_root.sh <segment>
EIC=/sphenix/user/tmengel/eic
source $EIC/setup.sh
SEG=$1; k=$(printf "%04d" $SEG); OUT=/sphenix/tg/tg01/jets/tmengel/eic/caen
echo "CAEN segment $SEG, host $(hostname), start $(date)"
cd $EIC/bic_testbeam/decoding
./caen_to_root $OUT/caen_3062_seg$SEG.root $EIC/data/beam_CAEN-00003062-$k.evt $SEG
rc=$?; echo "exit $rc, end $(date)"; exit $rc
