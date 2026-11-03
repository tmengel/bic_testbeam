#!/bin/bash
# condor job: one KCU segment of run 3062 -> /sphenix/tg/tg01/jets/tmengel/eic/kcu<board>/kcu<board>_3062_seg<k>.root (+ _head.root)
#   kcu_to_root.sh <board 0|1> <segment>
EIC=/sphenix/user/tmengel/eic
source $EIC/setup.sh
BOARD=$1; SEG=$2; RUN=00003062
DATA=$EIC/data; OUT=/sphenix/tg/tg01/jets/tmengel/eic/kcu$BOARD; PID=$((12001 + BOARD))
k=$(printf "%04d" $SEG); n=$(printf "%04d" $((SEG + 1)))
NEXT=$DATA/beam_ROC_$BOARD-$RUN-$n.evt; [ -f "$NEXT" ] || NEXT=-
echo "KCU$BOARD segment $SEG: $DATA/beam_ROC_$BOARD-$RUN-$k.evt, next $NEXT, host $(hostname), start $(date)"
cd $EIC/bic_testbeam/decoding
./kcu_to_root $OUT/kcu${BOARD}_3062_seg$SEG.root $DATA/beam_ROC_$BOARD-$RUN-$k.evt $SEG $NEXT 6 $PID
rc=$?; echo "exit $rc, end $(date)"; exit $rc
