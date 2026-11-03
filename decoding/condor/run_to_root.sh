#!/bin/bash
# condor job: one segment of one stream of any run -> /sphenix/tg/tg01/jets/tmengel/eic/<stream>/run<RUN>/<stream>_<RUN>_seg<k>.root
#   run_to_root.sh <stream kcu0|kcu1|caen|astropix> <run> <segment> [raw directory]
# The run type prefix (beam_/cosmics_) is taken from the file found; the KCU hookup epoch follows the run number
# (kcu_to_root: runs >= 3085 have the SFIL5 North summed boards on KCU0, tree sfil5n).
EIC=/sphenix/user/tmengel/eic
source $EIC/setup.sh
STREAM=$1; RUN=$2; SEG=$3; RAW=${4:-/gpfs01/eic/data1/BIC_Beamtest}
case $STREAM in
  kcu0) TAG=ROC_0; PID=12001 ;;
  kcu1) TAG=ROC_1; PID=12002 ;;
  caen) TAG=CAEN ;;
  astropix) TAG=AstroPix ;;
  *) echo "unknown stream $STREAM"; exit 1 ;;
esac
R8=$(printf "%08d" $RUN); k=$(printf "%04d" $SEG); n=$(printf "%04d" $((SEG + 1)))
IN=$(ls $RAW/*_$TAG-$R8-$k.evt 2>/dev/null | head -1)
[ -n "$IN" ] || { echo "no input for $STREAM run $RUN segment $SEG in $RAW"; exit 1; }
OUT=/sphenix/tg/tg01/jets/tmengel/eic/$STREAM/run$RUN; mkdir -p $OUT; B=$OUT/${STREAM}_${RUN}_seg$SEG
echo "$STREAM run $RUN segment $SEG: $IN -> $B.root, host $(hostname), start $(date)"
cd $EIC/bic_testbeam/decoding
case $STREAM in
  kcu0|kcu1)
    NEXT=$(ls $RAW/*_$TAG-$R8-$n.evt 2>/dev/null | head -1); [ -n "$NEXT" ] || NEXT=-
    ./kcu_to_root $B.root $IN $SEG $NEXT 6 $PID 0 $RUN ;;
  caen) ./caen_to_root $B.root $IN $SEG ;;
  astropix) ./astro_to_root $B.root $IN ;;
esac
rc=$?; echo "exit $rc, end $(date)"; exit $rc
