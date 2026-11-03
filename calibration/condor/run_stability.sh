#!/bin/bash
# condor job: response vs time of one fit-tree file;  run_stability.sh <kcu0|kcu1> <tag = RUN or RUN_partPP> <tmax s>
# -> /sphenix/tg/tg01/jets/tmengel/eic/<board>/stability/stability_<tag>.root  (hadd the parts, then calibration/stability/gain_stability_figs.py)
EIC=/sphenix/user/tmengel/eic
source $EIC/setup.sh
B=$1; TAG=$2; TMAX=${3:-45000}; T=/sphenix/tg/tg01/jets/tmengel/eic; OUT=$T/$B/stability; mkdir -p $OUT; cd $OUT
echo "start $(date) on $(hostname): $B $TAG"
case $B in
  kcu0) root -l -b -q "$EIC/bic_testbeam/calibration/stability/gain_stability.C+(\"kcu0\", \"$T/kcu0/fit_runs/$TAG/kcu0_sfil_fit_$TAG.root\", \"$EIC/bic_testbeam/data/kcu0_$TAG.root\", \"$OUT/stability_$TAG.root\", 0, $TMAX)" ;;
  kcu1) root -l -b -q "$EIC/bic_testbeam/calibration/stability/gain_stability.C+(\"kcu1\", \"$T/kcu1/fit_runs/$TAG/kcu1_sfil5_fit_$TAG.root\", \"\", \"$OUT/stability_$TAG.root\", 0, $TMAX)" ;;
esac
rc=$?; echo "exit $rc, end $(date)"; exit $rc
