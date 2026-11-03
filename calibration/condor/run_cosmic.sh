#!/bin/bash
# condor job: cosmic pulse finding on one KCU0 cosmic file;  run_cosmic.sh <run> <start> <nmax>   (nmax -1 = to the end)
# -> /sphenix/tg/tg01/jets/tmengel/eic/kcu0/cosmics/cosmic_<run>_<start>.root ; hadd the chunks, then calibration/cosmics/cosmic_mip_fit.py
EIC=/sphenix/user/tmengel/eic
source $EIC/setup.sh
RUN=$1; START=$2; NMAX=${3:--1}; OUT=/sphenix/tg/tg01/jets/tmengel/eic/kcu0/cosmics; mkdir -p $OUT; cd $OUT
echo "start $(date) on $(hostname): run $RUN from $START, $NMAX triggers"
root -l -b -q "$EIC/bic_testbeam/calibration/cosmics/cosmic_mip.C+(\"$EIC/bic_testbeam/data/kcu0_$RUN.root\", \"$OUT/cosmic_${RUN}_$START.root\", $START, $NMAX)"
rc=$?; echo "exit $rc, end $(date)"; exit $rc
