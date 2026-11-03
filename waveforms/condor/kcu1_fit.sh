#!/bin/bash
# condor job: kcu1_sfil5_fit.C over the whole KCU1 file -> /sphenix/tg/tg01/jets/tmengel/eic/kcu1/fit_full/
EIC=/sphenix/user/tmengel/eic
source $EIC/setup.sh
OUT=/sphenix/tg/tg01/jets/tmengel/eic/kcu1/fit_full
mkdir -p $OUT; cd $OUT
echo "start $(date) on $(hostname)"
root -l -b -q "$EIC/bic_testbeam/waveforms/kcu1_sfil5_fit.C+(-1, 0, Fit1Flags(), \"$EIC/bic_testbeam/data/kcu1_3062.root\", \"$OUT/kcu1_sfil5_fit_full.root\")"
rc=$?; echo "exit $rc, end $(date)"; ls -la $OUT; exit $rc
