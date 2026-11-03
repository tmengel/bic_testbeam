#!/bin/bash
# condor job: kcu0_sfil_fit.C over the whole KCU0 file -> /sphenix/tg/tg01/jets/tmengel/eic/kcu0/fit_full/
EIC=/sphenix/user/tmengel/eic
source $EIC/setup.sh
OUT=/sphenix/tg/tg01/jets/tmengel/eic/kcu0/fit_full
cd $OUT
echo "start $(date) on $(hostname)"
root -l -b -q "$EIC/bic_testbeam/waveforms/kcu0_sfil_fit.C+(-1, 0, FitFlags(), \"$EIC/bic_testbeam/data/kcu0_3062.root\", \"$OUT/kcu0_sfil_fit_full.root\")"
rc=$?; echo "exit $rc, end $(date)"; ls -la $OUT; exit $rc
