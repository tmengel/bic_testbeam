#!/bin/bash
# condor job: template pulse fit of one KCU file of any run;  run_fit.sh <kcu0|kcu1> <tag = RUN or RUN_partPP>
# kcu0 -> kcu0_sfil_fit.C (SFIL1-4), kcu1 -> kcu1_sfil5_fit.C (SFIL5 South); output and figures in .../<board>/fit_runs/<tag>/
# optional 3rd argument: fit output of another run whose templates replace those of channels with < 100 clean pulses (low-rate runs)
EIC=/sphenix/user/tmengel/eic
source $EIC/setup.sh
B=$1; TAG=$2; TPL=$3; IN=$EIC/bic_testbeam/data/${B}_$TAG.root; OUT=/sphenix/tg/tg01/jets/tmengel/eic/$B/fit_runs/$TAG
mkdir -p $OUT; cd $OUT
echo "start $(date) on $(hostname): $B $TAG"
case $B in
  kcu0) root -l -b -q "$EIC/bic_testbeam/waveforms/kcu0_sfil_fit.C+(-1, 0, []{ FitFlags f; f.templateFile = \"$TPL\"; return f; }(), \"$IN\", \"$OUT/kcu0_sfil_fit_$TAG.root\")" ;;
  kcu1) root -l -b -q "$EIC/bic_testbeam/waveforms/kcu1_sfil5_fit.C+(-1, 0, []{ Fit1Flags f; f.templateFile = \"$TPL\"; return f; }(), \"$IN\", \"$OUT/kcu1_sfil5_fit_$TAG.root\")" ;;
esac
rc=$?; echo "exit $rc, end $(date)"; ls -la $OUT; exit $rc
