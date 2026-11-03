#!/bin/bash
# condor job: per-cell DRS4 pedestal correction of one CAEN segment -> /sphenix/tg/tg01/jets/tmengel/eic/caen/calib/
EIC=/sphenix/user/tmengel/eic
source $EIC/setup.sh
echo "start $(date) on $(hostname), $1 $2"
python3 -W ignore $EIC/bic_testbeam/calibration/caen_calib/caen_apply.py $1 $2
rc=$?; echo "exit $rc, end $(date)"; exit $rc
