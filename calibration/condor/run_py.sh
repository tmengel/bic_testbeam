#!/bin/bash
# condor job: run a python script of the calibration thread in its own directory;  run_py.sh <script path> [args...]
EIC=/sphenix/user/tmengel/eic
source $EIC/setup.sh
S=$1; shift; cd $(dirname $S)
echo "start $(date) on $(hostname): $S $@"
python3 -W ignore $S "$@"
rc=$?; echo "exit $rc, end $(date)"; exit $rc
