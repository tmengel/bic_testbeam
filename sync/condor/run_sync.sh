#!/bin/bash
# condor job: decoding statistics (run_decode_stats.py) and stream synchronization (run_sync.py) of one run from the merged
# ROOT trees -> output/sync_runs/run<RUN>/
#   run_sync.sh <run> [decode]     (decode: statistics only)
EIC=/sphenix/user/tmengel/eic
source $EIC/setup.sh
echo "sync run $1, start $(date) on $(hostname)"
cd $EIC/bic_testbeam/sync
python3 -u run_decode_stats.py $1 || echo "decode stats failed"
if [ "$2" = decode ]; then rc=0; else python3 -u run_sync.py $1; rc=$?; fi; echo "exit $rc, end $(date)"; exit $rc
