#!/bin/bash
# condor job: cosmic muons through the v4 model (macros/cosmic.mac);  cosmic.sh <nevents> <seed> [subdir]  -> /sphenix/tg/tg01/jets/tmengel/eic/sim/cosmic[/subdir]/bic_cosmic_s<seed>.root
# environment EXTRA = ";"-separated macro commands inserted before /bic/output (e.g. a smaller start plane)
EIC=/sphenix/user/tmengel/eic
source $EIC/setup.sh
N=$1; SEED=$2; OUT=/sphenix/tg/tg01/jets/tmengel/eic/sim/cosmic${3:+/$3}; mkdir -p $OUT
X_EXTRA=""; if [ -n "$EXTRA" ]; then X_EXTRA=$(echo "$EXTRA" | tr ";" "\n" | sed -e "s#\$#\\\\n#" | tr -d "\n"); fi
cd $EIC/bic_testbeam/sim/build
M=$(mktemp /tmp/bic_cosmic_XXXX.mac)
sed -e "s#/bic/output bic_cosmic.root#${X_EXTRA}/bic/output $OUT/bic_cosmic_s${SEED}.root#" -e "s#/run/beamOn 1000#/run/beamOn $N#" cosmic.mac > $M
echo "cosmic: $N muons, seed $SEED, host $(hostname), start $(date)"
./bic_sim $M $SEED; rc=$?; rm -f $M; echo "exit $rc, end $(date)"; exit $rc
