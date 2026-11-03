#!/bin/bash
# condor job: bic_sim energy point;  scan.sh <energy MeV> <nevents> <seed> [beam x mm = 30 (under P6); 80 = under P5] [Baby BCAL y offset mm = 0] [beam y mm = -5] [output subdir, e.g. v4]
# environment EXTRA = ";"-separated macro commands inserted before /run/initialize (e.g. the v4b light model)
EIC=/sphenix/user/tmengel/eic
source $EIC/setup.sh
E=$1; N=$2; SEED=$3; X=${4:-30}; Y=${5:-0}; YB=${6:--5}; OUT=/sphenix/tg/tg01/jets/tmengel/eic/sim${7:+/$7}; mkdir -p $OUT
TAG=""; [ "$X" != "30" ] && TAG="_x$X"; [ "$Y" != "0" ] && TAG="${TAG}_bby$Y"; [ "$YB" != "-5" ] && TAG="${TAG}_yb$YB"
cd $EIC/bic_testbeam/sim/build
M=$(mktemp /tmp/bic_sim_XXXX.mac)
X_EXTRA=""; if [ -n "$EXTRA" ]; then X_EXTRA=$(echo "$EXTRA" | tr ";" "\n" | sed -e "s#\$#\\\\n#" | tr -d "\n"); fi
sed -e "s#/bic/gun/energy 4000#/bic/gun/energy $E#" -e "s#/bic/gun/x 30#/bic/gun/x $X#" -e "s#/bic/output bic_sim.root#/bic/output $OUT/bic_sim_e${E}${TAG}_s${SEED}.root#" -e "s#/run/beamOn 100#/run/beamOn $N#" -e "s#/run/initialize#/bic/geom/bbcalYOffset $Y mm\n${X_EXTRA}/run/initialize#" -e "s#/bic/gun/y -5#/bic/gun/y $YB#" run.mac > $M
echo "energy $E MeV, $N events, seed $SEED, host $(hostname), start $(date)"; cat $M
./bic_sim $M $SEED; rc=$?; rm -f $M; echo "exit $rc, end $(date)"; exit $rc
