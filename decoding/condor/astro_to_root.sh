#!/bin/bash
# condor job: the AstroPix file -> /sphenix/tg/tg01/jets/tmengel/eic/astropix/astropix_3062.root
EIC=/sphenix/user/tmengel/eic
source $EIC/setup.sh
OUT=/sphenix/tg/tg01/jets/tmengel/eic/astropix
cd $EIC/bic_testbeam/decoding
./astro_to_root $OUT/astropix_3062.root $EIC/data/beam_AstroPix-00003062-0000.evt && ln -sfn $OUT/astropix_3062.root $EIC/bic_testbeam/data/astropix_3062.root
