#!/bin/bash

bash /opt/sphenix/core/bin/sphenix_setup.sh  new
export ONLINE_MAIN=/eic/data/purschke/install
export LD_LIBRARY_PATH=$ONLINE_MAIN/lib:$LD_LIBRARY_PATH
PATH=$ONLINE_MAIN/bin:$PATH
export ROOT_INCLUDE_PATH=$ONLINE_MAIN/include:$ONLINE_MAIN/include/Event:$ONLINE_MAIN/include/pmonitor
echo "using the ePIC environment"

