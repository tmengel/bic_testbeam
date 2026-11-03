#!/bin/bash
# condor job: chain the CAEN segments in order (one trigger per RCDAQ event, nothing straddles)
EIC=/sphenix/user/tmengel/eic
source $EIC/setup.sh
OUT=/sphenix/tg/tg01/jets/tmengel/eic/caen
hadd -ff -f $OUT/caen_3062.root $OUT/caen_3062_seg0.root $OUT/caen_3062_seg1.root $OUT/caen_3062_seg2.root || exit 1
ln -sfn $OUT/caen_3062.root $EIC/bic_testbeam/data/caen_3062.root
python3 - <<'PY' > $OUT/caen_3062_check.txt 2>&1
import uproot, numpy as np
f = uproot.open('/sphenix/tg/tg01/jets/tmengel/eic/caen/caen_3062.root')
print({k: f[k].num_entries for k in f.keys() if hasattr(f[k], 'num_entries')})
h = f['header'].arrays(['evt', 'seg', 'event_counter', 'gtt', 'ok_size', 'ok_marker'], library='np')
print('triggers', len(h['evt']), 'ok_size', h['ok_size'].mean(), 'ok_marker', h['ok_marker'].mean())
print('event_counter steps', dict(zip(*np.unique(np.diff(h['event_counter'].astype(np.int64)), return_counts=True))))
print('gtt backward steps', int((np.diff(h['gtt']) < 0).sum()), 'range s', h['gtt'][0] * 8.533e-9, h['gtt'][-1] * 8.533e-9)
print('RCDAQ evt gaps', int((np.diff(h['evt']) != 1).sum()), 'per segment', dict(zip(*np.unique(h['seg'], return_counts=True))))
for k in f.keys():
    if k.startswith('info'): print(k, str(f[k]))
PY
cat $OUT/caen_3062_check.txt
