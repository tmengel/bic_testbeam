#!/bin/bash
# Wait for a run_to_root cluster to drain, check every segment log, and submit run_merge.job (plain condor_submit) for the
# stream/run pairs whose segments all exited 0. Failed or missing segments are listed and their pairs are not merged.
#   run_watch.sh <cluster> [segment list = run_to_root_list.txt]
CL=$1; LIST=${2:-run_to_root_list.txt}; LOGS=/sphenix/user/tmengel/eic/output/root_logs
cd /sphenix/user/tmengel/eic/bic_testbeam/decoding/condor
while condor_q $CL -af ClusterId 2>/dev/null | grep -q .; do
  h=$(condor_q $CL -af JobStatus 2>/dev/null | grep -c '^5$')
  [ "$h" -gt 0 ] && echo "$(date +%H:%M) $h job(s) held in cluster $CL"
  sleep 300
done
echo "$(date +%H:%M) cluster $CL drained"
python3 - "$LIST" "$LOGS" <<'PY'
import sys, collections, re
lst, logs = sys.argv[1], sys.argv[2]
ok, bad = collections.defaultdict(list), collections.defaultdict(list)
for line in open(lst):
    st, run, seg = line.split()
    try: txt = open('%s/%s_%s_seg%s.out' % (logs, st, run, seg)).read()
    except OSError: txt = ''
    (ok if re.search(r'^exit 0,', txt, re.M) else bad)[(st, run)].append(seg)
good = sorted(k for k in ok if k not in bad)
open('run_merge_ready.txt', 'w').write(''.join('%s %s\n' % k for k in good))
print('segments ok: %d, failed: %d' % (sum(map(len, ok.values())), sum(map(len, bad.values()))))
for k, v in sorted(bad.items()): print('FAILED %s %s segments %s' % (k[0], k[1], ' '.join(v)))
print('merging %d stream/run pairs' % len(good))
PY
[ -s run_merge_ready.txt ] && condor_submit -a "list = run_merge_ready.txt" run_merge.job | tail -1
