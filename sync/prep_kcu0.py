"""KCU0 trigger list: unique timestamps of complete 20-sample waveforms
(stock decoder, ts_dump kcu output) -> out/kcu0_ts.npy (KCU ticks, sorted)."""
import numpy as np, uproot, glob
A=[uproot.open(f)['kcu'].arrays(['ts','nsamp'],library='np') for f in sorted(glob.glob('out/beam_ROC_0-00003062-*.root'))]
ts=np.concatenate([a['ts'] for a in A]).astype(np.int64); ns=np.concatenate([a['nsamp'] for a in A])
t=np.unique(ts[(ns==20)&(ts>1e8)])            # drop fragments and the stale ts~0 buffer flush
np.save('out/kcu0_ts.npy',t); print("KCU0 complete triggers:",len(t))
