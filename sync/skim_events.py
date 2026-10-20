"""Reduced event tree for the repository (< 100 MB): triggers recorded by KCU0, KCU1 and CAEN, or with an
AstroPix match (log10 LR > 1). Same tree name and branches as out/events_3062.root
-> out/events_3062_skim.root"""
import numpy as np, uproot
t = uproot.open('out/events_3062.root')['events']
d = t.arrays(library='np')
sel = ((d['has_kcu0'] == 1) & (d['has_kcu1'] == 1) & (d['has_caen'] == 1)) | (d['astro_log10LR'] > 1)
print(f"{sel.sum()} of {len(sel)} triggers kept")
with uproot.recreate('out/events_3062_skim.root', compression=uproot.LZMA(6)) as f:
    f.mktree('events', {k: (v.dtype, v.shape[1:]) for k, v in d.items()})
    f['events'].extend({k: v[sel] for k, v in d.items()})
