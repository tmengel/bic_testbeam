"""Unwrap AstroPix fpga_ts (32 bit, 40 MHz) -> out/astro_unwrapped.npz.
Keeps layer==1, chip id<=8 (the rest are corrupt packets)."""
import numpy as np, uproot
a=uproot.open('out/astro_3062.root')['astro'].arrays(library='np')
ok=(a['layer']==1)&(a['id']<=8)
ts=a['fpga_ts'][ok].astype(np.int64)
d=np.diff(ts); w=np.cumsum(np.r_[0,(d<-2**31).astype(np.int64)-(d>2**31).astype(np.int64)])
np.savez('out/astro_unwrapped.npz',tu=ts+w*2**32,rt=a['rt'][ok],ev=a['evt'][ok],id=a['id'][ok],loc=a['loc'][ok],iscol=a['iscol'][ok],tot=a['tot'][ok])
print("AstroPix hits:",ok.sum(),"wraps:",w.max())
