"""Coarse AstroPix->KCU0 offset from beam on/off structure: CAEN (no readout loss,
mapped onto the KCU0 clock) defines beam-on; FFT cross-correlation with the
AstroPix cluster rate gives the lag."""
import numpy as np, uproot
tk=np.load('out/kcu0_ts.npy').astype(np.int64)
m=np.load('out/caen_kcu0_match.npz')
g=uproot.open('out/caen_3062.root')['caen'].arrays(['gtt'],library='np')['gtt'].astype(float)
p=np.polyfit(g[m['icaen']],tk[m['ikcu0']].astype(float),1); tc=np.polyval(p,g)
cl=np.load('out/astro_clusters.npz'); x=(cl['Fmin']+int(cl['coarse'])).astype(float)
def xcorr(a,b):                                   # r[lag]: b(t+lag) vs a(t)
    a=a-a.mean(); b=b-b.mean(); n=len(a)
    f=np.fft.irfft(np.conj(np.fft.rfft(a,2*n))*np.fft.rfft(b,2*n),2*n)
    return np.r_[f[-n+1:],f[:n]]/np.sqrt((a*a).sum()*(b*b).sum()), np.arange(-n+1,n)
for bwms,maxlag_s in ((10,30),(1,3)):
    bw=bwms*1e6/6.25; e=np.arange(tk[0]-40e9/6.25,tk[-1]+40e9/6.25,bw)
    rc=np.histogram(tc,e)[0].astype(float); ra=np.histogram(x,e)[0].astype(float)
    k=np.ones(max(1,int(100/bwms)))/max(1,int(100/bwms)); sm=np.convolve(rc,k,'same')
    onS=(sm>0.1*np.median(sm[sm>0])).astype(float)
    r,l=xcorr(onS,ra); lag_s=l*bwms/1e3
    s=np.abs(lag_s)<=maxlag_s; rs,ls=r[s],lag_s[s]; o=np.argsort(rs)[::-1][:5]
    print(f"bin {bwms:2d} ms: AstroPix time sits {np.round(ls[o],3)} s after KCU0 (r {np.round(rs[o],3)}), median r {np.median(rs):.3f}")
