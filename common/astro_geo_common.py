"""Shared inputs for the AstroPix geometry plots.

True-partner distributions are obtained WITHOUT any physics selection: every
phase-matched KCU0 candidate (10 MHz phase 189+-1) in the population-consistent
latency window gets weight +1, every control-band candidate (phase 60+-1, 120+-1)
gets -1/2. Weighted sums are then background-subtracted true-partner sums.
Candidates are common-trigger table rows (out/triggers_3062.npz): ki is a row there,
i_kcu0/i_kcu1/i_caen give the stream rows (-1 when that stream missed the trigger)."""
import numpy as np, uproot
def cluster_xy():
    """per AstroPix cluster: chip, first row and column hit (-1 if none); x, y in mm (nominal pitch)"""
    raw=uproot.open('out/astro_3062.root')['astro'].arrays(['layer','id','loc','iscol'],library='np')
    ok=(raw['layer']==1)&(raw['id']<=8)
    loc=raw['loc'][ok].astype(int); isc=raw['iscol'][ok].astype(int)
    fi=np.load('out/astro_cluster_firsthit.npy'); nh=len(loc); ncl=len(fi)
    cid=np.searchsorted(fi,np.arange(nh),side='right')-1
    A=np.load('out/astro_attr.npz'); chip=A['chip'].astype(int)
    row=np.full(ncl,-1); col=np.full(ncl,-1)
    for flag,arr in ((0,row),(1,col)):
        h=np.nonzero((isc==flag)&(loc<=34))[0]
        u,first=np.unique(cid[h],return_index=True); arr[u]=loc[h[first]]
    x=np.where(col>=0,chip*20.0+col*0.5+0.25,np.nan); y=np.where(row>=0,row*0.5+0.25,np.nan)
    return chip,row,col,x,y
def load():
    chip,row,col,_,_=cluster_xy(); fi=np.load('out/astro_cluster_firsthit.npy'); A=np.load('out/astro_attr.npz')
    R=np.load('out/astro_matcher_rows.npz'); cl,ki,band,lat=R['cl'],R['ki'],R['band'],R['lat']
    toend=A['to_end'][cl]
    keep=((lat>-26)&(lat<-7)&(toend<20))|((lat>-97)&(lat<-76)&(toend>13))
    cl,ki,band=cl[keep],ki[keep],band[keep]
    w=np.where(band==0,1.0,-0.5)
    TT=np.load('out/triggers_3062.npz')
    return dict(chip=chip,row=row,col=col,cl=cl,ki=ki,w=w,n=len(fi),i_kcu0=TT['i_kcu0'][ki],i_kcu1=TT['i_kcu1'][ki],
                i_caen=TT['i_caen'][ki],pad=TT['pad'][ki])
def wmean(x, w, o, edges):
    """background-subtracted mean of o in bins of x, with error"""
    b=np.digitize(x,edges)-1; nb=len(edges)-1; ok=(b>=0)&(b<nb)
    sw=np.bincount(b[ok],w[ok],nb); swo=np.bincount(b[ok],w[ok]*o[ok],nb)
    m=np.where(sw>0,swo/np.where(sw>0,sw,1),np.nan)
    var=np.bincount(b[ok],(w[ok]**2)*(o[ok]-m[b[ok]])**2,nb)/np.where(sw>0,sw,1)**2
    return m,np.sqrt(var),sw
