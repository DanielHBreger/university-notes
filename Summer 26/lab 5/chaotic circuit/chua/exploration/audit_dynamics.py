"""Recover circuit parameters from simultaneous voltage records, in SI units."""
from pathlib import Path
import json
import numpy as np
import pandas as pd
from scipy.signal import savgol_filter
from scipy.optimize import least_squares

HERE = Path(__file__).resolve().parent

def samples(path, nrows=100000):
    d = pd.read_csv(path, nrows=nrows).to_numpy()
    dt = np.median(np.diff(d[:, 0]))
    w = max(7, int(round(12e-6 / dt)) | 1)
    v = savgol_filter(d[:, 1:3], w, 3, axis=0)
    dv = savgol_filter(d[:, 1:3], w, 3, deriv=1, delta=dt, axis=0)
    ddv = savgol_filter(d[:, 1:3], w, 3, deriv=2, delta=dt, axis=0)
    ix = np.arange(w, len(d)-w, max(1, int(round(4e-6/dt))))
    return v[ix], dv[ix], ddv[ix]

def main():
    rows=[]; arrays=[]
    for sweep in ['forward','back']:
        tab=pd.read_csv(HERE/f'{sweep}_rpot.csv')
        tab=tab[tab.status.eq('ok')]
        selected=tab.iloc[np.unique(np.linspace(0,len(tab)-1,32).astype(int))]
        for rec in selected.itertuples():
            v,dv,ddv=samples(HERE/sweep/rec.filename)
            r=992+rec.rpot_ohm
            # Fit C1 and one line on the negative shoulder, away from knees.
            mask=(v[:,0]>-4.8)&(v[:,0]<-1.6)
            X=np.column_stack([dv[:,0]*1e-9,v[:,0]*1e-3,np.ones(len(v))*1e-3])
            target=(v[:,1]-v[:,0])/r
            fit=np.linalg.lstsq(X[mask],target[mask],rcond=None)[0]
            rows.append(dict(sweep=sweep,filename=rec.filename,rpot=rec.rpot_ohm,
                             c1_nf=fit[0],gb_ms=fit[1],intercept_ma=fit[2],
                             rmse_ua=np.std((target-X@fit)[mask])*1e6))
            arrays.append((sweep,rec.filename,rec.rpot_ohm,v,dv,ddv))
    pd.DataFrame(rows).to_csv(HERE/'dynamics_audit.csv',index=False)
    np.savez_compressed(HERE/'dynamics_samples.npz',**{
        f'{sw}__{fn}__{rp}':np.column_stack([v,dv,ddv]) for sw,fn,rp,v,dv,ddv in arrays})
    print(pd.DataFrame(rows).round(3).to_string(index=False))
    # Node 2: v1'-v2' = R*C2*v2'' + (R/L)*v2 + (rL/L)*(R*C2*v2'-(v1-v2)).
    v=np.concatenate([a[3] for a in arrays]); dv=np.concatenate([a[4] for a in arrays]); ddv=np.concatenate([a[5] for a in arrays])
    R=np.concatenate([np.full(len(a[3]),992+a[2]) for a in arrays])
    def fun(p):
        c,l,rl=p[0]*1e-9,p[1]*1e-3,p[2]
        res=dv[:,0]-dv[:,1]-R*c*ddv[:,1]-R/l*v[:,1]-rl/l*(R*c*dv[:,1]-(v[:,0]-v[:,1]))
        # Independent channel DC offsets are nuisance parameters per record.
        start=0
        for a in arrays:
            end=start+len(a[3]); res[start:end]-=res[start:end].mean(); start=end
        return res/10000
    fit=least_squares(fun,[100,18,2],bounds=([50,10,0],[150,30,100]))
    print('C2 nF, L mH, rL ohm',fit.x,'residual',np.std(fit.fun))

if __name__=='__main__': main()
