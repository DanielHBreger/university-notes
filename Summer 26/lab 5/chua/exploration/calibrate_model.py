"""Estimate a continuous diode model and passive dynamics from retained sweeps.

Fits local circuit equations, not the bifurcation diagram. Back sweep records
are reserved for validation. The capacitor equation is integrated over 40 us
windows to avoid differentiating quantisation noise.
"""
from pathlib import Path
import json
import numpy as np
import pandas as pd
from scipy.signal import resample_poly
from scipy.integrate import cumulative_trapezoid
from scipy.optimize import least_squares

HERE=Path(__file__).resolve().parent

def read_sample(path):
    d=pd.read_csv(path,nrows=160000).to_numpy()
    dt=float(np.median(np.diff(d[:,0])))
    dec=max(1,round(2e-6/dt))
    v=resample_poly(d[:,1:3],1,dec)[50:-50]
    return dt*dec,v

def prepare():
    data={}
    for sw in ['forward','back']:
        tab=pd.read_csv(HERE/f'{sw}_rpot.csv'); tab=tab[tab.status.eq('ok')]
        # Sample the entire dial, including the outer cycle at low resistance.
        inds=np.unique(np.r_[np.linspace(0,len(tab)-1,40).astype(int),np.arange(len(tab)-15,len(tab))])
        for row in tab.iloc[inds].itertuples():
            dt,v=read_sample(HERE/sw/row.filename)
            data[f'{sw}__{row.filename}__{row.rpot_ohm}__{dt}']=v
    np.savez_compressed(HERE/'calibration_samples.npz',**data)

def load():
    data=[]
    for key,v in np.load(HERE/'calibration_samples.npz').items():
        sw,fn,r,dt=key.split('__');data.append((sw,fn,float(r),float(dt),v))
    return data

def basis(v,knots):
    return np.column_stack([np.ones(len(v)),v]+[np.maximum(v-k,0) for k in knots])

def integrated_data(data,knots):
    xx=[];yy=[];groups=[]
    for sw,fn,r,dt,v in data:
        n=max(1,round(40e-6/dt))
        ix=np.arange(0,len(v)-n,10)
        integ=cumulative_trapezoid(basis(v[:,0],knots),dx=dt,axis=0,initial=0)
        x=np.column_stack([(v[ix+n,0]-v[ix,0])*1e-9,(integ[ix+n]-integ[ix])*1e-3])
        ir=cumulative_trapezoid((v[:,1]-v[:,0])/(992+r),dx=dt,initial=0)
        y=ir[ix+n]-ir[ix]
        # Equal weighting of each trace, regardless of sampling rate/length.
        scale=np.sqrt(len(ix))*40e-6*1e-3
        xx.append(x/scale);yy.append(y/scale);groups.extend([f'{sw}/{fn}']*len(ix))
    return np.concatenate(xx),np.concatenate(yy),np.array(groups)

def fit_diode(data):
    # Inner and saturation knees are optimized; every segment joins exactly.
    train=[r for r in data if r[0]=='forward' and int(r[1][5:-4])%2==0]
    def fun(k,ret=False):
        x,y,groups=integrated_data(train,k)
        coeff=np.linalg.lstsq(x,y,rcond=None)[0]
        return (coeff,x@coeff-y) if ret else x@coeff-y
    opt=least_squares(fun,[-6.7,-1.,1.,5.8],bounds=([-7.5,-1.6,.4,5.],[ -5.8,-.4,1.6,6.7]),diff_step=1e-3)
    coeff,res=fun(opt.x,True)
    print('knots',opt.x,'C1, intercept, slopes/hinges',coeff,'cost',np.linalg.norm(res))
    print('slopes mS',np.cumsum(coeff[2:]))
    return opt.x,coeff

def fit_passive(data):
    freq=[];f1=[];f2=[];rr=[]
    for sw,fn,r,dt,v in data:
        if sw!='forward' or int(fn[5:-4])%2: continue
        # Hann window suppresses record-edge leakage; retain spectral peaks.
        z=np.fft.rfft((v-v.mean(axis=0))*np.hanning(len(v))[:,None],axis=0)
        f=np.fft.rfftfreq(len(v),dt)
        mask=(f>1500)&(f<15000)&(abs(z[:,1])>0.04*max(abs(z[:,1])))
        freq.extend(f[mask]);f1.extend(z[mask,0]);f2.extend(z[mask,1]);rr.extend([992+r]*sum(mask))
    w=2j*np.pi*np.array(freq);a=np.array(f1);b=np.array(f2);rr=np.array(rr)
    def fun(p):
        c,l,rl,gain=p[0]*1e-9,p[1]*1e-3,p[2],p[3]
        h=gain/(1+rr*(w*c+1/(rl+w*l)))
        res=(a*h-b)/np.linalg.norm(b)
        return np.r_[res.real,res.imag]
    fit=least_squares(fun,[100,18,2,1],bounds=([70,12,0,.95],[130,25,50,1.05]),xtol=1e-11,gtol=1e-11,ftol=1e-11)
    print('C2 nF, L mH, rL ohm, channel gain',fit.x,'relative residual',np.linalg.norm(fit.fun))
    return fit.x

def main():
    if not (HERE/'calibration_samples.npz').exists(): prepare()
    data=load()
    k,c=fit_diode(data); passive=fit_passive(data)
    np.savez(HERE/'calibration_draft.npz',knots=k,coeff=c,passive=passive)
    x,y,groups=integrated_data(data,k)
    rows=[]
    for name in np.unique(groups):
        m=groups==name
        rows.append(dict(record=name,relative_integral_error=np.linalg.norm((x@c-y)[m])/np.linalg.norm(y[m])))
    pd.DataFrame(rows).to_csv(HERE/'calibration_validation.csv',index=False)
    print(pd.DataFrame(rows).to_string(index=False))

if __name__=='__main__':main()
