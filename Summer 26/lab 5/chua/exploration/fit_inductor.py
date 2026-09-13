"""Diagnose amplitude dependence of the measured inductor branch."""
import numpy as np
from scipy.signal import savgol_filter
from scipy.integrate import cumulative_trapezoid
from scipy.optimize import least_squares
from calibrate_model import load

def regress(data,c2=100):
    xx=[];yy=[];tags=[]
    for sw,fn,r,dt,v in data:
        w=max(7,round(35e-6/dt)|1)
        vs=savgol_filter(v,w,4,axis=0)
        dv=savgol_filter(v,w,4,deriv=1,delta=dt,axis=0)
        il=(vs[:,0]-vs[:,1])/(992+r)-c2*1e-9*dv[:,1]
        q=il*1e3
        n=round(60e-6/dt);idx=np.arange(w,len(v)-w-n,15)
        ints=cumulative_trapezoid(np.column_stack([il,il*q*q,np.ones(len(il))]),dx=dt,axis=0,initial=0)
        target=cumulative_trapezoid(vs[:,1],dx=dt,initial=0)
        x=np.column_stack([(il[idx+n]-il[idx])*1e-3,(il[idx+n]**3-il[idx]**3)*1e3/3,ints[idx+n]-ints[idx]])
        y=target[idx+n]-target[idx]
        scale=np.sqrt(len(idx))*60e-6
        xx.append(x/scale);yy.append(y/scale);tags.extend([f'{sw}/{fn}']*len(idx))
    x=np.concatenate(xx);y=np.concatenate(yy);tags=np.array(tags)
    train=np.array([z.startswith('forward/') and int(z.split('trace')[1][:-4])%2==0 for z in tags])
    c=np.linalg.lstsq(x[train],y[train],rcond=None)[0]
    print('C2',c2,'L0 mH, L2 mH/mA2, R0 ohm, R2 ohm/mA2, offset V',c,'relative residual',np.linalg.norm((x@c-y)[train])/np.linalg.norm(y[train]),flush=True)
    for label,mask in [('small',np.array([z.split('/')[1] not in [f'trace{k}.csv' for k in range(406,419)] for z in tags])&train),('back',~np.array([z.startswith('forward/') for z in tags]))]:
        print(label,np.linalg.norm((x@c-y)[mask])/np.linalg.norm(y[mask]),flush=True)
    return c

if __name__=='__main__':
    d=load()
    for c2 in [85,95,100,105,110]:regress(d,c2)
