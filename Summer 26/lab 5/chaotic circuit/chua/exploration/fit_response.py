"""Test a finite-response nonlinear element against integrated KCL.

tau*dj/dt+j=g(v1), j=(v2-v1)/R-C1*dv1/dt.
Thus integral(ir)+tau*delta(ir) = C1*delta(v1) + tau*C1*delta(dv1) + integral(g).
"""
import numpy as np
import sys
from scipy.integrate import cumulative_trapezoid
from scipy.signal import savgol_filter
from scipy.optimize import least_squares
from calibrate_model import load,basis

def fit():
    data=[d for d in load() if d[0]=='forward' and int(d[1][5:-4])%2==0]
    if '--core' in sys.argv:
        data=[d for d in data if 330<d[2]<865]
    prep=[]
    for sw,fn,r,dt,v in data:
        w=max(9,round(25e-6/dt)|1);v=savgol_filter(v,w,4,axis=0)
        dv=savgol_filter(v[:,0],w,4,deriv=1,delta=dt)
        n=round(40e-6/dt);ix=np.arange(w,len(v)-n-w,12)
        ir=(v[:,1]-v[:,0])/(992+r)
        integ=cumulative_trapezoid(ir,dx=dt,initial=0)
        scale=np.sqrt(len(ix))*40e-6*1e-3
        prep.append((v[:,0],dt,ix,n,integ[ix+n]-integ[ix],v[ix+n,0]-v[ix,0],dv[ix+n]-dv[ix],ir[ix+n]-ir[ix],scale))
    def fun(p,ret=False):
        c,tau=p[:2];k=p[2:];xx=[];yy=[]
        for v,dt,ix,n,integ,dx,ddx,dir,scale in prep:
            b=cumulative_trapezoid(basis(v,k),dx=dt,axis=0,initial=0)
            xx.append((b[ix+n]-b[ix])*1e-3/scale)
            yy.append((integ+tau*1e-6*dir-c*1e-9*dx-tau*c*1e-15*ddx)/scale)
        x=np.concatenate(xx);y=np.concatenate(yy);coef=np.linalg.lstsq(x,y,rcond=None)[0]
        return (coef,x@coef-y) if ret else x@coef-y
    opt=least_squares(fun,[10.8,1.5,-6.7,-1.,.95,6.],bounds=([6,0,-7.5,-1.5,.5,5],[15,15,-5.5,-.5,1.5,7]),diff_step=1e-3,max_nfev=50)
    c,res=fun(opt.x,True);print('C1 tau knots',opt.x,'coeff',c,'cost',np.linalg.norm(res),flush=True)
    suffix='_core' if '--core' in sys.argv else ''
    np.savez(f'chua/response_fit{suffix}.npz',params=opt.x,coeff=c)

if __name__=='__main__':fit()
