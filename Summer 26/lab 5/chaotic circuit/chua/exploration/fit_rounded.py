"""Calibrate rounded diode knees using integrated capacitor KCL."""
import numpy as np
from scipy.integrate import cumulative_trapezoid
from scipy.optimize import least_squares
from calibrate_model import load
from fit_curve import curve
from check_models import show

def current(v,p):
    c1,ga,gn,gp,bn,bp,i0,w=p
    def soft(x):return .5*(x+np.sqrt(x*x+w*w))
    return (ga*v+i0+(gp-ga)*soft(v-bp)-(gn-ga)*soft(-v-bn))*1e-3

def make_g(p):
    v=np.linspace(-6.7,5.8,2501);i=current(v,p)
    v=np.r_[-40,v,40];i=np.r_[i[0]+.004*(-40+6.7),i,i[-1]+.004*(40-5.8)]
    def g(x):return np.interp(x,v,i)
    g.knots=(v,i);g.segments=[(-40,-6.7,.004,0),(-6.7,-1,p[2]*1e-3,0),(-1,1,p[1]*1e-3,0),(1,5.8,p[3]*1e-3,0),(5.8,40,.004,0)]
    return g

def fit():
    data=[d for d in load() if d[0]=='forward' and int(d[1][5:-4])%2==0 and 330<d[2]<865]
    prep=[]
    for sw,fn,r,dt,v in data:
        n=round(30e-6/dt);ix=np.arange(0,len(v)-n,12)
        integ=cumulative_trapezoid((v[:,1]-v[:,0])/(992+r),dx=dt,initial=0)
        prep.append((v[:,0],dt,ix,n,(integ[ix+n]-integ[ix]),(v[ix+n,0]-v[ix,0]),np.sqrt(len(ix))*30e-6*1e-3))
    def fun(p):
        res=[]
        for v,dt,ix,n,y,dx,scale in prep:
            gg=cumulative_trapezoid(current(v,p),dx=dt,initial=0)
            res.extend((p[0]*1e-9*dx+gg[ix+n]-gg[ix]-y)/scale)
        return res
    opt=least_squares(fun,[11.5,-.77,-.413,-.413,1.02,.93,0,.1],bounds=([9,-1,-.48,-.48,.5,.5,-.1,.002],[15,-.7,-.35,-.35,1.5,1.5,.1,1.5]),diff_step=1e-4,ftol=1e-9)
    print(opt.x,'cost',np.linalg.norm(opt.fun),flush=True)
    np.savez('chua/rounded_fit.npz',params=opt.x)
    return opt.x

if __name__=='__main__':
    p=fit();show('rounded',make_g(p),p[0])
