"""Fit a resolved continuous current curve using integrated KCL."""
import numpy as np
from scipy.integrate import cumulative_trapezoid
from scipy.interpolate import BSpline
from calibrate_model import load
from check_models import show

def fit():
    data=load()
    # Cubic B-splines describe rounded knees without discontinuous bridges.
    knots=np.r_[[-8]*4,np.arange(-7.8,7.2,.2),[7.4]*4]
    xx=[];yy=[]
    for sw,fn,r,dt,v in data:
        if sw!='forward' or int(fn[5:-4])%2:continue
        n=round(30e-6/dt);ix=np.arange(0,len(v)-n,12)
        b=BSpline.design_matrix(v[:,0],knots,3,extrapolate=True).toarray()
        integ=cumulative_trapezoid(b,dx=dt,axis=0,initial=0)
        x=np.column_stack([(v[ix+n,0]-v[ix,0])*1e-9,(integ[ix+n]-integ[ix])*1e-3])
        ir=cumulative_trapezoid((v[:,1]-v[:,0])/(992+r),dx=dt,initial=0)
        y=ir[ix+n]-ir[ix]
        scale=np.sqrt(len(ix))*30e-6*1e-3
        # More weight on the core attractor than the saturation-only records.
        if r<315:scale*=3
        xx.append(x/scale);yy.append(y/scale)
    x=np.concatenate(xx);y=np.concatenate(yy)
    smooth=np.diff(np.eye(x.shape[1]-1),n=2,axis=0)*.015
    smooth=np.column_stack([np.zeros(len(smooth)),smooth])
    coeff=np.linalg.lstsq(np.vstack([x,smooth]),np.r_[y,np.zeros(len(smooth))],rcond=None)[0]
    print('c1',coeff[0],'resid',np.linalg.norm(x@coeff-y),flush=True)
    vk=np.linspace(-8,7.4,1541);ik=BSpline(knots,coeff[1:]*1e-3,3)(vk)
    np.savez('chua/curve_fit.npz',v=vk,i=ik,c1=coeff[0])
    return vk,ik,coeff[0]

def curve(vk,ik):
    # Linear extension beyond the observed range.
    gn=(ik[1]-ik[0])/(vk[1]-vk[0]);gp=(ik[-1]-ik[-2])/(vk[-1]-vk[-2])
    v=np.r_[-40,vk,40];i=np.r_[ik[0]+gn*(-40-vk[0]),ik,ik[-1]+gp*(40-vk[-1])]
    def g(x):return np.interp(x,v,i)
    g.knots=(v,i);g.segments=[(-40,0,gn,0),(0,40,gp,0)]
    return g

if __name__=='__main__':
    vk,ik,c1=fit();show('curve',curve(vk,ik),c1)
