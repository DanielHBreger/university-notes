"""Allow weak shoulder curvature that a five-line fit cannot resolve."""
import numpy as np
from scipy.integrate import cumulative_trapezoid
from scipy.optimize import least_squares
from calibrate_model import load
from check_models import show

def basis(v,k):
    left=np.minimum(v-k[0],0);right=np.maximum(v-k[1],0)
    return np.column_stack([np.ones(len(v)),v,left,right,left**2,right**2,left**3,right**3])

def make_g(k,c):
    v=np.linspace(-6.7,5.8,2501);i=basis(v,k)@c*1e-3
    v=np.r_[-40,v,40];i=np.r_[i[0]+.004*(-40+6.7),i,i[-1]+.004*(40-5.8)]
    def g(x):return np.interp(x,v,i)
    g.knots=(v,i);g.segments=[(-40,-6.7,.004,0),(-6.7,5.8,-.000413,0),(5.8,40,.004,0)]
    return g

def fit():
    data=[d for d in load() if d[0]=='forward' and int(d[1][5:-4])%2==0 and 330<d[2]<900]
    prep=[]
    for sw,fn,r,dt,v in data:
        n=round(30e-6/dt);ix=np.arange(0,len(v)-n,12)
        integ=cumulative_trapezoid((v[:,1]-v[:,0])/(992+r),dx=dt,initial=0)
        prep.append((v[:,0],dt,ix,n,integ[ix+n]-integ[ix],v[ix+n,0]-v[ix,0],np.sqrt(len(ix))*30e-6*1e-3))
    def fun(k,ret=False):
        xx=[];yy=[]
        for v,dt,ix,n,y,dx,scale in prep:
            b=cumulative_trapezoid(basis(v,k),dx=dt,axis=0,initial=0)
            xx.append(np.column_stack([dx*1e-9,(b[ix+n]-b[ix])*1e-3])/scale)
            yy.append(y/scale)
        x=np.concatenate(xx);y=np.concatenate(yy)
        regular=np.diag([0,0,0,0,0,1,1,3,3])
        c=np.linalg.lstsq(np.vstack([x,regular]),np.r_[y,np.zeros(9)],rcond=None)[0]
        return (c,x@c-y) if ret else np.r_[x@c-y,regular@c]
    opt=least_squares(fun,[-1.02,.90],bounds=([-1.5,.5],[-.5,1.5]),diff_step=1e-4)
    c,res=fun(opt.x,True);print('knots',opt.x,'C1,coeff',c,'cost',np.linalg.norm(res),flush=True)
    np.savez('chua/curvature_fit.npz',knots=opt.x,coeff=c)
    return opt.x,c

if __name__=='__main__':
    k,c=fit();show('curvature',make_g(k,c[1:]),c[0])
