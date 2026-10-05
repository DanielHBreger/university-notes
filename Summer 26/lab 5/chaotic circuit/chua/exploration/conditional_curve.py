"""Estimate g(v) from stationary conditional means; mean dv/dt at fixed v is zero."""
import numpy as np
from scipy.interpolate import UnivariateSpline
from scipy.integrate import cumulative_trapezoid
from calibrate_model import load
from fit_curve import curve
from check_models import show

def fit():
    data=[d for d in load() if d[0]=='forward' and int(d[1][5:-4])%2==0]
    vs=np.concatenate([d[4][:,0] for d in data])
    ir=np.concatenate([(d[4][:,1]-d[4][:,0])/(992+d[2]) for d in data])
    ix=np.floor((vs+8)/.1).astype(int)
    bins=np.arange(ix.max()+1);count=np.bincount(ix,minlength=len(bins))
    x=np.bincount(ix,weights=vs)/np.maximum(count,1)
    y=np.bincount(ix,weights=ir)/np.maximum(count,1)
    m=count>100
    # Saturation knees need finer curvature than shoulders; bounded weights
    # prevent turning-point dwelling from overwhelming the fit.
    spl=UnivariateSpline(x[m],y[m]*1e3,w=np.minimum(np.sqrt(count[m]/100),4),s=.04)
    vk=np.arange(-7.5,6.7,.01);ik=spl(vk)*1e-3
    g=curve(vk,ik)
    cs=[]
    for sw,fn,r,dt,v in data:
        n=round(40e-6/dt);idx=np.arange(0,len(v)-n,10)
        integ=cumulative_trapezoid((v[:,1]-v[:,0])/(992+r)-g(v[:,0]),dx=dt,initial=0)
        dx=v[idx+n,0]-v[idx,0];y=integ[idx+n]-integ[idx]
        c=np.dot(dx,y)/np.dot(dx,dx)*1e9
        cs.append((r,c))
    print('C1 per record',cs,flush=True)
    c1=np.median([c for r,c in cs if r>700]);print('c1',c1,flush=True)
    np.savez('chua/conditional_fit.npz',v=vk,i=ik,c1=c1)
    return g,c1

if __name__=='__main__':
    g,c1=fit();show('conditional',g,c1)
