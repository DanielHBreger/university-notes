"""Development comparison of circuit models against sweep transition ranges."""
import numpy as np
import simulate as s
from calibrate_model import load,fit_diode,integrated_data
from scipy.optimize import least_squares

def from_coeff(k,c):
    slopes=np.cumsum(c[2:])*1e-3
    bs=np.r_[-40,k,40]
    currents=(c[1]+c[2]*bs+sum(h*np.maximum(bs-b,0) for h,b in zip(c[3:],k)))*1e-3
    def g(v):return np.interp(v,bs,currents)
    g.knots=(bs,currents)
    g.segments=[(a,b,m,ia-m*a) for a,b,m,ia in zip(bs[:-1],bs[1:],slopes,currents[:-1])]
    return g

def show(tag,g,c1=11.5,c2=100,l=18,rl=0):
    s.set_components(c1*1e-9,c2*1e-9,l*1e-3)
    rs=np.arange(250,931,5.)
    f=s.sweep_continuation(g,rs,rl,.5e-6,.03,.025,'down')
    b=s.sweep_continuation(g,rs,rl,.5e-6,.03,.025,'up')
    np.savez(f'chua/check_{tag}.npz',r=rs,forward=np.array(f,dtype=object),back=np.array(b,dtype=object))
    print(tag,'hopf',s.hopf_point(g,rl),'back large until',max([r for r,m in zip(rs,b) if len(m) and min(m)>5],default=0),'forward large until',max([r for r,m in zip(rs,f) if len(m) and min(m)>5],default=0),flush=True)
    print([(r,s.describe(f[list(rs).index(r)])) for r in [330,400,600,675,700,725,750,775,800,900]],flush=True)

if __name__=='__main__':
    data=load()
    k,c=fit_diode([r for r in data if r[2]>330])
    # Independently fit saturation branches holding the inner three lines fixed.
    k0=k.copy();c0=c.copy()
    def fun(p,ret=False):
        kk=np.array([p[0],k0[1],k0[2],p[1]])
        cc=c0.copy();cc[1]=c0[1]-c0[3]*k0[0]+p[2]*kk[0]
        cc[2]=c0[2]+c0[3]-p[2];cc[3]=p[2];cc[-1]=p[3]
        x,y,_=integrated_data([r for r in data if r[0]=='forward' and r[2]<315],kk)
        return (kk,cc) if ret else x@cc-y
    opt=least_squares(fun,[-6.7,6.,-4.2,4.6],bounds=([-7.5,5.,-8,1],[-5.5,7.,-1,8]))
    k,c=fun(opt.x,True);np.savez('chua/check_fit.npz',knots=k,coeff=c)
    print('hybrid fit',k,c,flush=True)
    show('hybrid',from_coeff(k,c),c[0])
    show('old_sequential',s.make_g(symmetric=True,bp_scale=.8))
