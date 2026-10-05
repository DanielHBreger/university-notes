"""Validate a finite-response diode in true sequential sweeps."""
import numpy as np
from check_models import from_coeff
from response_integration import hold
import simulate as s
import sys

def main():
    suffix='_core' if '--core' in sys.argv else ''
    fit=np.load(f'chua/response_fit{suffix}.npz');p=fit['params'];g=from_coeff(p[2:],np.r_[p[0],fit['coeff']])
    if suffix:
        vs=np.array([-40.,-6.7,p[3],p[4],5.8,40.]);ii=g(vs)
        ii[0]=ii[1]+.004*(vs[0]-vs[1]);ii[-1]=ii[-2]+.004*(vs[-1]-vs[-2])
        g.knots=(vs,ii)
    rs=np.arange(250,931,5.);out=[]
    for direction in ['down','up']:
        y3=s.sweep_starts(g,[max(rs) if direction=='down' else min(rs)],0)[0 if direction=='down' else 1][:,0]
        y=np.r_[y3,g(y3[0])];mout=[None]*len(rs)
        order=np.arange(len(rs))[::-1] if direction=='down' else np.arange(len(rs))
        for ix in order:
            y,m=hold(y,992+rs[ix],0.,p[0]*1e-9,100e-9,18e-3,p[1]*1e-6,*g.knots,.25e-6,120000,100000)
            mout[ix]=m
        out.append(mout)
    np.savez(f'chua/check_response{suffix}.npz',r=rs,forward=np.array(out[0],dtype=object),back=np.array(out[1],dtype=object))
    for direction,mout in zip(['forward','back'],out):
        print(direction,'large until',max([r for r,m in zip(rs,mout) if len(m) and min(m)>5],default=0),flush=True)
        print([(r,s.describe(mout[list(rs).index(r)])) for r in [330,400,600,675,700,725,750,775,800,900]],flush=True)

if __name__=='__main__':main()
