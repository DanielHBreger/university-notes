"""Sensitivity check: keep the inner law fixed and vary saturation only."""
import numpy as np
import simulate as s
from integration import hold
from check_models import from_coeff

def model(bn,bp,gs):
    k=np.array([-bn,-1.02,.93,bp]);sl=np.array([gs,-.413,-.770,-.413,gs])
    c=np.r_[11.45,.36-(sl[0]-sl[1])*k[0],sl[0],np.diff(sl)]
    return from_coeff(k,c)

def main():
    s.set_components(11.45e-9,100e-9,18e-3)
    for bn in [6.2,6.7,7.2]:
        for bp in [5.3,5.8,6.3]:
            for gs in [1.5,4.,8.]:
                g=model(bn,bp,gs);y=np.array([7.,6.,0.]);top=0;lowamp=None
                for r in np.arange(500,861,10):
                    y,m=hold(y,992+r,0.,s.C1,s.C2,s.L,*g.knots,1e-6,25000,10000)
                    if r==500:lowamp=(min(m),max(m)) if len(m) else(0,0)
                    if len(m) and min(m)>4.5:top=r
                    else:break
                print(bn,bp,gs,'crisis',top,'max at500',lowamp,flush=True)

if __name__=='__main__':main()
