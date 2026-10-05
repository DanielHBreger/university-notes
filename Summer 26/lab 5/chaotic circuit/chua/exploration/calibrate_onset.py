"""Test identifiability using the observed onset and first doubling."""
import numpy as np
from scipy.optimize import least_squares
from integration import hold
from scan_saturation import model
import simulate as s

def components(c2,rl,gb=-.413,rhopf=905,period=328):
    R=992+rhopf;w=2*np.pi/(period*1e-6);G=gb*1e-3
    def residual(p):
        c1,l=p[0]*1e-9,p[1]*1e-3
        A=np.array([[(-1/R-G)/c1,1/R/c1,0],[1/R/(c2*1e-9),-1/R/(c2*1e-9),-1/(c2*1e-9)],[0,1/l,-rl/l]])
        ev=np.linalg.eigvals(A);z=ev[np.argmax(ev.imag)]
        return [z.real/w,(z.imag-w)/w]
    fit=least_squares(residual,[11.5,18],bounds=([5,5],[30,40]),gtol=1e-12,xtol=1e-12,ftol=1e-12)
    return fit.x

def main():
    g=model(6.7,5.8,4)
    # measured lower equilibrium boundary constrains the central conductance
    for c2 in [80,100,120]:
        for rl in [0,10,30]:
            c1,l=components(c2,rl);s.set_components(c1*1e-9,c2*1e-9,l*1e-3)
            y=s.sweep_starts(g,[910],rl)[0][:,0]
            r1=0
            for r in np.arange(910,699,-5):
                y,m=hold(y,992+r,rl,s.C1,s.C2,s.L,*g.knots,1e-6,60000,20000)
                if len(m)>10 and np.mean(abs(m[-50:][1:]-m[-50:][:-1]))>.035:
                    r1=r;break
            print('C2 rL C1 L',c2,rl,c1,l,'R1',r1,flush=True)

if __name__=='__main__':main()
