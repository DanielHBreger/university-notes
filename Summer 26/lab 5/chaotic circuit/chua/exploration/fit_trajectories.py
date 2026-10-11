"""Estimate one circuit model from short forward-record trajectories.

Back records and odd-numbered forward records are not fitted. Resistance is
fixed by the independently measured divider. Initial inductor currents are
nuisance parameters, since only the capacitor voltages were recorded.
"""
import json
import numpy as np
from scipy.optimize import least_squares
from calibrate_model import load
from scipy.signal import savgol_filter
from integration import trajectory

def diode(p):
    # p: C1 nF, C2 nF, L mH, rL ohm, Ga/Gbn/Gbp mS, Bn/Bp V, I0 mA.
    _,_,_,_,ga,gn,gp,bn,bp,i0=p
    vs=np.array([-40.,-6.75,-bn,bp,6.,40.])
    currents=np.empty(6);currents[2]=-ga*bn+i0;currents[3]=ga*bp+i0
    currents[1]=currents[2]+gn*(vs[1]-vs[2]);currents[0]=currents[1]+4*(vs[0]-vs[1])
    currents[4]=currents[3]+gp*(vs[4]-vs[3]);currents[5]=currents[4]+4*(vs[5]-vs[4])
    return vs,currents*1e-3

def windows():
    data=[r for r in load() if r[0]=='forward' and int(r[1][5:-4])%2==0]
    selected=[]
    for r in [380,450,550,650,680,720,740,755,770,870]:
        d=min(data,key=lambda a:abs(a[2]-r))
        for start in [2500,6500]:
            sw,fn,r,dt,v=d
            # 0.4 ms, downsampled for fit weights. Integration uses 0.5 us.
            ix=np.arange(start,start+round(.0004/dt),5)
            time=(ix-start)*dt
            dv=savgol_filter(v[:,1],21,4,deriv=1,delta=dt)
            selected.append((f'{sw}/{fn}',r,time,v[ix],dv[start]))
    return selected

def fit():
    win=windows()
    p0=np.array([11.5,100,18,3,-.77,-.413,-.413,1.02,.93,0])
    currents=[((v[0,0]-v[0,1])/(992+r)-100e-9*dv)*1e3 for name,r,t,v,dv in win]
    lo=[9,75,14,0,-.85,-.46,-.46,.7,.6,-.08]+[-25]*len(win)
    hi=[14,130,25,70,-.70,-.38,-.36,1.4,1.4,.08]+[25]*len(win)
    count=0
    def residual(params):
        nonlocal count
        p=params[:10];vk,ik=diode(p);res=[]
        for (name,r,t,v,dv),il in zip(win,params[10:]):
            y0=np.r_[v[0],il*1e-3]
            dt=.5e-6;steps=int(np.ceil(t[-1]/dt))
            pred=trajectory(y0,992+r,p[3],p[0]*1e-9,p[1]*1e-9,p[2]*1e-3,vk,ik,dt,steps,0)
            vals=np.column_stack([np.interp(t,np.arange(steps+1)*dt,pred[:,j]) for j in range(2)])
            res.extend(((vals-v)/[1,.5]).ravel())
        count+=1
        if count%100==0:print(count,'RMS',np.sqrt(np.mean(np.square(res))),p,flush=True)
        return res
    opt=least_squares(residual,np.r_[p0,currents],bounds=(lo,hi),x_scale='jac',diff_step=1e-4,max_nfev=100,ftol=1e-8)
    print('final',opt.x[:10], 'RMS',np.sqrt(np.mean(opt.fun**2)),opt.message,flush=True)
    np.savez('chua/trajectory_fit.npz',params=opt.x[:10],nuisance=opt.x[10:])
    print('records', [w[0] for w in win],flush=True)

if __name__=='__main__':fit()
