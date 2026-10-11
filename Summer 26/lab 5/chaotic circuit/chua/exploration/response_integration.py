"""RK4 for an optional first-order response of the nonlinear element."""
import numpy as np
from numba import njit

@njit(cache=True)
def deriv(a,b,c,d,r,rl,c1,c2,l,tau,vk,ik):
    ir=(b-a)/r
    return (ir-d)/c1,(-ir-c)/c2,(b-rl*c)/l,(np.interp(a,vk,ik)-d)/tau

@njit(cache=True)
def step(a,b,c,d,h,r,rl,c1,c2,l,tau,vk,ik):
    a1,b1,c_1,d1=deriv(a,b,c,d,r,rl,c1,c2,l,tau,vk,ik)
    a2,b2,c_2,d2=deriv(a+h/2*a1,b+h/2*b1,c+h/2*c_1,d+h/2*d1,r,rl,c1,c2,l,tau,vk,ik)
    a3,b3,c_3,d3=deriv(a+h/2*a2,b+h/2*b2,c+h/2*c_2,d+h/2*d2,r,rl,c1,c2,l,tau,vk,ik)
    a4,b4,c_4,d4=deriv(a+h*a3,b+h*b3,c+h*c_3,d+h*d3,r,rl,c1,c2,l,tau,vk,ik)
    return a+h/6*(a1+2*a2+2*a3+a4),b+h/6*(b1+2*b2+2*b3+b4),c+h/6*(c_1+2*c_2+2*c_3+c_4),d+h/6*(d1+2*d2+2*d3+d4)

@njit(cache=True)
def hold(y,r,rl,c1,c2,l,tau,vk,ik,h,nsettle,ncollect):
    a,b,c,d=y
    for k in range(nsettle):a,b,c,d=step(a,b,c,d,h,r,rl,c1,c2,l,tau,vk,ik)
    mx=np.empty(ncollect//2+1);count=0;p2=a;p1=a
    for k in range(ncollect):
        a,b,c,d=step(a,b,c,d,h,r,rl,c1,c2,l,tau,vk,ik)
        if p1>p2 and p1>=a:mx[count]=p1;count+=1
        p2=p1;p1=a
    return np.array([a,b,c,d]),mx[:count]

@njit(cache=True)
def trajectory(y,r,rl,c1,c2,l,tau,vk,ik,h,nsteps,nskip):
    a,b,c,d=y;out=np.empty((nsteps-nskip+1,4))
    for k in range(nsteps+1):
        if k>=nskip:out[k-nskip]=a,b,c,d
        if k<nsteps:a,b,c,d=step(a,b,c,d,h,r,rl,c1,c2,l,tau,vk,ik)
    return out
