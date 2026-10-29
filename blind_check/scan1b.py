import numpy as np
from idx import chir
def IA_n(r,a1,a2,th,lx,n,M=80,N=1024):
    phi=2*np.pi*np.arange(N)/N; lz=1j*phi; tot=0
    for m in range(-M,M+1):
        f=np.exp(-m*1j*th)*np.exp((-m-n)*lz)*chir(lz+1j*a1,m,r,lx)*chir(lz+1j*a2,m,r,lx)
        tot+=f.mean()
    return tot
def IB_n(r,a1,a2,th,lx,n,c,M=80,N=1024):
    lab=1j*th-0.5j*(a1+a2); lwB=1j*(a1-a2)-lab/2; R=(1-r)/2
    phi=2*np.pi*np.arange(N)/N; lz=1j*phi; tot=0
    for m in range(-M,M+1):
        f=np.exp(-m*lwB)*np.exp((m+c*n)*lz)*chir(lz+lab,m+n,R,lx)*chir(-lz,-m,R,lx)
        tot+=f.mean()
    return tot
pts=[(0.3,0.7,-1.9,0.8,np.log(0.35)),(0.3,0.7+0.01,-1.9,0.8,np.log(0.35)),(0.3,0.7,-1.9+0.01,0.8,np.log(0.35)),
     (0.3,0.7,-1.9,0.81,np.log(0.35)),(0.3,0.7,-1.9,0.8,np.log(0.36)),(0.31,0.7,-1.9,0.8,np.log(0.35))]
for n in [1,2,-1]:
  for c in [-1.5,-0.5,0.5,1.5]:
    R=[IA_n(*p,n)/IB_n(*p,n,c) for p in pts]
    L=np.log(R)
    d=[(L[i]-L[0])/0.01 for i in range(1,4)]
    dx=(L[4]-L[0])/(np.log(0.36)-np.log(0.35)); dr=(L[5]-L[0])/0.01
    print(n,c,'ratio',np.round(R[0],6),'dlog/d(a1,a2,th)/i:',np.round(np.array(d)/1j,4),'dlog/dlx',np.round(dx,4),'dlog/dr',np.round(dr,4))
print('----C_n')
for n in [-3,-2,-1,1,2,3,4]:
  for p in pts[:1]+[(0.41,2.3,0.4,-2.2,np.log(0.5)),(0.17,-0.6,1.1,2.9,np.log(0.25))]:
    r,a1,a2,th,lx=p
    mono=np.exp(n*(0.5j*th-0.25j*a1+0.75j*a2)+n*(3*r-1)/4*lx)
    print(n, np.round(IA_n(*p,n,M=150)/IB_n(*p,n,0.5,M=150)/mono,12))
