"""Phase 1(b): background topological flux n for A.
A:  I_A(n) = sum_m oint  omega^{-m} z^{-m-n} Chir[z y1,m,r] Chir[z y2,m,r]
B (frame a=omega (y1y2)^{-1/2} on Q, nothing on Qt, R_Q=R_Qt=(1-r)/2):
    I_B(n) = sum_m oint  wB'^{-m} z^{m+n/2} Chir[z a, m+n, (1-r)/2] Chir[1/z, -m, (1-r)/2],
    wB' = (y1/y2) a^{-1/2}
Claim: I_A(n) = omega^{n/2} y1^{-n/4} y2^{3n/4} x^{n(3r-1)/4} I_B(n)."""
import numpy as np
from idx import chir
def IA_n(r,a1,a2,th,lx,n,M=150,N=1024):
    phi=2*np.pi*np.arange(N)/N; lz=1j*phi; tot=0
    for m in range(-M,M+1):
        tot+=(np.exp(-m*1j*th)*np.exp((-m-n)*lz)*chir(lz+1j*a1,m,r,lx)*chir(lz+1j*a2,m,r,lx)).mean()
    return tot
def IB_n(r,a1,a2,th,lx,n,M=150,N=1024,weyl=False):
    if weyl: a1,a2=a2,a1
    la=1j*th-0.5j*(a1+a2); lwB=1j*(a1-a2)-la/2; R=(1-r)/2
    phi=2*np.pi*np.arange(N)/N; lz=1j*phi; tot=0
    for m in range(-M,M+1):
        tot+=(np.exp(-m*lwB)*np.exp((m+n/2)*lz)*chir(lz+la,m+n,R,lx)*chir(-lz,-m,R,lx)).mean()
    return tot*np.exp(n*(0.5j*th-0.25j*a1+0.75j*a2)+n*(3*r-1)/4*lx)
if __name__=="__main__":
    pts=[(0.3,0.7,-1.9,0.8,np.log(0.35)),(0.41,2.3,0.4,-2.2,np.log(0.5)),(0.17,-0.6,1.1,2.9,np.log(0.25))]
    for p in pts:
        print('point',p[:4],'x=%.3f'%np.exp(p[4]))
        for n in [-3,-2,-1,0,1,2,3,5]:
            A=IA_n(*p,n); B=IB_n(*p,n); Bw=IB_n(*p,n,weyl=True)
            print(f"  n={n:+d} I_A={A:.14f}  |I_A-RHS|={abs(A-B):.1e}  (Weyl image y1<->y2: {abs(A-Bw):.1e})")
