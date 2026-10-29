import numpy as np, sys
from idx import chir
r,mu,kap,th=0.3,0.25,-2.0,0.8
eps=float(sys.argv[1]); n=int(round(-kap/eps))
lx=-eps/2; R=(1-r)/2; la=1j*th; a1,a2=eps*mu,-eps*mu
lwB=1j*(a1-a2)-la/2
N=4096; phi=2*np.pi*np.arange(N)/N; lz=1j*phi
pref=np.exp(n*(0.5j*th-0.25j*a1+0.75j*a2)+n*(3*r-1)/4*lx)
cs={}
for m in range(-3*n,2*n+1):
    f=np.exp(-m*lwB)*np.exp((m+n/2)*lz)*chir(lz+la,m+n,R,lx)*chir(-lz,-m,R,lx)
    cs[m]=f.mean()*pref
tot=sum(cs.values())
print('eps',eps,'n',n,'I_B=',tot)
for m in range(-3*n,2*n+1, max(1,n//10)):
    print(m, f"{abs(cs[m]):.3e}", np.round(cs[m],6))
