"""Cross-check: brute-force I_A (original variables) at the same (eps,n) as the B chain evaluation."""
import numpy as np, time
from idx import chir
from phase2_IB import chain
p=(0.3,0.25,-2.0,0.8); r,mu,kap,th=p
for n in [120,240]:
    eps,c=chain(*p,n,4096*2); IB=sum(c.values())
    lx=-eps/2; N=int(2**np.ceil(np.log2(240/eps)))
    phi=2*np.pi*np.arange(N)/N; lz=1j*phi; IA=0; t=time.time()
    for m in range(-100,101):
        f=np.exp(-1j*th*m)*np.exp((-m-n)*lz)*chir(lz+1j*eps*mu,m,r,lx)*chir(lz-1j*eps*mu,m,r,lx)
        IA+=f.mean()
    print(f"n={n} eps={eps:.5f}: I_A={IA:.12f}  I_B(chain)={IB:.12f}  |diff|={abs(IA-IB):.1e}  ({time.time()-t:.0f}s)")
