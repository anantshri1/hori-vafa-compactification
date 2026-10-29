import numpy as np
from phase2_IB import chain
from zfun import ZB_rho, mp
p=(0.3,0.25,-2.0,0.8); r=p[0]
Zt=complex(ZB_rho(*p)); Ztp=complex(ZB_rho(p[0],p[1],p[2],p[3]+mp.pi))
for n in [959,961,1919]:
    eps,c=chain(*p,n,int(2**np.ceil(np.log2(160/0.02)))*(2 if n>1000 else 1)); v=sum(c.values())/eps**(2*r-1)
    print(f"n={n} (odd) eps={eps:.5f}: I_B/eps^(2r-1)={v:.8f}   Z(th)={Zt.real:.8f}  Z(th+pi)={Ztp.real:.8f}")
