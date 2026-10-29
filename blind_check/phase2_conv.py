import numpy as np, time, sys
from phase2_IB import chain
from zfun import ZB_rho, mp
pts=[(0.3,0.25,-2.0,0.8),(0.4,0.6,-1.0,2.0),(0.2,0.1,-3.0,1.3),(0.45,0.8,-0.5,2.7)]
which=[int(a) for a in sys.argv[1:]] or [0,1]
for i in which:
    p=pts[i]; Z=complex(ZB_rho(*p)); r=p[0]
    print('point',p,' Z_2d =',Z.real, flush=True)
    prev=None
    for n in [60,120,240,480,960,1920,3840]:
        t=time.time()
        eps,_=chain(*p,n,16)  # just to get eps
        N=int(2**np.ceil(np.log2(80/eps)))
        _,c=chain(*p,n,N); I=sum(c.values())
        _,c2=chain(*p,n,2*N); I2=sum(c2.values())
        ratio=I2/eps**(2*r-1)
        print(f"  n={n:5d} eps={eps:.5f} N={N} |I(N)-I(2N)|={abs(I-I2):.1e}  I_B/eps^(2r-1)={ratio.real:.10f}{ratio.imag:+.2e}i  err={abs(ratio-Z):.3e}  err/eps={abs(ratio-Z)/eps:.4f}  ({time.time()-t:.1f}s)",flush=True)
