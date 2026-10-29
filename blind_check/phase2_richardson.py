import numpy as np, sys
from phase2_IB import chain
from zfun import ZB_rho
pts=[(0.3,0.25,-2.0,0.8),(0.4,0.6,-1.0,2.0),(0.2,0.1,-3.0,1.3),(0.45,0.8,-0.5,2.7)]
for p in pts:
    Z=complex(ZB_rho(*p)); r=p[0]; E=[];V=[]
    for n in [240,480,960,1920,3840]:
        eps,_=chain(*p,n,16); N=int(2**np.ceil(np.log2(160/eps)))
        _,c=chain(*p,n,N); E.append(eps); V.append(sum(c.values())/eps**(2*r-1))
    E=np.array(E); V=np.array(V)
    for deg in [1,2,3]:
        A=np.vander(E[-(deg+1):],deg+1,increasing=True)
        coef=np.linalg.solve(A,V[-(deg+1):])
        print(p,f"poly-in-eps extrapolation deg {deg}: Z_extrap={coef[0]:.12f}  |Z_extrap - Z_2d|={abs(coef[0]-Z):.2e}  O(eps) coeff={coef[1]:.5f}",flush=True)
