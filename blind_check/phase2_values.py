"""Phase 2(c): limiting 2d expression at the four points by two independent methods:
 (i) Bessel integral from the B-side limit, (ii) CP^1 Coulomb-branch residue sum (A-side GLSM)."""
from zfun import *
mp.mp.dps = 25
pts=[(0.3,0.25,-2.0,0.8),(0.4,0.6,-1.0,2.0),(0.2,0.1,-3.0,1.3),(0.45,0.8,-0.5,2.7)]
for p in pts:
    b = ZB_rho(*p); a = Z_CP1(*p, M=int(90/abs(p[2]))+20, K=int(90/abs(p[2]))+20)
    r,mu,kap,th = p
    print(p, ' Z_B(Bessel)=', mp.nstr(b.real,16), ' Z_CP1(residues)=', mp.nstr(a.real,16),
          ' |diff|=', mp.nstr(abs(a-b),3), ' bare integral (no e^{kap r/2}, no 4):', mp.nstr((b/(4*mp.exp(kap*r/2))).real,16), flush=True)
