"""LOG 14: Coulomb limit of the planar mirror, zero gauge flux, eta-flux -n (n ~ 1/eps, the image of the electric Higgs-limit
background). Measure |integrand| along rays sigma = s v + sigma0 (z_a = q^{-i sigma_a}) to find non-compact (continuum) directions.
Chirals (charge on L,T,B,R; R): B(0,0,1,0;r, eta-flux), L-B, T-L, -T, R-B, T-R (r), B-T (2-2r)."""
import numpy as np
from idx import qpoch
r = 0.3
CH = [((0,0,1,0), r, 1), ((1,0,-1,0), r, 0), ((-1,1,0,0), r, 0), ((0,-1,0,0), r, 0),
      ((0,0,-1,1), r, 0), ((0,1,0,-1), r, 0), ((0,-1,1,0), 2-2*r, 0)]
def logmod(sig, eps, n):
    x = np.exp(-eps/2); q = x*x; tot = 0.0
    for Q, R, eta in CH:
        nu = -n*eta  # eta-flux -n
        ze = q**(-1j*np.dot(Q, sig))
        v = (ze*x**(R-1))**(-abs(nu)/2)*qpoch(x**(2-R+abs(nu))/ze, q)/qpoch(ze*x**(R+abs(nu)), q)
        tot += np.log(abs(v))
    return tot
eps = 0.02; n = 50
s0 = np.array([0.13, -0.21, 0.37, 0.05])
dirs = {'(1,1,1,1)': (1,1,1,1), '(1,0,0,0)': (1,0,0,0), '(0,1,0,0)': (0,1,0,0), '(0,0,1,0)': (0,0,1,0),
        '(1,1,1,0)': (1,1,1,0), '(1,0,1,1)': (1,0,1,1), '(1,2,3,4) generic': (1,2,3,4)}
ss = np.array([4., 8., 16., 32.])
for name, v in dirs.items():
    L = [logmod(s*np.array(v, float) + s0, eps, n) for s in ss]
    slopes = np.diff(L)/np.diff(np.log(ss))
    print(f"{name:20s} power-law exponents: " + " ".join(f"{p:+.3f}" for p in slopes), flush=True)
