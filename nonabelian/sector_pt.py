"""LOG 35: pointwise test of the Bessel-regime LG on exact 3d mirror sectors (qpoch parts; zero points/monomials exact).
Fluxes m = (L,T,B,R) = (-2,-1,-3,-2) nu* + delta, eta-flux n = 4 nu*, x^{nu*} ~ eps. Prediction: product over the six light
chirals of exp[(zeta x^{R+|nu|} - x^{2-R+|nu|}/zeta)/(1-q)] (with zeta -> zeta^{-1} for nu<0), vertical chiral -> 1 + O(eps)."""
import numpy as np
from idx import qpoch
r = 0.3
CH = [((0,0,1,0), 1, r), ((1,0,-1,0), 0, r), ((-1,1,0,0), 0, r), ((0,-1,0,0), 0, r),
      ((0,0,-1,1), 0, r), ((0,1,0,-1), 0, r), ((0,-1,1,0), 0, 2-2*r)]
rng = np.random.default_rng(1)
eta = np.exp(-0.9j)
for eps in (0.2, 0.1, 0.05):
    x = np.exp(-eps/2); q = x*x; ns = int(round(2*np.log(1/eps)/eps)); n = 4*ns
    errs, verts = [], []
    for trial in range(8):
        d = rng.integers(-3, 4, 4); m = np.array([-2, -1, -3, -2])*ns + d
        phi = rng.uniform(-np.pi, np.pi, 4)
        ex = 1; pr = 1
        for Q, haseta, R in CH:
            nu = int(np.dot(Q, m) + (n if haseta else 0))
            zeta = (-1)**nu*np.exp(1j*np.dot(Q, phi))*(eta**-1 if haseta else 1)
            a = abs(nu)
            val = qpoch(x**(2 - R + a)/zeta, q)/qpoch(zeta*x**(R + a), q)
            ex *= val
            if Q == (0,-1,1,0): verts.append(abs(val - 1)); continue
            pr *= np.exp((zeta*x**(R + a) - x**(2 - R + a)/zeta)/(1 - q))
        errs.append(abs(ex/pr - 1))
    print(f"eps={eps} nu*={ns}: max |exact/LG - 1| = {max(errs):.2e}; vertical chiral |val-1| max = {max(verts):.2e}", flush=True)
