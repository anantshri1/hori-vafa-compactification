"""LOG 28: planar-mirror monomials vs the Gr(2,4) ladder (EHX/BCFKvS) W = X + Y/X + Z/X + W/Y + W/Z + q/W.
Our 7 chirals (charges on L,T,B,R) give monomials w^{Q}: B, L/B, T/L, 1/T, R/B, T/R, B/T. Count critical points of
W_c = ladder + c X/W for c = 0 and generic c (homotopy-free: random-start Newton in log variables, dedupe)."""
import numpy as np
from scipy.optimize import fsolve
q = 1.3*np.exp(0.4j)
def grad(lv, c):
    X, Y, Z, W = np.exp(lv)
    # x d/dx of W_c
    return np.array([X - (Y+Z)/X + c*X/W, Y/X - W/Y, Z/X - W/Z, W/Y + W/Z - q/W - c*X/W])
def roots(c, tries=4000, rng=np.random.default_rng(0)):
    sols = []
    for _ in range(tries):
        z0 = rng.normal(0, 1.5, 4) + 1j*rng.uniform(-np.pi, np.pi, 4)
        f = lambda v: np.concatenate([grad(v[:4] + 1j*v[4:], c).real, grad(v[:4] + 1j*v[4:], c).imag])
        v, info, ier, _ = fsolve(f, np.concatenate([z0.real, z0.imag]), full_output=True)
        if ier == 1 and np.linalg.norm(f(v)) < 1e-10:
            p = np.exp(v[:4] + 1j*v[4:])
            if all(np.linalg.norm(p - s) > 1e-6*np.linalg.norm(p) for s in sols): sols.append(p)
    return sols
for c in (0.0, 0.7*np.exp(0.3j), 1.0, -1.0):
    s = roots(c); print(f"c={c}: {len(s)} critical points; W values:", np.round(sorted([p[0]+(p[1]+p[2])/p[0]+p[3]/p[1]+p[3]/p[2]+q/p[3]+c*p[0]/p[3] for p in s], key=lambda z: z.real), 4), flush=True)
