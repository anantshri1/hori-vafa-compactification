"""LOG 31: ladder with twisted masses entering as linear (log) terms: w d/dw [W_ladder + mu.log w] = 0."""
import numpy as np
from scipy.optimize import fsolve
rng = np.random.default_rng(11)
lad = np.array([(1,0,0,0),(-1,1,0,0),(-1,0,1,0),(0,-1,0,1),(0,0,-1,1),(0,0,0,-1)], float)
def count(E, c, mu, tries=6000):
    sols = []
    def g(v):
        l = v[:4] + 1j*v[4:]; r = E.T @ (c*np.exp(E @ l)) + mu
        return np.concatenate([r.real, r.imag])
    for _ in range(tries):
        v0 = np.concatenate([rng.normal(0, 2, 4), rng.uniform(-np.pi, np.pi, 4)])
        with np.errstate(all='ignore'):
            v, info, ier, _ = fsolve(g, v0, full_output=True)
        if ier == 1 and np.linalg.norm(g(v)) < 1e-10 and np.all(np.abs(v[:4]) < 25):
            p = np.exp(v[:4] + 1j*v[4:])
            if all(np.linalg.norm(p - s) > 1e-6*(1 + np.linalg.norm(p)) for s in sols): sols.append(p)
    return len(sols)
c = np.array([1, 1, 1, 1, 1, 1.3*np.exp(0.4j)])
for mu in [np.zeros(4), np.array([.31, -.12, .47, -.2])*1j, rng.normal(size=4) + 1j*rng.normal(size=4)]:
    print("ladder, unit coeffs, mu =", np.round(mu, 2), "-> #crit", count(lad, c, mu), flush=True)
lad7 = np.vstack([lad, (1, 0, 0, -1)])
print("ladder+X/W, mu generic ->", count(lad7, np.append(c, .7), rng.normal(size=4) + 1j*rng.normal(size=4)), flush=True)
