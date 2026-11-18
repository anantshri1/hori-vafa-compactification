import numpy as np, sys
from scipy.optimize import fsolve
def count(E, c, tries=3000, seed=0):
    E = np.array(E, float); rng = np.random.default_rng(seed); sols = []
    def g(v):
        l = v[:4] + 1j*v[4:]; t = c*np.exp(E @ l); gr = E.T @ t
        return np.concatenate([gr.real, gr.imag])
    for _ in range(tries):
        v0 = np.concatenate([rng.normal(0, 1.5, 4), rng.uniform(-np.pi, np.pi, 4)])
        with np.errstate(all='ignore'):
            v, info, ier, _ = fsolve(g, v0, full_output=True)
        if ier == 1 and np.linalg.norm(g(v)) < 1e-9 and np.all(np.abs(v[:4]) < 30):
            p = np.exp(v[:4] + 1j*v[4:])
            if all(np.linalg.norm(p - s) > 1e-5*(1+np.linalg.norm(p)) for s in sols): sols.append(p)
    return len(sols)
lad = [(1,0,0,0),(-1,1,0,0),(-1,0,1,0),(0,-1,0,1),(0,0,-1,1),(0,0,0,-1)]          # (X,Y,Z,W) ~ (B,L,R,T)
# chart b=0 k=(-2,-1,-1,-2) in (L,T,B,R): 1/B, B/L, T/L, 1/T, B/R, T/R
ch0 = [(0,0,-1,0),(-1,0,1,0),(-1,1,0,0),(0,-1,0,0),(0,0,1,-1),(0,1,0,-1)]
rng = np.random.default_rng(5)
for name, E in [("ladder", lad), ("ladder+X/W", lad + [(1,0,0,-1)]), ("chart b=0 6-term", ch0)]:
    for trial in range(2):
        c = rng.normal(size=len(E)) + 1j*rng.normal(size=len(E))
        print(name, "generic coeffs -> #crit", count(E, c, seed=trial), flush=True)
print("ladder unit coeffs (q=1.3e^{0.4i}) ->", count(lad, np.array([1,1,1,1,1,1.3*np.exp(0.4j)])))
