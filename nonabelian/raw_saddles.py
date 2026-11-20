"""LOG 46: self-consistent saddles of the piecewise raw Bessel-regime W (all 2^7 branch patterns).
Chiral i, Y_i = log(w^{Q_i}) + log(x^{bg_i}) (B carries the eta background n with x^n = eps^4):
  branch + (|e^{Y_i}| < 1, nu > 0):  T_i = x^R f_i e^{Y_i}/eps
  branch - (|e^{Y_i}| > 1, nu < 0):  T_i = -x^{2-R} e^{-Y_i}/(f_i eps)
W = sum_i T_i + mu.log w ; critical: sum_i Q_ia (+-) T_i + mu_a = 0. Keep a solution only if every |e^{Y_i}| is on the
side its branch assumes. Count all consistent saddles; expect 6, all in the ladder pattern."""
import numpy as np, itertools
from scipy.optimize import fsolve
Q = np.array([(0,0,1,0),(1,0,-1,0),(-1,1,0,0),(0,-1,0,0),(0,0,-1,1),(0,1,0,-1),(0,-1,1,0)], float)   # (L,T,B,R)
nm = ['B','L-B','T-L','-T','R-B','T-R','B-T']; R = np.array([.3]*6 + [1.4])
rng = np.random.default_rng(2)
ph = np.exp(1j*rng.uniform(-np.pi, np.pi, 7)); ph[6] = 1/(ph[1]*ph[2]); ph[5] = 1/(ph[4]*ph[6])   # triangle constraints
M = rng.normal(0, .3, 4); mu = np.array([M[1]-M[0], M[2]-M[1], M[2]-M[1], M[3]-M[2]])[[0, 1, 2, 3]]  # (L,T,B,R): T=B column
mu = np.array([M[1]-M[0], M[2]-M[1], M[2]-M[1], M[3]-M[2]])
for eps in (0.1, 0.03):
    x = np.exp(-eps/2); bg = np.zeros(7); bg[0] = 4*np.log(eps)   # log x^n = 4 log eps
    found = []
    for br in itertools.product((1, -1), repeat=7):
        br = np.array(br)
        def T(l):
            Y = Q @ l + bg
            return np.where(br > 0, x**R*ph*np.exp(Y)/eps, -x**(2-R)*np.exp(-Y)/(ph*eps))
        def g(v):
            l = v[:4] + 1j*v[4:]; r = Q.T @ (br*T(l)) + mu
            return np.concatenate([r.real, r.imag])
        for _ in range(150):
            v0 = np.concatenate([rng.normal(np.log(eps)*rng.uniform(0, 4, 4), 2), rng.uniform(-np.pi, np.pi, 4)])
            with np.errstate(all='ignore'):
                v, info, ier, _ = fsolve(g, v0, full_output=True)
            if ier != 1 or np.linalg.norm(g(v)) > 1e-9*(1 + np.abs(T(v[:4]+1j*v[4:])).max()): continue
            Y = Q @ v[:4] + bg.real
            if np.all(np.where(br > 0, Y < 0, Y > 0)):
                if all(np.linalg.norm(v[:4] - f[1]) > 1e-6 for f in found):
                    found.append((tuple(br), v[:4].copy()))
    pats = {}
    for b, _ in found: pats[b] = pats.get(b, 0) + 1
    print(f"eps={eps}: {len(found)} self-consistent saddles; by branch pattern (B,L-B,T-L,-T,R-B,T-R,B-T):", pats, flush=True)
