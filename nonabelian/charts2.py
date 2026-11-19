"""LOG 45: Bessel-regime charts. Variables (X,Y,Z,W) = (w_B,w_L,w_R,w_T). Always active: B -> X, -T -> 1/W.
Triangle T1 = {L-B: Y/X, T-L: W/Y, B-T: X/W}, T2 = {R-B: Z/X, T-R: W/Z, B-T}. In each triangle one chiral drops.
Ladder: B-T drops (in both). Others: B-T active, one of {L-B,T-L} and one of {R-B,T-R} drop. Count critical points of
each chart's LG with masses mu.log w (mu_B = mu_T), generic coefficients 1 and q."""
import numpy as np, itertools
from scipy.optimize import fsolve
rng = np.random.default_rng(7)
mon = {'B': (1,0,0,0), 'L-B': (-1,1,0,0), 'R-B': (-1,0,1,0), 'T-L': (0,-1,0,1), 'T-R': (0,0,-1,1), '-T': (0,0,0,-1), 'B-T': (1,0,0,-1)}
q = 1.3*np.exp(0.4j)
def count(names, mu, tries=5000):
    E = np.array([mon[n] for n in names], float); c = np.array([q if n == '-T' else 1.0 for n in names])
    sols = []
    def g(v):
        l = v[:4] + 1j*v[4:]; r = E.T @ (c*np.exp(E @ l)) + mu
        return np.concatenate([r.real, r.imag])
    for _ in range(tries):
        v0 = np.concatenate([rng.normal(0, 2, 4), rng.uniform(-np.pi, np.pi, 4)])
        with np.errstate(all='ignore'):
            v, info, ier, _ = fsolve(g, v0, full_output=True)
        if ier == 1 and np.linalg.norm(g(v)) < 1e-11 and np.all(np.abs(v[:4]) < 25):
            p = np.exp(v[:4] + 1j*v[4:])
            if all(np.linalg.norm(p - s) > 1e-7*(1 + np.linalg.norm(p)) for s in sols): sols.append(p)
    return len(sols)
a, bL, bR = rng.normal(size=3)*.6 + 1j*rng.normal(size=3)*.6
mu = np.array([a, bL, bR, a])
charts = {'ladder (B-T drops)': ['B', 'L-B', 'T-L', 'R-B', 'T-R', '-T']}
for d1, d2 in itertools.product(['L-B', 'T-L'], ['R-B', 'T-R']):
    charts[f'B-T active, drop {d1},{d2}'] = [n for n in ['B', 'L-B', 'T-L', 'R-B', 'T-R', '-T', 'B-T'] if n not in (d1, d2)]
for k, v in list(charts.items())[:0]:
    print(f"{k:35s} monomials {v}: #crit = {count(v, mu)}", flush=True)
print("--- alternative chart with the reduction's coefficients: invariant c_B c_{-T}/c_BT ~ x^n/eps ~ eps^3 (ladder q = O(1)) ---")
for eps in (0.1, 0.03, 0.01):
    names = ['B', 'T-L', 'T-R', '-T', 'B-T']
    E = np.array([mon[n] for n in names], float); c = np.array([1, 1, 1, q, eps**-3])   # c_BT large <=> invariant ~ eps^3
    sols = []
    def g(v):
        l = v[:4] + 1j*v[4:]; r = E.T @ (c*np.exp(E @ l)) + mu
        return np.concatenate([r.real, r.imag])
    for _ in range(3000):
        v0 = np.concatenate([rng.normal(0, 4, 4), rng.uniform(-np.pi, np.pi, 4)])
        with np.errstate(all='ignore'):
            v, info, ier, _ = fsolve(g, v0, full_output=True)
        if ier == 1 and np.linalg.norm(g(v)) < 1e-10 and np.all(np.abs(v[:4]) < 60):
            if all(np.linalg.norm(v[:4] - s) > 1e-6 for s in sols): sols.append(v[:4])
    print(f"eps={eps}: {len(sols)} critical points; max |log|w|| = {max(np.abs(s).max() for s in sols) if sols else None:.2f}; "
          f"monomial magnitudes |B/T| at them: {[f'{np.exp(s[0]-s[3]):.2e}' for s in sols]}", flush=True)
