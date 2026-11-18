"""LOG 31: Bethe-vacuum count (3d Witten index) in the symmetric quantisation: chiral factor (X^{1/2}/(1-X))^{Q_a},
CS z^{K}, FI. Electric U(2) K = -J, 4 fundamentals; mirror: abelian 4-node quiver (refined.py data). Generic fugacities."""
import numpy as np
from scipy.optimize import fsolve
rng = np.random.default_rng(3)
def solve(F, nvar, tries, scale=2.0):
    sols = []
    for _ in range(tries):
        v0 = np.concatenate([rng.normal(0, scale, nvar), rng.uniform(-np.pi, np.pi, nvar)])
        g = lambda v: (lambda r: np.concatenate([r.real, r.imag]))(F(v[:nvar] + 1j*v[nvar:]))
        with np.errstate(all='ignore'):
            v, info, ier, _ = fsolve(g, v0, full_output=True)
        if ier == 1 and np.linalg.norm(g(v)) < 1e-9 and np.all(np.abs(v[:nvar]) < 25):
            z = np.exp(v[:nvar] + 1j*v[nvar:])
            sols.append(z)
    return sols
def dedupe(sols, key):
    out = []
    for s in sols:
        k = key(s)
        if all(np.linalg.norm(k - o) > 1e-6*(1 + np.linalg.norm(k)) for o in out): out.append(k)
    return out
y = np.exp(rng.normal(0, .3, 4) + 1j*rng.uniform(-np.pi, np.pi, 4)); y /= np.prod(y)**.25
w = 0.7*np.exp(0.9j)
def F_el(l):   # log z_a
    z = np.exp(l); s = z[0]*z[1]
    return np.array([np.log(-(z[a]**2/s)/w) - np.sum(np.log(1 - y*z[a])) for a in range(2)])  # mod 2 pi i handled below
def F_el_poly(l):
    z = np.exp(l); s = z[0]*z[1]
    return np.array([-(z[a]**2) - w*s*np.prod(1 - y*z[a]) for a in range(2)])/(1+np.abs(z)**4)
sols = solve(F_el_poly, 2, 4000)
sols = [s for s in sols if abs(s[0] - s[1]) > 1e-5*(1 + abs(s[0]))]          # Weyl: exclude z1 = z2
el = dedupe(sols, lambda s: np.sort_complex(s))
print("electric U(2)_{(0,-2)}[4,0] Bethe vacua (unordered, z1 != z2):", len(el), flush=True)
# ---- mirror: abelian quiver, K' = K + 1/2 sum_i Q_i Q_i^T (integer by parity), chiral factor (1 - X_i)^{-Q_ia}
import optc_defs as O
Qs = np.array([c[0] for c in [(O.one(O.B, 1), {}, 0)] + O.rest], float)
K = np.array(O.Kmat(1), float); Kp = K + 0.5*Qs.T @ Qs
assert np.allclose(Kp, np.round(Kp)), Kp
print("mirror K' =\n", Kp)
f = np.exp(rng.normal(0, .3, len(Qs)) + 1j*rng.uniform(-np.pi, np.pi, len(Qs)))   # generic chiral fugacities
t = np.exp(rng.normal(0, .3, 4) + 1j*rng.uniform(-np.pi, np.pi, 4))               # generic FIs
def F_mir(l):
    z = np.exp(l); X = f*np.exp(Qs @ l)
    lhs = np.exp(Kp @ l)*np.prod((1 - X)[:, None]**(-Qs), axis=0)
    return (lhs - t)/(1 + np.abs(lhs))
ms = dedupe(solve(F_mir, 4, 6000), lambda s: s)
print("mirror Bethe vacua:", len(ms), flush=True)
