"""Compare the 6 critical values of ladder + mu.log w with Gr(2,4) equivariant GLSM critical values (several sign/map
conventions tried; report the best pairwise-difference match)."""
import numpy as np, itertools
from scipy.optimize import fsolve
rng = np.random.default_rng(2)
lad = np.array([(1,0,0,0),(-1,1,0,0),(-1,0,1,0),(0,-1,0,1),(0,0,-1,1),(0,0,0,-1)], float)  # vars (X,Y,Z,W)=(B,L,R,T)
def crit(c, mu, tries=5000):
    out = []
    def g(v):
        l = v[:4] + 1j*v[4:]; r = lad.T @ (c*np.exp(lad @ l)) + mu
        return np.concatenate([r.real, r.imag])
    for _ in range(tries):
        v0 = np.concatenate([rng.normal(0, 2, 4), rng.uniform(-np.pi, np.pi, 4)])
        with np.errstate(all='ignore'):
            v, info, ier, _ = fsolve(g, v0, full_output=True)
        if ier == 1 and np.linalg.norm(g(v)) < 1e-11 and np.all(np.abs(v[:4]) < 25):
            l = v[:4] + 1j*v[4:]
            if all(np.linalg.norm(np.exp(l) - np.exp(s)) > 1e-6*(1 + np.linalg.norm(np.exp(l))) for s in out): out.append(l)
    return out
m = np.array([0.31, -0.12, 0.47, -0.66]) + 0.05j*np.array([1, -2, 0.5, 0.3]); q = 1.3*np.exp(0.4j)
c = np.array([1, 1, 1, 1, 1, q])
def glsm_vals(m, q, s):
    r = np.roots(np.poly(m) - np.array([0, 0, 0, 0, s*q]))
    vals = []
    for i, j in itertools.combinations(range(4), 2):
        sig = np.array([r[i], r[j]])
        vals.append(sum(-np.log(s*q)*sa + sum((sa - mj)*(np.log(sa - mj) - 1) for mj in m) for sa in sig))
    return np.array(vals), r
best = (np.inf, None)
for perm in itertools.permutations(range(4)):
    mp = m[list(perm)]
    for sgn in (+1, -1):
        mu = sgn*np.array([mp[1]-mp[0], mp[2]-mp[1], mp[2]-mp[1], mp[3]-mp[2]])  # (B? ) assigned below
        for order in [(0,1,2,3)]:
            # node order in lad vars: X=B, Y=L, Z=R, W=T ; column FIs: L: m2-m1, T,B: m3-m2, R: m4-m3
            muv = sgn*np.array([mp[2]-mp[1], mp[1]-mp[0], mp[3]-mp[2], mp[2]-mp[1]])
            L = crit(c, muv, tries=1500)
            if len(L) != 6: continue
            Wl = np.array([np.sum(c*np.exp(lad @ l)) + muv @ l for l in L])
            for s in (+1, -1):
                Wg, _ = glsm_vals(mp, q, s)
                # compare sorted pairwise derivative-free invariant: W values mod 2 pi i (log branches) -> use exp-free:
                # use the multiset of e^{W} differences is branch-sensitive; compare d W / d q numerically instead:
                dq = 1e-6
                Wl2 = []
                for l in crit(c*np.array([1,1,1,1,1,1+dq]), muv, tries=1500):
                    Wl2.append(np.sum(c*np.array([1,1,1,1,1,1+dq])*np.exp(lad @ l)) + muv @ l)
                if len(Wl2) != 6: continue
                # q dW/dq on LG side = term q/W at the critical point (exact): use it directly
                tl = np.sort_complex(np.array([c[5]*np.exp(lad[5] @ l) for l in L]))
                Wgp, r = glsm_vals(mp, q*(1+dq), s); tg = np.sort_complex(-(Wgp - Wg)/dq)
                err = np.min([np.linalg.norm(np.sort_complex(k*tl) - tg) for k in (1, -1, 2, -2, .5, -.5)])
                if err < best[0]: best = (err, (perm, sgn, s))
print("best match of q dW/dq spectra (6 vacua):", best)
