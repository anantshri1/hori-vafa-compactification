"""LOG 32: normalisation-free test. Gr(2,4) equivariant QC: q dW/dq at the 6 vacua = -(r_i + r_j), r = 4 roots of
prod(sigma - m_j) = +-q. Pair sums of 4 numbers <=> the 6 values split into 3 complementary pairs with EQUAL sums.
Test on the LG side (ladder + mu.log w): t = value of the q/W term at each critical point. Also: the 4 numbers a_i are then
recoverable and must satisfy prod(a - m'_j) = const for some masses m' (affine in mu)."""
import numpy as np, itertools
from scipy.optimize import fsolve
rng = np.random.default_rng(4)
lad = np.array([(1,0,0,0),(-1,1,0,0),(-1,0,1,0),(0,-1,0,1),(0,0,-1,1),(0,0,0,-1)], float)
def crit(c, mu, tries=6000):
    out = []
    def g(v):
        l = v[:4] + 1j*v[4:]; r = lad.T @ (c*np.exp(lad @ l)) + mu
        return np.concatenate([r.real, r.imag])
    for _ in range(tries):
        v0 = np.concatenate([rng.normal(0, 2, 4), rng.uniform(-np.pi, np.pi, 4)])
        with np.errstate(all='ignore'):
            v, info, ier, _ = fsolve(g, v0, full_output=True)
        if ier == 1 and np.linalg.norm(g(v)) < 1e-12 and np.all(np.abs(v[:4]) < 25):
            p = np.exp(v[:4] + 1j*v[4:])
            if all(np.linalg.norm(p - s) > 1e-7*(1 + np.linalg.norm(p)) for s in out): out.append(p)
    return out
def best_pairing(t):
    best = (np.inf, None)
    idx = list(range(6))
    for p in [((0,a),) for a in range(1,6)]:
        rest = [i for i in idx if i not in p[0]]
        for b in rest[1:]:
            pr2 = (rest[0], b); rest2 = [i for i in rest if i not in pr2]
            pairs = [p[0], pr2, tuple(rest2)]
            s = [t[i] + t[j] for i, j in pairs]
            err = (abs(s[0]-s[1]) + abs(s[1]-s[2]))/np.mean(np.abs(t))
            if err < best[0]: best = (err, pairs, np.mean(s))
    return best
q = 1.3*np.exp(0.4j); c = np.array([1, 1, 1, 1, 1, q])
for trial in range(0):
    mu = rng.normal(size=4)*0.6 + 1j*rng.normal(size=4)*0.6
    P = crit(c, mu)
    t = np.array([q*np.exp(lad[5] @ np.log(p)) for p in P])
    err, pairs, S = best_pairing(t) if len(t) == 6 else (None, None, None)
    print(f"mu={np.round(mu,2)}: #crit={len(P)}  pair-sum test rel.err={err:.2e}" if err is not None else f"#crit={len(P)}", flush=True)
    if err is not None and err < 1e-6:
        # recover a_i: a1+a2 = t_p1, a3+a4 = S - ..., solve 4 numbers from the 6 pair sums
        A = []; T = {}
        for (i, j) in pairs: pass
        print("   complementary pairs:", pairs, " common sum S =", np.round(S, 6))
print("--- constrained masses (B and T in the same column) ---", flush=True)
for rel in [(1, 1), (1, -1), (1, 0), (0, 1)]:
    for trial in range(2):
        a, bL, bR = rng.normal(size=3)*0.6 + 1j*rng.normal(size=3)*0.6
        mu = np.array([rel[0]*a, bL, bR, rel[1]*a])   # (B, L, R, T)
        P = crit(c, mu, tries=4000)
        t = np.array([q*np.exp(lad[5] @ np.log(p)) for p in P])
        if len(t) == 6:
            err, pairs, S = best_pairing(t); print(f"mu_B={rel[0]}a, mu_T={rel[1]}a: #crit 6, pair-sum rel.err {err:.2e}", flush=True)
        else: print(f"mu_B={rel[0]}a, mu_T={rel[1]}a: #crit {len(t)}", flush=True)
