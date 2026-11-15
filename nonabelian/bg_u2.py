"""Sector map with background topological flux n for U(2)_(0,-2)[4,0] (electric bg on 'w') vs the planar mirror
(bg on 'eta', value scanned). A correct map gives a ratio that is a monomial in the fugacities; we test that the ratio at
two fugacity points differs only by (w1/w2)^a with a in (1/4)Z, and report |ratio| and the best a."""
import numpy as np, refined as RF
RF.x = 0.05; x = RF.x; r = RF.r
from idx import theory_index
pts = [((np.exp(0.4j), np.exp(1.1j), np.exp(-0.6j)), np.exp(0.9j)), ((np.exp(-0.8j), np.exp(0.5j), np.exp(1.7j)), np.exp(-0.7j))]
def el(y, e, n):
    ch = [((1, 0), {f'y{j}': 1}, 1-r) for j in range(4)] + [((0, 1), {f'y{j}': 1}, 1-r) for j in range(4)]
    return theory_index(2, ch, [[-1, -1], [-1, -1]], {}, {0: ['w'], 1: ['w']},
                        {**{f'y{j}': y[j] for j in range(4)}, 'w': e*x**(2*r)}, x, M=8, Ngrid=40, vector=2, weyl=2, bg={'w': n})
import idx
orig = idx.theory_index
def mir(y, e, nb):
    def patched(*a, **k): k['bg'] = {'eta': nb}; return orig(*a, **k)
    RF.theory_index = patched
    v = RF.mirror(y, e, 1, M=2, Ng=20); RF.theory_index = orig; return v
for n in (2,):
    E = [el(list(y3) + [1/np.prod(y3)], e, n) for y3, e in pts]
    print('E', E, flush=True)
    for nb in (n, -n, n//2, -n//2):
        Mv = [mir(list(y3) + [1/np.prod(y3)], e, nb) for y3, e in pts]
        rat = [a/b for a, b in zip(E, Mv)]
        q = rat[0]/rat[1]; w1, w2 = pts[0][1], pts[1][1]
        best = min(((abs(q/(w1/w2)**(k/4) - 1), k/4) for k in range(-16, 17)))
        print(f"n={n} eta-flux {nb:+d}: ratios {rat[0]:.6f} {rat[1]:.6f}; best w-power {best[1]} resid {best[0]:.2e}", flush=True)
