"""LOG 30: Bessel-regime charts of the planar mirror. Flux m_a = k_a nu* + O(1/eps), nu* = (2/eps) ln(1/eps); chiral i survives
iff |k.Q_i + b_i| = 1 (b = background in units of nu*), giving monomial w^{sign * Q_i}. Kouchnirenko: #critical points
= 4! Vol(Newton polytope) for generic coefficients when 0 is interior."""
import numpy as np, itertools
from math import factorial
from scipy.spatial import ConvexHull
CH = {'B': ((0,0,1,0), 'eta'), 'L-B': ((1,0,-1,0), None), 'T-L': ((-1,1,0,0), None), '-T': ((0,-1,0,0), None),
      'R-B': ((0,0,-1,1), None), 'T-R': ((0,1,0,-1), None), 'B-T': ((0,-1,1,0), None)}   # order (L,T,B,R)
def interior(P):
    try:
        h = ConvexHull(P)
    except Exception:
        return False, 0.0
    return bool(np.all(h.equations[:, -1] < -1e-9)), h.volume
for b in (0, 1, -1):
    charts = {}
    for k in itertools.product(range(-3, 4), repeat=4):
        mons = []
        for nm, (Q, fl) in CH.items():
            v = int(np.dot(k, Q)) + (b if fl else 0)
            if abs(v) == 1: mons.append((nm, tuple(int(np.sign(v))*np.array(Q))))
        P = np.array([m for _, m in mons], float)
        if len(P) < 5 or np.linalg.matrix_rank(P) < 4: continue
        ok, vol = interior(P)
        key = tuple(sorted(n for n, _ in mons))
        charts.setdefault(key, (k, ok, round(factorial(4)*vol), [m for m in mons]))
    print(f"background b = {b}: {len(charts)} distinct charts")
    for key, (k, ok, nv, mons) in sorted(charts.items(), key=lambda t: -len(t[0])):
        print(f"  k={k} terms={len(key)} 0-interior={ok} #crit={nv if ok else 'n/a'}: " + ", ".join(f"{n}{m}" for n, m in mons))
