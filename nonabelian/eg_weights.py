"""Tangent weights (additive, in units where X = e^{2 pi i w}) at fixed points: toric (4 non-pole chirals) and Gr(2,4)."""
import numpy as np, itertools, eg as G
M = G.M; mu = G.mu; nm = ['X','Y/X','Z/X','W/Y','W/Z','q/W']
m = np.linalg.lstsq(G.E.T, mu, rcond=None)[0]; eta = np.array([1.0, 0.37])
def red(w): return np.round(w - np.round(w.real), 4)
print("toric fixed points (weights mod 1):")
for P in itertools.combinations(range(6), 2):
    QP = G.Qg[:, list(P)]; d = np.linalg.det(QP)
    if abs(d) < 1e-9 or np.any(np.linalg.solve(QP, eta) < 0): continue
    s = np.linalg.solve(QP.T, -m[list(P)])
    w = sorted([red(G.Qg[:, i] @ s + m[i]) for i in range(6) if i not in P], key=lambda z: (z.real, z.imag))
    print(" ", (nm[P[0]], nm[P[1]]), w)
print("Gr(2,4) vacua (tangent weights M_k - M_j, j in vacuum, k not):")
for j1, j2 in itertools.combinations(range(4), 2):
    w = sorted([red(M[k] - M[j]) for j in (j1, j2) for k in range(4) if k not in (j1, j2)], key=lambda z: (z.real, z.imag))
    print(" ", (j1, j2), w)
