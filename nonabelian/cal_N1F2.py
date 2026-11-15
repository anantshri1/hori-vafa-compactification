"""Calibration: thesis planar-mirror data at N=1, F=2 (Unitary_Mirrors_1.tex eq. mirpartfun) translated to the index
engine, against the electric U(1) with [2,0] (eq. u_n_chiral). Scans the discrete translation choices (FI sign, R-mixing
powers) and reports the best; the abelian CS pair is already verified, so this fixes the dictionary for N=2."""
import numpy as np, itertools
from idx import theory_index
x, r = 0.3, 0.3
pts = [(np.exp(0.4j), np.exp(0.9j)), (np.exp(1.3j), np.exp(-0.7j)), (np.exp(2.2j), np.exp(1.6j))]  # (y, w)
def el(y, w, p):
    return theory_index(1, [((1,), {'y': 1}, r), ((1,), {'y': -1}, r)], [[-1]], {}, {0: ['om']},
                        {'y': y, 'om': w*x**p}, x, M=14, Ngrid=64)
def mir(y, w, t, s):
    return theory_index(1, [((1,), {'eta': -1}, 0), ((-1,), {}, 0)], [[1]], {(0, 'eta'): -0.5}, {0: ['fi']},
                        {'eta': w*x**t, 'fi': y**(2*s)}, x, M=14, Ngrid=64)
res = []
for p in [0, 1-r, -(1-r), (1-r)/2, -(1-r)/2]:
    E = [el(y, w, p) for y, w in pts]
    for t in [0, 1-r, -(1-r), (1-r)/2, -(1-r)/2, r, -r]:
        for s in (1, -1):
            for we in (1, -1):
                Mv = [mir(y, w**we, t, s) for y, w in pts]
                rat = [a/b for a, b in zip(E, Mv)]
                spread = max(abs(z/rat[0]-1) for z in rat)
                res.append((spread, p, t, s, we, rat[0]))
res.sort(key=lambda z: z[0])
for z in res[:5]: print(f"spread {z[0]:.1e}: electric omega=w x^{z[1]:.2f}; mirror eta-fug=w^{z[4]} x^{z[2]:.2f}, FI=y^{2*z[3]}; ratio {z[5]:.6f}")
