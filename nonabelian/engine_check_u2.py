"""Engine check for the non-abelian (U(2)) part: Aharony duality U(2)_0 with 2 flavours <-> WZ model of mesons M_ij and
monopoles V_+-: M fugacity a^2 y_i yt_j, R=2r;  V_+- fugacity a^{-2} w^{+-1}, R = 2(1-r) - 1."""
import numpy as np
from idx import theory_index, qpoch
x, r = 0.15, 0.3
def free(zeta, R):
    return qpoch(x**(2-R)/zeta, x*x)/qpoch(zeta*x**R, x*x)
for a, y, yt, w in [(np.exp(0.3j), np.exp(0.5j), np.exp(-0.4j), np.exp(1.1j)), (np.exp(-0.7j), np.exp(1.2j), np.exp(0.6j), np.exp(-0.3j))]:
    ys = [y, 1/y]; yts = [yt, 1/yt]
    ch = []
    for i in range(2):
        for c in range(2):
            Q = [0, 0]; Q[c] = 1; ch.append((tuple(Q), {'a': 1, f'y{i}': 1}, r))
            Q = [0, 0]; Q[c] = -1; ch.append((tuple(Q), {'a': 1, f't{i}': 1}, r))
    fl = {'a': a, 'y0': ys[0], 'y1': ys[1], 't0': yts[0], 't1': yts[1], 'w': w}
    E = theory_index(2, ch, [[0, 0], [0, 0]], {}, {0: ['w'], 1: ['w']}, fl, x, M=12, Ngrid=64, vector=2, weyl=2)
    D = np.prod([free(a*a*yi*tj, 2*r) for yi in ys for tj in yts]) * free(a**-2*w, 1-2*r) * free(a**-2/w, 1-2*r)
    print(E, D, abs(E/D-1))
