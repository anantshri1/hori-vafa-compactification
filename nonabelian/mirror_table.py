"""Planar mirror of U(2)_{(0,-2)} with [4,0], built from the thesis figure (ChirPlan_SQCDdual) and its caption:
nodes L (column 1), T and B (column 2, top/bottom), R (column 3); diagonal chirals R=r, vertical chiral R=2-2r;
self CS +1, mixed CS -1/2 on diagonal links and +1/2 on the vertical link; flavour eta attached to B, the other flavour to T;
column-2 FIs carry the R-mixing +-iQ/4 (here an x-power +-c on the FI fugacity, scanned). Orientation from the closed
triangles and the meson path eta -> B -> L -> T -> flavour. Compared with the thesis reference index (unrefined, r=1/5)."""
import numpy as np, itertools, sys
from idx import theory_index
r = 0.2
def table(x, e): return 1 + e*x**0.8 + e**2*x**1.6 - 16*x**2 + e**3*x**2.4 + e**4*x**3.2 + (e**5 + 88)*x**4
L, T, B, R = 0, 1, 2, 3
def mirror(x, eta, o, c, t, M=2, Ng=12):
    def f(head, tail): Q = [0]*4; Q[head] += o; Q[tail] -= o; return tuple(Q)
    def one(node, sgn): Q = [0]*4; Q[node] = o*sgn; return tuple(Q)
    ch = [(one(B, 1), {'eta': -o}, r), (f(L, B), {}, r), (f(T, L), {}, r), (one(T, -1), {}, r),
          (f(R, B), {}, r), (f(T, R), {}, r), (f(B, T), {}, 2 - 2*r)]
    K = np.eye(4).tolist()
    for a, b, v in [(L, B, -.5), (L, T, -.5), (R, B, -.5), (R, T, -.5), (T, B, .5)]:
        K[a][b] += v; K[b][a] += v
    fl = {'eta': eta*x**t, 'fT': x**c, 'fB': x**(-c), 'one': 1.0}
    return theory_index(4, ch, K, {(B, 'eta'): -.5}, {L: ['one'], T: ['fT'], B: ['fB'], R: ['one']}, fl, x, M=M, Ngrid=Ng)
if __name__ == "__main__":
    x = 0.05
    res = []
    for o, es, c, t in itertools.product((1, -1), (1, -1), (0.5, -0.5, 1, -1, 0), (0,)):
        err = max(abs(mirror(x, e**es, o, c, t) - table(x, e)) for e in (np.exp(0.7j), np.exp(-1.9j)))
        res.append((err, o, es, c, t)); print(f"o={o:+d} eta^{es:+d} c={c} t={t}: max|M-table| = {err:.2e}", flush=True)
