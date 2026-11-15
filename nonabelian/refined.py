"""Refined test of U(2)_{(0,-2)} [4,0] <-> planar mirror: electric with flavour fugacities y_1..y_4 (prod = 1),
FI = eta x^{2r}; mirror with column FIs (y_{I+1}/y_I)^s, column-2 R-mixing x^{+-r/2}, eta-fugacity eta^{-1}."""
import numpy as np, sys
from idx import theory_index
r = 0.3; x = float(sys.argv[1]) if len(sys.argv) > 1 else 0.05
L, T, B, R = 0, 1, 2, 3
def electric(y, e):
    ch = [((1, 0), {f'y{j}': 1}, 1-r) for j in range(4)] + [((0, 1), {f'y{j}': 1}, 1-r) for j in range(4)]
    return theory_index(2, ch, [[-1, -1], [-1, -1]], {}, {0: ['w'], 1: ['w']},
                        {**{f'y{j}': y[j] for j in range(4)}, 'w': e*x**(2*r)}, x, M=8, Ngrid=40, vector=2, weyl=2)
def mirror(y, e, s, M=int(sys.argv[2]) if len(sys.argv) > 2 else 2, Ng=int(sys.argv[3]) if len(sys.argv) > 3 else 20):
    def f(head, tail): Q = [0]*4; Q[head] += 1; Q[tail] -= 1; return tuple(Q)
    def one(node, sgn): Q = [0]*4; Q[node] = sgn; return tuple(Q)
    ch = [(one(B, 1), {'eta': -1}, r), (f(L, B), {}, r), (f(T, L), {}, r), (one(T, -1), {}, r),
          (f(R, B), {}, r), (f(T, R), {}, r), (f(B, T), {}, 2 - 2*r)]
    K = np.eye(4).tolist()
    for a, b, v in [(L, B, -.5), (L, T, -.5), (R, B, -.5), (R, T, -.5), (T, B, .5)]:
        K[a][b] += v; K[b][a] += v
    rat = [(y[i+1]/y[i])**s for i in range(3)]
    fl = {'eta': 1/e, 'f1': rat[0], 'fT': rat[1]*x**(r/2), 'fB': rat[1]*x**(-r/2), 'f3': rat[2]}
    return theory_index(4, ch, K, {(B, 'eta'): -.5}, {L: ['f1'], T: ['fT'], B: ['fB'], R: ['f3']}, fl, x, M=M, Ngrid=Ng)
if __name__ == "__main__":
    for y3, e in [((np.exp(0.4j), np.exp(1.1j), np.exp(-0.6j)), np.exp(0.9j)), ((np.exp(-0.8j), np.exp(0.5j), np.exp(1.7j)), np.exp(-0.7j))]:
        y = list(y3) + [1/np.prod(y3)]
        E = electric(y, e)
        for s in (1, -1):
            m = mirror(y, e, s)
            print(f"e={e:.3f} s={s:+d}: E={E:.8f} M={m:.8f} |E/M-1|={abs(E/m-1):.2e}", flush=True)
