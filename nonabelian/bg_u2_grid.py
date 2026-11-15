"""Settle the 1% residual of the flux-2 sector map: same two points, finer grids (mirror Ng=26, electric Ng=56)."""
import numpy as np, refined as RF, bg_u2 as G, idx
from idx import theory_index
orig = idx.theory_index
x, r = G.x, G.r
def el(y, e, n, Ng):
    ch = [((1, 0), {f'y{j}': 1}, 1-r) for j in range(4)] + [((0, 1), {f'y{j}': 1}, 1-r) for j in range(4)]
    return theory_index(2, ch, [[-1, -1], [-1, -1]], {}, {0: ['w'], 1: ['w']},
                        {**{f'y{j}': y[j] for j in range(4)}, 'w': e*x**(2*r)}, x, M=10, Ngrid=Ng, vector=2, weyl=2, bg={'w': n})
def mir(y, e, nb, Ng):
    def patched(*a, **k): k['bg'] = {'eta': nb}; return orig(*a, **k)
    RF.theory_index = patched
    v = RF.mirror(y, e, 1, M=2, Ng=Ng); RF.theory_index = orig; return v
if __name__ == "__main__":
   for yang, wang in [([0.4, 1.1, -0.6], 0.9), ([-0.8, 0.5, 1.7], -0.7)]:
      y3 = np.exp(1j*np.array(yang)); y = list(y3) + [1/np.prod(y3)]; e = np.exp(1j*wang)
      E = el(y, e, 2, 56)
      R = E/mir(y, e, -2, 26); print('C =', R/(e*(y[0]*y[1])**2), flush=True)
