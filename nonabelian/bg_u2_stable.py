"""Flux-2 sector map in a numerically stable regime (r=0.45, x=0.03): electric with M=6, Ng=96 (M=5,6 agree to 3e-5;
M>=7 is destroyed by cancellation), mirror with M=2 and 3. Prints C = ratio/(w (y1 y2)^2) at two points."""
import numpy as np, refined as RF, idx
from idx import theory_index
orig = idx.theory_index
RF.x, RF.r = 0.03, 0.45; x, r = RF.x, RF.r
def el(y, e):
    ch = [((1, 0), {f'y{j}': 1}, 1-r) for j in range(4)] + [((0, 1), {f'y{j}': 1}, 1-r) for j in range(4)]
    return theory_index(2, ch, [[-1, -1], [-1, -1]], {}, {0: ['w'], 1: ['w']}, {**{f'y{j}': y[j] for j in range(4)}, 'w': e*x**(2*r)},
                        x, M=6, Ngrid=96, vector=2, weyl=2, bg={'w': 2})
def mir(y, e, M):
    def patched(*a, **k): k['bg'] = {'eta': -2}; return orig(*a, **k)
    RF.theory_index = patched
    v = RF.mirror(y, e, 1, M=M, Ng=20); RF.theory_index = orig; return v
for yang, wang in [([0.4, 1.1, -0.6], 0.9), ([-0.8, 0.5, 1.7], -0.7)]:
    y3 = np.exp(1j*np.array(yang)); y = list(y3) + [1/np.prod(y3)]; e = np.exp(1j*wang)
    E = el(y, e)
    for M in (2, 3):
        print(M, 'C =', E/mir(y, e, M)/(e*(y[0]*y[1])**2), flush=True)
