"""Monomial check of the flux-2 sector map with the mirror flux range enlarged (M=3): ratio / (w (y1 y2)^2) at two points."""
import numpy as np, refined as RF, bg_u2 as G
import idx
orig = idx.theory_index
def mir(y, e, nb, M):
    def patched(*a, **k): k['bg'] = {'eta': nb}; return orig(*a, **k)
    RF.theory_index = patched
    v = RF.mirror(y, e, 1, M=M, Ng=20); RF.theory_index = orig; return v
for yang, wang in [([0.4, 1.1, -0.6], 0.9), ([-0.8, 0.5, 1.7], -0.7)]:
    y3 = np.exp(1j*np.array(yang)); y = list(y3) + [1/np.prod(y3)]; e = np.exp(1j*wang)
    E = G.el(y, e, 2)
    for M in (2, 3):
        R = E/mir(y, e, -2, M); print(M, R/(e*(y[0]*y[1])**2), flush=True)
