"""Fix the U(2) CS conventions against the thesis reference index (N=2, F=4, r=1/5): scan the sign of the CS matrix,
the flux-dependent sign factor and the FI orientation. The zero-flux sector already matches (1 - 16x^2 + ...)."""
import numpy as np, itertools
from idx import theory_index
r = 0.2
def table(x, e): return 1 + e*x**0.8 + e**2*x**1.6 - 16*x**2 + e**3*x**2.4 + e**4*x**3.2 + (e**5 + 88)*x**4
ch = [((1, 0), {}, 1-r)]*4 + [((0, 1), {}, 1-r)]*4
signs = {'1': lambda m: 1, '(-1)^sum m': lambda m: (-1)**int(m[0]+m[1]), '(-1)^m1m2': lambda m: (-1)**int(m[0]*m[1]),
         '(-1)^(sum m)^2/..': lambda m: (-1)**int((m[0]+m[1])*(m[0]+m[1]+1)//2)}
res = []
for l, (sn, sf), fs in itertools.product((1, -1), signs.items(), (1, -1)):
    err = 0
    for x in (0.05, 0.1):
        for eta in (np.exp(0.7j), np.exp(-1.9j)):
            E = theory_index(2, ch, [[l, l], [l, l]], {}, {0: ['w'], 1: ['w']}, {'w': eta**fs}, x, M=4, Ngrid=48, vector=2, weyl=2, signfn=sf)
            err = max(err, abs(E - table(x, eta)))
    res.append((err, l, sn, fs)); print(f"K={l:+d}*J  sign {sn:18s} FI^{fs:+d}: max|E-table| = {err:.2e}", flush=True)
