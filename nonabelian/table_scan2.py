"""Add a gauge-R mixing x^{c sum m} (equivalently an x-power on the FI fugacity) to the U(2) electric side and scan c."""
import numpy as np, itertools
from idx import theory_index
r = 0.2
def table(x, e): return 1 + e*x**0.8 + e**2*x**1.6 - 16*x**2 + e**3*x**2.4 + e**4*x**3.2 + (e**5 + 88)*x**4
ch = [((1, 0), {}, 1-r)]*4 + [((0, 1), {}, 1-r)]*4
res = []
for l, fs, c in itertools.product((1, -1), (1, -1), np.round(np.arange(-2.0, 2.01, 0.1), 2)):
    err = 0
    for x in (0.05,):
        for eta in (np.exp(0.7j), np.exp(-1.9j)):
            E = theory_index(2, ch, [[l, l], [l, l]], {}, {0: ['w'], 1: ['w']}, {'w': eta**fs*x**c}, x, M=4, Ngrid=40, vector=2, weyl=2)
            err = max(err, abs(E - table(x, eta)))
    res.append((err, l, fs, c))
res.sort()
for z in res[:6]: print(f"err {z[0]:.2e}: K={z[1]:+d}J, FI=eta^{z[2]:+d} x^{z[3]}")
