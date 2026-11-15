"""Electric side vs the thesis reference index (Unitary_Mirrors_1.tex, Table tab:SQCD_index): N=2, F=4, r=1/5, chirals at
trial R = 1 - r, flavour fugacities 1:  1 + eta x^{4/5} + eta^2 x^{8/5} - 16 x^2 + eta^3 x^{12/5} + eta^4 x^{16/5} + (eta^5+88) x^4 + O(x^{21/5})."""
import numpy as np
from idx import theory_index
r = 0.2
def table(x, e):
    return 1 + e*x**0.8 + e**2*x**1.6 - 16*x**2 + e**3*x**2.4 + e**4*x**3.2 + (e**5 + 88)*x**4
def electric(x, eta, fi_sign):
    ch = [((1, 0), {}, 1-r)]*4 + [((0, 1), {}, 1-r)]*4
    return theory_index(2, ch, [[-1, -1], [-1, -1]], {}, {0: ['w'], 1: ['w']}, {'w': eta**fi_sign}, x, M=10, Ngrid=48, vector=2, weyl=2)
for x in (0.05, 0.1):
    for eta in (np.exp(0.7j), np.exp(-1.9j)):
        T = table(x, eta)
        for s in (1, -1):
            E = electric(x, eta, s)
            print(f"x={x} eta={eta:.3f} FI sign {s:+d}: E={E:.8f} table={T:.8f} |diff|={abs(E-T):.2e}", flush=True)
