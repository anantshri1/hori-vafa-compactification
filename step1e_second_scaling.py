"""
Step 1e: second scaling. Topological background flux n = round(kappa/eps), w = e^{i theta} fixed.
Exact identity (step 1d):  I_E(w, n) = (w x^{-R})^{n/2} Chir(w, n, R_T).
Question: in this scaling, where does the electric side's weight sit (fixed flux m? which m?)
and what does the dual reduce to?
"""
import sys, numpy as np
sys.path.insert(0, '.')
from index_lib import chir_ratio, free_chiral
import mpmath as mp

R, theta, kappa = 0.2, 0.7, 1.0
RT = (1 - R) / 2
w = np.exp(1j * theta)

for eps in [0.4, 0.2, 0.1, 0.05]:
    x = np.exp(-eps / 2); q = x * x
    n = int(round(kappa / eps))
    N = 8192
    z = np.exp(2j * np.pi * (np.arange(N) + 0.5) / N)
    M = int(40 / eps)
    S = {}
    for m in range(-M, M + 1):
        S[m] = w ** (-m) * x ** ((1 - R) * abs(m) / 2) * np.mean(
            z ** (-(m + abs(m)) // 2 - n) * chir_ratio(z, m, R, x))
    tot = sum(S.values())
    dual = (w * x ** (-R)) ** (n / 2) * free_chiral(w, n, RT, x)
    ms = np.array(sorted(S)); a = np.array([abs(S[m]) for m in ms])
    top = ms[np.argsort(-a)[:5]]
    fixed = abs(sum(S[m] for m in ms if abs(m) <= 3)) / abs(tot)
    # dual: leading phase predicted by Li2 asymptotics of the Pochhammer ratio
    u = x ** abs(n)
    phase_pred = 2 * float(mp.im(mp.polylog(2, u * np.exp(1j * theta)))) / eps
    print(f"eps={eps:<5} n={n:<3} |I_E/dual-1|={abs(tot/dual-1):.1e}  |I_E|={abs(tot):.4f}  "
          f"share |m|<=3: {fixed:.3f}  top |S_m| at m={list(top)}  "
          f"arg(dual)={np.angle(dual):+.3f}  [2 Im Li2/eps mod 2pi={(phase_pred+np.pi)%(2*np.pi)-np.pi:+.3f}]")
