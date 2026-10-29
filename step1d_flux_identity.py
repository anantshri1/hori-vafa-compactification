"""
Step 1d (prep for the second scaling): exact 3d identity with topological background flux n.
Electric: Z_cl = w^{-m} z^{-k m - n} (k = 1/2, integer sign = +1 as fixed in step 1a);
grouped z-power:  z^{-(m+|m|)/2 - n}.
Dual candidates: Chir(+-w^{+-1}, +-n, R_T), principal-branch zero-point.
"""
import sys, numpy as np
sys.path.insert(0, '.')
from index_lib import chir_ratio, free_chiral

R, x, N, M = 0.2, 0.3, 256, 80
RT = (1 - R) / 2
z = np.exp(2j * np.pi * (np.arange(N) + 0.5) / N)
ratio = {m: chir_ratio(z, m, R, x) for m in range(-M, M + 1)}

def I_E(w, n):
    return sum(w ** (-m) * x ** ((1 - R) * abs(m) / 2) *
               np.mean(z ** (-(m + abs(m)) // 2 - n) * ratio[m]) for m in ratio)

ws = np.exp(1j * np.array([0.3, 1.1, 2.5]))
cands = {}
for sw in [1, -1]:
    for pw in [1, -1]:
        for sn in [1, -1]:
            cands[f"Chir({'-' if sw<0 else ''}w^{pw}, {'-' if sn<0 else ''}n)"] = (sw, pw, sn)

for n in [1, 2, -1, -3]:
    res = []
    for name, (sw, pw, sn) in cands.items():
        err = max(abs(I_E(w, n) / free_chiral(sw * w ** pw, sn * n, RT, x) - 1) for w in ws)
        res.append((err, name))
    res.sort()
    print(f"n = {n:>2}: best {res[0][1]:<16} err {res[0][0]:.1e}   (next best {res[1][0]:.1e})")
