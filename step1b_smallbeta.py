"""
Step 1b: small-circle limit of U(1)_{1/2} + chiral <-> free chiral.

x = e^{-eps/2}, q = e^{-eps}; topological fugacity scaled as w = q^{-i mu} (so the dual
free chiral is in its 2d 'Higgs' scaling, step 0 test C with nu = 0).

Dual side prediction (step 0):  Chir(w, 0, R_T) ~ (1-q)^{R_T-1-2i mu} G(R_T/2 - i mu)/G(1 - R_T/2 + i mu).

Electric side: brute force, decomposed by flux sector, to see WHERE the 2d physics sits.
Sign convention: the one fixed in step 1a (no (-1)^{m/2}; integer sign = +1).
"""
import numpy as np
from math import gamma as _g
from scipy.special import gamma as cgamma
from index_lib import chir_ratio, free_chiral

R, mu = 0.2, 0.35
RT = (1 - R) / 2

def electric_sectors(eps, N, M):
    x = np.exp(-eps / 2)
    z = np.exp(2j * np.pi * (np.arange(N) + 0.5) / N)
    w = np.exp(1j * eps * mu)            # w = q^{-i mu}
    S = {}
    for m in range(-M, M + 1):
        integ = z ** (-(m + abs(m)) // 2) * chir_ratio(z, m, R, x)
        S[m] = w ** (-m) * x ** ((1 - R) * abs(m) / 2) * np.mean(integ)
    return S, x, w

print(f"R = {R}, mu = {mu}, R_T = {RT}")
for eps, N, M in [(0.4, 512, 400), (0.2, 1024, 800), (0.1, 2048, 1600), (0.05, 4096, 3200)]:
    S, x, w = electric_sectors(eps, N, M)
    tot = sum(S.values())
    dual = free_chiral(w, 0, RT, x)
    q = x * x
    pred2d = (1 - q) ** (RT - 1 - 2j * mu) * cgamma(RT / 2 - 1j * mu) / cgamma(1 - RT / 2 + 1j * mu)
    ms = np.array(sorted(S))
    absS = np.array([abs(S[m]) for m in ms])
    # where is the weight?  weighted mean of |m| and share from |m| <= 3
    mean_absm = np.sum(np.abs(ms) * absS) / np.sum(absS)
    small = abs(sum(S[m] for m in ms if abs(m) <= 3)) / abs(tot)
    pos = sum(abs(S[m]) for m in ms if m > 0) / np.sum(absS)
    print(f"\neps = {eps}")
    print(f"  |I_E / I_dual - 1|          = {abs(tot/dual - 1):.2e}   (3d identity, numerics check)")
    print(f"  |I_dual / 2d prediction - 1| = {abs(dual/pred2d - 1):.2e}   (dual side -> 2d free chiral)")
    print(f"  weighted mean |m|           = {mean_absm:8.2f}    eps * <|m|> = {eps*mean_absm:.3f}")
    print(f"  |sum over |m|<=3| / |total| = {small:.3f}")
    print(f"  share of sum|S_m| from m > 0 = {pos:.2e}")
    print(f"  tail |S_-M| / |S_0|          = {abs(S[-M])/abs(S[0]):.1e}")
