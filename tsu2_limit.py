"""
T[SU(2)] small-circle test.  Exact self-mirror (verified to 1e-15):  I(Y, W, a) = I(W, Y, 1/a).

Scaling:  Y = q^{-i mu}, a = q^{-i nu}  (-> 1),   W = e^{i theta/2} fixed  (w = W^2 = e^{i theta}).
  Copy A = I(Y, W, a):     hypers in Higgs scaling, FI fixed (theta only).
  Copy B = I(W, Y, 1/a):   flavour fugacity W fixed, topological Y^2 -> 1: Coulomb limit.

Copy A prediction, sector by sector (step 0 formula for each chiral at z = (-1)^m q^{-i sigma}):
  S_m(eps) ~ (eps/2pi) * eps^{2r-2-4i nu} * int dsigma prod_i G(r/2 - i s_i + |m|/2)/G(1 - r/2 + i s_i + |m|/2)
             * G((2-2r)/2 + 2i nu)/G(1 - (2-2r)/2 - 2i nu)
  s_i = sigma+mu+nu, sigma-mu+nu, -sigma-mu+nu, -sigma+mu+nu.  (This is the Benini-Cremonesi S^2 partition
  function of the N=(2,2) T*CP^1 GLSM sector m; the |m| vs BC signs cancel pairwise here.)
"""
import sys, numpy as np, mpmath as mp
from tsu2_lib import qinf

r = 0.5
mu, nu, theta = 0.3, 0.15, 0.9


def sectors(Y, a, x, M, N):
    """Flux-sector integrals of T[SU(2)] (without w^{-m}), streaming over |m| to save memory."""
    q = x * x
    z = np.exp(2j * np.pi * (np.arange(N) + 0.5) / N)
    fugs = [z * Y * a, z / Y * a, a / (z * Y), a * Y / z]
    adj = qinf(a ** 2 * x ** (2 - (2 - 2 * r)), q) / qinf(a ** -2 * x ** (2 - 2 * r), q)
    S = {}
    for p in (0, 1):
        s = (-1.0) ** p
        cn = [s / f * x ** (2 - r + p) for f in fugs]
        cd = [s * f * x ** (r + p) for f in fugs]
        Pn = [qinf(c, q) for c in cn]; Pd = [qinf(c, q) for c in cd]
        fn = [np.ones(N, complex) for _ in fugs]; fd = [np.ones(N, complex) for _ in fugs]
        j, mabs = 0, p
        while mabs <= M:
            prod = np.ones(N, complex)
            for i in range(4):
                prod *= (Pn[i] / fn[i]) / (Pd[i] / fd[i])
            val = adj * a ** (-2 * mabs) * x ** (2 * (1 - r) * mabs) * np.mean(prod)
            S[mabs] = val
            if mabs:
                S[-mabs] = val          # integrand depends on |m| only (no CS)
            for i in range(4):
                fn[i] = fn[i] * (1 - cn[i] * q ** j); fd[i] = fd[i] * (1 - cd[i] * q ** j)
            j += 1; mabs += 2
    return S


def Z2d_sector(m, eps):
    """Predicted copy-A sector from the 2d T*CP^1 GLSM (including the eps prefactor)."""
    mp.mp.dps = 20
    def integrand(sig):
        val = mp.mpf(1)
        for si in (sig + mu + nu, sig - mu + nu, -sig - mu + nu, -sig + mu + nu):
            val *= mp.gamma(r/2 - 1j*si + abs(m)/2) / mp.gamma(1 - r/2 + 1j*si + abs(m)/2)
        return val
    I = mp.quad(integrand, [-mp.inf, -5, 0, 5, mp.inf])
    Rphi = 2 - 2*r
    adj = mp.gamma(Rphi/2 + 2j*nu) / mp.gamma(1 - Rphi/2 - 2j*nu)
    return complex(eps/(2*mp.pi) * eps**(2*r - 2 - 4j*nu) * I * adj)


if __name__ == "__main__":
    W = np.exp(1j * theta / 2); w = W * W
    for eps, N in [(0.2, 1024), (0.1, 2048), (0.05, 4096)]:
        x = np.exp(-eps / 2); q = x * x
        Y = q ** (-1j * mu); a = q ** (-1j * nu)
        M = int(70 / eps)
        SA = sectors(Y, a, x, M, N)
        SB = sectors(W, 1 / a, x, M, N)
        IA = sum(w ** (-m) * SA[m] for m in SA)
        IB = sum((Y * Y) ** (-m) * SB[m] for m in SB)
        def diag(S, top, I):
            share = abs(sum(top ** (-m) * S[m] for m in S if abs(m) <= 3)) / abs(I)
            ab = {m: abs(S[m]) for m in S}
            mean = sum(abs(m) * v for m, v in ab.items()) / sum(ab.values())
            return share, mean
        shA, mA = diag(SA, w, IA); shB, mB = diag(SB, Y * Y, IB)
        print(f"\neps = {eps}:  |I_A/I_B - 1| = {abs(IA/IB - 1):.1e}  (self-mirror at small eps)")
        print(f"  copy A (Higgs):   share |m|<=3 = {shA:.3f}   eps<|m|> = {eps*mA:.3f}")
        print(f"  copy B (Coulomb): share |m|<=3 = {shB:.3f}   eps<|m|> = {eps*mB:.3f}")
        for m in [0, 1, 2]:
            Zp = Z2d_sector(m, eps)
            print(f"  copy A sector m={m}:  S_m / 2d GLSM sector - 1 = {abs(SA[m]/Zp - 1):.2e}")
