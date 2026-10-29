"""
Small circle for A = U(1)_1 [2,0]  <->  B = U(1)_{-1} [1,1].
Exact map (verified 1e-13): wB = y^2, a = wA^{1/2}, rB = (1 - rA)/2.

Limit H: y = q^{-i mu} (A's SU(2) mass, B's FI scaled), wA = e^{i theta} fixed.
  A predicted: Higgs -> CP^1 GLSM sectors
     S^A_m ~ (1-q)^{2rA-2} (1/2pi) int dsigma (1-q)^{-4i sigma} prod_pm G(rA/2 - i(sig pm mu) + |m|/2)/G(1 - rA/2 + i(sig pm mu) + |m|/2)
     [(1-q)^{-4i sigma}: total gauge charge 2 -> FI running like log eps (CP^1 asymptotic freedom)]
  B predicted: Coulomb -> flux continuum.
Limit C (reverse): a = q^{-i nu} (B's axial mass, A's FI scaled), y = e^{i phi} fixed.
  B in its Higgs limit would be the 2d U(1) GLSM with charges (+1,-1): expected unbounded.
"""
import numpy as np, mpmath as mp
from tsu2_lib import qinf

def sectors(chirals, zpow, x, M, N, pref):
    """chirals: list of (fug(z), gauge charge e, R). Sector m: prod Chir ratios with flux e*m,
       times z^{zpow(m)} and pref(m); constant term in z."""
    q = x * x
    z = np.exp(2j*np.pi*(np.arange(N)+0.5)/N)
    S = {}
    for p in (0, 1):
        s = (-1.0)**p
        data = []
        for fug, e, R in chirals:
            f = fug(z)
            cn = s / f * x**(2 - R + p); cd = s * f * x**(R + p)
            data.append([cn, cd, qinf(cn, q), qinf(cd, q), np.ones(N, complex), np.ones(N, complex)])
        j, mabs = 0, p
        while mabs <= M:
            prod = np.ones(N, complex)
            for d in data:
                prod *= (d[2]/d[4]) / (d[3]/d[5])
            for m in ([mabs, -mabs] if mabs else [0]):
                S[m] = pref(m) * np.mean(z**zpow(m) * prod)
            for d in data:
                d[4] = d[4]*(1 - d[0]*q**j); d[5] = d[5]*(1 - d[1]*q**j)
            j += 1; mabs += 2
    return S

def theory_A(y, rA, x, M, N):
    ch = [(lambda z: z*y, 1, rA), (lambda z: z/y, 1, rA)]
    return sectors(ch, lambda m: -(m + abs(m)), x, M, N, lambda m: x**((1-rA)*abs(m)))

def theory_B(a, rB, x, M, N):
    ch = [(lambda z: z*a, 1, rB), (lambda z: a/z, -1, rB)]
    return sectors(ch, lambda m: m, x, M, N, lambda m: a**(-abs(m)) * x**((1-rB)*abs(m)))

def tot(S, w):
    return sum(w**(-m) * S[m] for m in S)

def diag(S, w, I):
    share = abs(sum(w**(-m)*S[m] for m in S if abs(m) <= 3)) / abs(I)
    ab = {m: abs(S[m]) for m in S}
    mean = sum(abs(m)*v for m, v in ab.items()) / sum(ab.values())
    neg = sum(v for m, v in ab.items() if m < 0) / sum(ab.values())
    return share, mean, neg

mp.mp.dps = 15
def glsm_sector(m, mu, rA, Lq):
    """(1/2pi) int dsig e^{-4 i sig Lq} prod Gamma ratios;  Lq = log(1-q)."""
    def f(sig):
        v = mp.e**(-4j*sig*Lq)
        for si in (sig + mu, sig - mu):
            v *= mp.gamma(rA/2 - 1j*si + abs(m)/2) / mp.gamma(1 - rA/2 + 1j*si + abs(m)/2)
        return v
    return complex(mp.quad(f, mp.linspace(-60, 60, 41)) / (2*mp.pi))

if __name__ == "__main__":
    rA, mu, theta = 0.3, 0.25, 0.8
    rB = (1 - rA)/2
    print(f"Limit H: y = q^(-i mu) (mu={mu}), wA = e^(i theta) (theta={theta}), rA={rA}, rB={rB}")
    for eps, N in [(0.2, 1024), (0.1, 2048), (0.05, 4096)]:
        x = np.exp(-eps/2); q = x*x
        y = q**(-1j*mu); wA = np.exp(1j*theta)
        M = int(70/eps)
        SA = theory_A(y, rA, x, M, N); SB = theory_B(wA**0.5, rB, x, M, N)
        IA = tot(SA, wA); IB = tot(SB, y**2)
        dA = diag(SA, wA, IA); dB = diag(SB, y**2, IB)
        print(f"\n eps={eps}: |I_A/I_B - 1| = {abs(IA/IB-1):.1e}")
        print(f"   A: fixed-flux share {dA[0]:.3f}  eps<|m|> {eps*dA[1]:.3f}  weight in m<0 {dA[2]:.3f}")
        print(f"   B: fixed-flux share {dB[0]:.3f}  eps<|m|> {eps*dB[1]:.3f}  weight in m<0 {dB[2]:.3f}")
        Lq = np.log(1 - q)
        for m in [0, 1, -1, 2]:
            G = (1 - q)**(2*rA - 2) * glsm_sector(m, mu, rA, Lq)
            print(f"   A sector m={m:>2}: |S_m / GLSM_m - 1| = {abs(SA[m]/G - 1):.2e}")
