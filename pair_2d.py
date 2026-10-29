"""Correct phase (nA > 0). Derived B-side asymptotics: chirals at flux N+-m, u ~ eps, (c;q) ~ exp(-c/(1-q)):
   S^B_m -> pref * J0( 4 t sqrt( sinh^2(rho/2) + sin^2(theta/2) ) ),  rho = eps m,  t = x^N/(1-q),
   pref = a^{-N} x^{(1-rB) N},  N = nA/2.
=> I_B -> pref (1/eps) int drho e^{-2 i mu rho} J0(4 t sqrt(sinh^2(rho/2)+sin^2(theta/2))).
A-side: I_A -> eps^{2rA-1} sum_m e^{-i theta m} (-1)^{n m} Z_m^{CP1}(kappa),  kappa = -eps nA - 4 log(1-q).
Exact identity: I_A = (wA x^{3rA-1})^{nA/4} I_B.  So the 2d prediction is
   eps^{2rA-1} Z_CP1 = contact * pref / eps * FourierJ0 .
"""
import numpy as np, mpmath as mp
from scipy.special import j0
from pair_limit import sectors
from pair_Bside import B_sectors
mp.mp.dps = 15
rA, mu, theta, kt = 0.3, 0.25, 0.8, 2.0
rB = (1 - rA)/2

def ZCP1(kap, MM=30):
    def Zm(m):
        f = lambda s: mp.e**(1j*s*kap)*mp.fprod(mp.gamma(rA/2-1j*si+abs(m)/2)/mp.gamma(1-rA/2+1j*si+abs(m)/2) for si in (s+mu, s-mu))
        return complex(mp.quad(f, mp.linspace(-80, 80, 81))/(2*mp.pi))
    return {m: Zm(m) for m in range(-MM, MM+1)}

for eps, N in [(0.2, 1024), (0.1, 2048), (0.05, 4096)]:
    x = np.exp(-eps/2); q = x*x; Lq = np.log(1-q)
    nA = 2*int(round((kt - 4*Lq)/(2*eps))); Nn = nA//2
    kap = -eps*nA - 4*Lq
    y = q**(-1j*mu); wA = np.exp(1j*theta); a = np.exp(1j*theta/2)
    M = int(70/eps)
    # exact A and B
    ch = [(lambda z: z*y, 1, rA), (lambda z: z/y, 1, rA)]
    SA = sectors(ch, lambda m: -(m+abs(m)) - nA, x, M, N, lambda m: x**((1-rA)*abs(m)))
    IA = sum(wA**(-m)*SA[m] for m in SA)
    SB = B_sectors(a, Nn, x, min(M, Nn - 1), N)
    IB = sum((y*y)**(-m)*SB[m] for m in SB)
    contact = (wA*x**(3*rA - 1))**(nA/4)
    t = x**Nn/(1-q); pref = a**(-Nn)*x**((1-rB)*Nn)
    # (1) B sector profile
    errs = []
    for m in [0, int(0.5/eps), int(1.5/eps)]:
        rho = eps*m
        pred = pref*j0(4*t*np.sqrt(np.sinh(rho/2)**2 + np.sin(theta/2)**2))
        errs.append((rho, abs(SB[m]/pred - 1)))
    # (2) B total vs Fourier integral
    rhos = np.linspace(-40, 40, 160001); dr = rhos[1]-rhos[0]
    F = np.sum(np.exp(-2j*mu*rhos)*j0(4*t*np.sqrt(np.sinh(rhos/2)**2 + np.sin(theta/2)**2)))*dr
    IB_pred = pref*F/eps
    # (3) A total vs CP1 GLSM
    Z = ZCP1(kap)
    IA_pred = eps**(2*rA-1)*sum(wA**(-m)*(-1)**(nA*m)*Z[m] for m in Z)
    print(f"eps={eps:<5} nA={nA:<4} kap={kap:+.3f} t={t:.3f}  identity {abs(IA/(contact*IB)-1):.0e}")
    print("   B sector vs J0 profile: " + "  ".join(f"rho={r:.2f}:{e:.1e}" for r, e in errs))
    print(f"   |I_B / Fourier-J0 - 1| = {abs(IB/IB_pred-1):.2e}    |I_A / (eps^(2rA-1) Z_CP1) - 1| = {abs(IA/IA_pred-1):.2e}")
    print(f"   2d identity check: |eps^(2rA-1) Z_CP1 / (contact pref FourierJ0 / eps) - 1| = {abs(IA_pred/(contact*IB_pred)-1):.2e}")
