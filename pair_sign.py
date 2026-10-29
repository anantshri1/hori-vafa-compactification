"""Which sign of A's background flux is the Higgs (CP^1) phase?  3d vacua of U(1)_1 [2,0] with FI xi:
xi of one sign -> Higgs vacuum (CP^1); other sign -> topological vacua at sigma = -xi/k_eff (k_eff = 2 side).
The n<0 run sat at m0 = -n/2 > 0 (k_eff = 2 side): topological phase. Try n > 0 and compare fixed-flux sectors
with the Gamma-function GLSM, with BOTH signs of the log-running term."""
import numpy as np, mpmath as mp
from pair_limit import sectors
mp.mp.dps = 15
rA, mu, theta = 0.3, 0.25, 0.8

def glsm(m, kap):
    f = lambda s: mp.e**(1j*s*kap)*mp.fprod(mp.gamma(rA/2-1j*si+abs(m)/2)/mp.gamma(1-rA/2+1j*si+abs(m)/2) for si in (s+mu, s-mu))
    return complex(mp.quad(f, mp.linspace(-80, 80, 81))/(2*mp.pi))

for kap2d in [2.0, -2.0]:
    print(f"\ntarget 2d FI kappa = {kap2d}:  n = (kappa - 4 log(1-q))/eps  (so n > 0)")
    for eps, N in [(0.2, 1024), (0.1, 2048), (0.05, 4096)]:
        x = np.exp(-eps/2); q = x*x; Lq = np.log(1-q)
        n = int(round((kap2d - 4*Lq)/eps))
        y = q**(-1j*mu); wA = np.exp(1j*theta)
        ch = [(lambda z: z*y, 1, rA), (lambda z: z/y, 1, rA)]
        S = sectors(ch, lambda m: -(m+abs(m)) - n, x, int(70/eps), N, lambda m: x**((1-rA)*abs(m)))
        I = sum(wA**(-m)*S[m] for m in S)
        share = abs(sum(wA**(-m)*S[m] for m in S if abs(m) <= 6))/abs(I)
        ab = {m: abs(S[m]) for m in S}; mpk = max(ab, key=ab.get)
        k1 = -eps*n - 4*Lq         # z^{-n} at z=(-1)^m q^{-i sig} gives q^{i sig n} = e^{-i eps sig n}
        k2 = eps*n - 4*Lq
        G1 = {m: (-1)**(n*m)*(1-q)**(2*rA-2)*glsm(m, k1) for m in (0, 1, -1)}
        print(f"  eps={eps:<5} n={n:<5} share|m|<=6={share:.3f} peak m={mpk}   "
              f"S_m/GLSM_m(kap={k1:+.2f}): " + "  ".join(f"m={m}:{S[m]/G1[m]:.3f}" for m in (0, 1, -1)))
