"""
T[SU(2)] test, round 2.
 (1) copy A vs 2d T*CP^1 GLSM with the CORRECTED prefactor eps^{2r-2-4i nu} (the measure eps is already
     inside it; round 1 multiplied by eps twice -> uniform 1/eps mismatch).  Sector by sector AND total.
 (2) copy B: is the flux sum a continuum?  Scaling collapse of S^B at rho = eps*m across eps.
"""
import numpy as np, mpmath as mp
from tsu2_limit import sectors, r, mu, nu, theta

mp.mp.dps = 20
Rphi = 2 - 2 * r
adj2d = complex(mp.gamma(Rphi/2 + 2j*nu) / mp.gamma(1 - Rphi/2 - 2j*nu))

def Zm(m):
    f = lambda sig: mp.fprod(mp.gamma(r/2 - 1j*si + abs(m)/2) / mp.gamma(1 - r/2 + 1j*si + abs(m)/2)
                             for si in (sig+mu+nu, sig-mu+nu, -sig-mu+nu, -sig+mu+nu))
    return complex(mp.quad(f, [-mp.inf, -5, 0, 5, mp.inf]) / (2*mp.pi)) * adj2d

MM = 40
Z = {m: Zm(m) for m in range(0, MM + 1)}
w = np.exp(1j * theta)
Z2d_tot = Z[0] + sum(2 * np.cos(theta * m) * Z[m] for m in range(1, MM + 1))   # Z_m even in m
print(f"2d T*CP^1 GLSM: Z = {Z2d_tot:.6f};  |Z_40|/|Z_0| = {abs(Z[40]/Z[0]):.1e} (m-sum tail)")

W = np.exp(1j * theta / 2)
prof = {}
for eps, N in [(0.2, 1024), (0.1, 2048), (0.05, 4096)]:
    x = np.exp(-eps/2); q = x*x
    Y = q**(-1j*mu); a = q**(-1j*nu)
    M = int(70/eps)
    SA = sectors(Y, a, x, M, N)
    SB = sectors(W, 1/a, x, M, N)
    IA = sum(w**(-m) * SA[m] for m in SA)
    pref = eps**(2*r - 2 - 4j*nu)
    print(f"\neps={eps}:")
    for m in [0, 1, 2, 5]:
        print(f"  copy A sector m={m}:  |S_m/(pref*Z_m) - 1| = {abs(SA[m]/(pref*Z[m]) - 1):.2e}")
    print(f"  copy A total:   |I/(pref*Z_2d) - 1|      = {abs(IA/(pref*Z2d_tot) - 1):.2e}")
    prof[eps] = (SB, pref)

# (2) scaling collapse of copy B's sectors: is eps^{-p} S^B(rho/eps) eps-independent?
print("\ncopy B profile |S^B_m| at fixed rho = eps*m (should collapse up to an eps power if flux is a continuum):")
rhos = [0.4, 0.8, 1.2, 1.6, 2.4, 3.2]
for eps, (SB, pref) in prof.items():
    vals = [abs(SB[int(round(rh/eps))]) for rh in rhos]
    norm = vals[2]
    print(f"  eps={eps:<5}  |S^B| / |S^B(rho=1.2)| at rho={rhos}: " + "  ".join(f"{v/norm:.4f}" for v in vals)
          + f"   |S^B(1.2)| = {norm:.3e}")
