"""
T[SU(2)] copy B (Coulomb limit): derived continuum profile and the resulting 2d identity.

Derivation (Li2 terms cancel pairwise between Q and Qt; each chiral with fugacity zeta0 q^{i nu} leaves
(1 - s u/zeta0)^{(r-1)/2 + i nu} (1 - s u zeta0)^{(r-1)/2 + i nu}, u = x^{|m|} = e^{-rho/2}):
   S^B_m  ->  eps^{-4 i nu} * A_B * g(rho),   rho = eps |m|,
   g(rho) = e^{-2 i nu rho} u^{2(1-r)} oint dz/(2 pi i z) [(1-uzW)(1-uz/W)(1-u/(zW))(1-uW/z)]^{r-1+2i nu}
   A_B    = G(R_phi/2 - 2i nu)/G(1 - R_phi/2 + 2i nu)       (adjoint, a_B^{-2} = a^2 = q^{-2i nu})
Then  I_B ~ eps^{-1-4i nu} A_B int_{-inf}^{inf} drho e^{-2 i mu rho} g(|rho|)   and the self-mirror identity
I_A = I_B predicts the 2d statement
   Z_{S^2}[T*CP^1 GLSM](mu, nu, theta) * A_A = A_B * int drho e^{-2i mu rho} g(|rho|)        (*)
with A_A = G(R_phi/2 + 2i nu)/G(1 - R_phi/2 - 2i nu) the copy-A adjoint factor.
"""
import numpy as np, mpmath as mp
from tsu2_limit import sectors, r, mu, nu, theta

Rphi = 2 - 2*r
A_A = complex(mp.gamma(Rphi/2 + 2j*nu) / mp.gamma(1 - Rphi/2 - 2j*nu))
A_B = complex(mp.gamma(Rphi/2 - 2j*nu) / mp.gamma(1 - Rphi/2 + 2j*nu))
W = np.exp(1j*theta/2)
expo = r - 1 + 2j*nu

def g(rho, N=8192):
    u = np.exp(-rho/2)
    z = np.exp(2j*np.pi*(np.arange(N)+0.5)/N)
    prod = (1-u*z*W)*(1-u*z/W)*(1-u/(z*W))*(1-u*W/z)
    return np.exp(+2j*nu*rho) * u**(2*(1-r)) * np.mean(prod**expo)

# (1) sector-level check of the derived profile against the exact copy-B sectors
print("(1) exact copy-B sector vs derived Coulomb-gas profile  eps^{-4i nu} A_B g(rho):")
for eps, N in [(0.2, 1024), (0.1, 2048), (0.05, 4096)]:
    x = np.exp(-eps/2); q = x*x
    a = q**(-1j*nu)
    SB = sectors(W, 1/a, x, int(70/eps), N)
    errs = []
    for rho in [0.4, 1.2, 2.4, 4.0]:
        m = int(round(rho/eps)); rh = eps*m
        errs.append(abs(SB[m] / (eps**(-4j*nu) * A_B * g(rh)) - 1))
    print(f"  eps={eps:<5} rel. errors at rho = 0.4, 1.2, 2.4, 4.0: " + "  ".join(f"{e:.1e}" for e in errs))

# (2) the 2d identity (*): left = 2d GLSM (Benini-Cremonesi), right = Fourier transform of the Coulomb gas
mp.mp.dps = 20
def Zm(m):
    f = lambda sig: mp.fprod(mp.gamma(r/2 - 1j*si + abs(m)/2) / mp.gamma(1 - r/2 + 1j*si + abs(m)/2)
                             for si in (sig+mu+nu, sig-mu+nu, -sig-mu+nu, -sig+mu+nu))
    return complex(mp.quad(f, [-mp.inf, -5, 0, 5, mp.inf]) / (2*mp.pi))
for MM in [40, 160]:
    Zt = Zm(0) + sum(2*np.cos(theta*m)*Zm(m) for m in range(1, MM+1))
    print(f"\n(2) 2d GLSM side, flux cutoff {MM}: Z*A_A = {Zt*A_A:.6f}")
drho = 0.002
rhos = np.arange(drho/2, 60, drho)
gv = np.array([g(rh, N=4096) for rh in rhos])
right = A_B * 2 * np.sum(np.cos(2*mu*rhos) * gv) * drho          # g(|rho|) even in rho
print(f"    Coulomb-gas side: A_B * int e^(-2i mu rho) g(|rho|) = {right:.6f}")
