"""Higgs limit, flux sectors (m1, m2) of U(2)_{(0,-2)} [4,0] with background topological flux n ~ kappa/eps, vs the
Gr(2,4) GLSM sectors on S^2. 3d sector: z_a = (-1)^{m_a} q^{-i sigma_a}, Anant's chirals, W-boson factor with (-1)^{m_ab},
CS z^{(K m)}, FI w^{-sum m} z_a^{-n}. 2d sector (Benini-Cremonesi): (1/2) int d^2 sigma prod_{a<b}((sigma_ab)^2 + m_ab^2/4)
prod_{a,j} G(R/2 - i(sigma_a+mu_j) - m_a/2)/G(1 - R/2 + i(sigma_a+mu_j) - m_a/2) e^{i kappa sum sigma - i theta sum m}.
We print the ratio 3d/2d per sector; the reduction holds if the ratio is the same (sector-independent) as eps -> 0."""
import numpy as np
from scipy.special import loggamma
from idx import qpoch
r = 0.3; R = 1 - r
mu = np.array([0.31, -0.12, 0.47, -0.66]); theta = 0.7
def grids(L=7.0, N=700):
    s = np.linspace(-L, L, N); return list(np.meshgrid(s, s, indexing="ij")) + [s[1]-s[0]]
def sector3d(eps, m, n):
    x = np.exp(-eps/2); q = x*x
    S1, S2, ds = grids()
    val = np.ones_like(S1, dtype=complex); logmon = np.zeros_like(S1, dtype=complex)
    for a, S in enumerate((S1, S2)):
        nu = m[a]
        for mj in mu:
            zeta = (-1)**nu*q**(-1j*(S + mj))
            # Chir[zeta,nu,R] with zeta already containing (-1)^nu: (zeta x^{R-1})^{-|nu|/2} (zeta^{-1}(-1)^nu x^{2-R+|nu|})/(zeta(-1)^nu x^{R+|nu|})
            ze = q**(-1j*(S + mj))           # zeta*(-1)^nu
            val *= qpoch(x**(2 - R + abs(nu))/ze, q)/qpoch(ze*x**(R + abs(nu)), q)
            logmon += -abs(nu)/2*(np.log(ze) + (R - 1)*np.log(x)) - abs(nu)/2*1j*np.pi*nu  # zero point incl. (-1)^{nu} phase
    d = abs(m[0]-m[1])
    z12 = (-1)**(m[0]-m[1])*q**(-1j*(S1-S2))
    val *= (1 - (-1)**(m[0]-m[1])*z12*x**d)*(1 - (-1)**(m[0]-m[1])/z12*x**d)*x**(-d)
    Km = -(m[0] + m[1])
    for a, S in enumerate((S1, S2)):
        za = (-1)**m[a]*q**(-1j*S)
        val *= za**(Km - n)
    val *= np.exp(-1j*theta*(m[0] + m[1]))
    return 0.5*np.sum(val*np.exp(logmon))*(eps*ds/(2*np.pi))**2
def sector2d(m, kappa):
    S1, S2, ds = grids()
    lv = np.zeros_like(S1, dtype=complex)
    for a, S in enumerate((S1, S2)):
        for mj in mu:
            a_ = S + mj
            lv += loggamma(R/2 - 1j*a_ + abs(m[a])/2) - loggamma(1 - R/2 + 1j*a_ + abs(m[a])/2)
    val = np.exp(lv)*((S1-S2)**2 + (m[0]-m[1])**2/4)*np.exp(1j*kappa*(S1 + S2) - 1j*theta*(m[0] + m[1]))
    return 0.5*np.sum(val)*ds**2
if __name__ == "__main__":
    kappa_target = -1.0
    for eps in (0.2, 0.1, 0.05):
        # one-loop running: total charge F=4 per Cartan -> kappa = -eps n - 4 log(1-q)-like shift; choose n even near target
        n = int(round((-kappa_target - 8*np.log(1 - np.exp(-eps)))/eps/2))*2
        kap = -eps*n - 8*np.log(1 - np.exp(-eps))   # 4 chirals x (1-q)^{-2 i sigma} per Cartan
        out = []
        for m in [(0, 0), (1, 0), (1, 1), (2, 0)]:
            out.append(sector3d(eps, m, n)/sector2d(m, kap))
        print(f"eps={eps} n={n} kappa={kap:.3f}: ratio/ratio(0,0):", " ".join(f"{z/out[0]:.4f}" for z in out), f" ratio(0,0)={out[0]:.3e}", flush=True)
