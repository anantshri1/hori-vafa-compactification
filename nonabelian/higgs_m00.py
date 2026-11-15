"""Higgs limit, zero-flux sector of U(2)_{(0,-2)} with [4,0]: the 3d sector (contour integral over the maximal torus)
vs the 2d U(2) GLSM sector with 4 fundamentals (Gr(2,4)): Vandermonde eps^2 sigma_12^2 and Gamma ratios, with the
(1-q)^{R-1-2i(sigma+mu)} factors of the basic asymptotics kept. Masses y_j = q^{-i mu_j}, sum mu_j = 0."""
import numpy as np
from scipy.special import loggamma
from idx import qpoch
r = 0.3; R = 1 - r
mu = np.array([0.31, -0.12, 0.47, -0.66])
def sector3d(eps, N=1200, L=None):
    x = np.exp(-eps/2); q = x*x
    # z_a = q^{-i sigma_a}: sigma in a window where the integrand lives (|sigma| << 1/eps), uniform grid
    L = L or 8.0
    s = np.linspace(-L, L, N); ds = s[1]-s[0]
    S1, S2 = np.meshgrid(s, s, indexing='ij')
    val = np.ones_like(S1, dtype=complex)
    for S in (S1, S2):
        for m_ in mu:
            zeta = q**(-1j*(S + m_))
            val *= qpoch(x**(2-R)/zeta, q)/qpoch(zeta*x**R, q)
    z12 = q**(-1j*(S1-S2))
    val *= (1 - z12)*(1 - 1/z12)
    # measure: dz/(2 pi i z) = eps dsigma/(2 pi) per Cartan; Weyl 1/2
    return 0.5*np.sum(val)*(eps*ds/(2*np.pi))**2
def sector2d(eps, N=1200, L=8.0):
    s = np.linspace(-L, L, N); ds = s[1]-s[0]
    S1, S2 = np.meshgrid(s, s, indexing='ij')
    lv = np.zeros_like(S1, dtype=complex)
    for S in (S1, S2):
        for m_ in mu:
            a = S + m_
            lv += (R - 1 - 2j*a)*np.log(1 - np.exp(-eps)) + loggamma(R/2 - 1j*a) - loggamma(1 - R/2 + 1j*a)
    val = np.exp(lv) * eps**2*(S1-S2)**2
    return 0.5*np.sum(val)*(eps*ds/(2*np.pi))**2
for eps in (0.2, 0.1, 0.05):
    a, b = sector3d(eps), sector2d(eps)
    print(f"eps={eps}: 3d {a:.6e}  2d {b:.6e}  |ratio-1| {abs(a/b-1):.2e}", flush=True)
