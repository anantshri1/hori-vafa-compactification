"""
Numerical 3d N=2 superconformal index in Anant's convention (thesis eq. 1.41 /
3dSCI notebook `Chir`), double precision, brute force:

  * flux sum truncated at |m| <= M
  * contour integral over |z| = 1 by the trapezoid rule with N points
    (exponentially accurate for integrands analytic in an annulus)
  * q-Pochhammers as finite products with a tail cutoff

Grouping rule (step 0, finding 2): all half-integer powers of a gauge fugacity
(chiral zero-points + CS term) are combined into ONE integer power before
evaluation, so the integrand is single-valued on the unit circle.
"""
import numpy as np


def qpoch(a, q, tol=1e-17, kmax=200000):
    """(a; q)_inf = prod_{k>=0} (1 - a q^k), vectorised over a, q real in (0,1)."""
    a = np.asarray(a, dtype=complex)
    out = np.ones_like(a)
    qk = 1.0
    for _ in range(kmax):
        term = a * qk
        out *= 1.0 - term
        qk *= q
        if qk * np.max(np.abs(a)) < tol:
            break
    return out


def chir_ratio(z, m, R, x):
    """Pochhammer ratio of Chir[z, m, R] (zero-point factor NOT included)."""
    s = (-1.0) ** m
    q = x * x
    return qpoch(s * x ** (2 - R + abs(m)) / z, q) / qpoch(s * z * x ** (R + abs(m)), q)


def free_chiral(w, n, R, x):
    """Full Chir[w, n, R] for a free (ungauged) chiral, |w| = 1, principal branch."""
    return (w * x ** (R - 1)) ** (-abs(n) / 2) * chir_ratio(np.asarray(w, complex), n, R, x)
