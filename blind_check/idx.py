"""Index library: S^2 x S^1 superconformal index, conventions as specified.
Fugacities are passed as LOGS (complex) so that half-integer powers are unambiguous.
x = exp(lx)."""
import numpy as np

def poch(a, q, tol=1e-18):
    # (a;q)_inf, a may be an array
    out = np.ones_like(a, dtype=complex)
    term = a.astype(complex)
    while True:
        out = out * (1 - term)
        if np.max(np.abs(term)) < tol:
            break
        term = term * q
    return out

def chir(lzeta, nu, R, lx):
    """Chir[zeta,nu,R] with zeta = exp(lzeta) (array), integer nu."""
    x = np.exp(lx); q = x * x
    s = (-1.0) ** nu
    zeta = np.exp(lzeta)
    pref = np.exp(-abs(nu) / 2 * (lzeta + (R - 1) * lx))
    num = poch(s / zeta * x ** (2 - R + abs(nu)), q)
    den = poch(s * zeta * x ** (R + abs(nu)), q)
    return pref * num / den

def u1_index(chirals, k, lomega, n, lx, sign, M=40, N=512, extra=None):
    """U(1)_k with chirals list of (charge, lfug, flavourflux, R).
    Z_cl = sign(m) * omega^{-m} z^{-k m - n}.
    sign: function m -> +-1.
    extra(m) : optional extra classical factor."""
    phi = 2 * np.pi * np.arange(N) / N
    lz = 1j * phi
    tot = 0j
    for m in range(-M, M + 1):
        f = np.ones(N, dtype=complex) * sign(m) * np.exp(-m * lomega) * np.exp((-k * m - n) * lz)
        for (c, lf, s, R) in chirals:
            f = f * chir(c * lz + lf, c * m + s, R, lx)
        if extra is not None:
            f = f * extra(m)
        tot += np.mean(f)
    return tot
