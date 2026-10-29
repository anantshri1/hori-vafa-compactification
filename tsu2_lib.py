"""
T[SU(2)] superconformal index in Anant's convention (thesis eq. 1.41 / `Chir`).

N=2 content: U(1) gauge, Q_1, Q_2 (charge +1), Qt_1, Qt_2 (charge -1), neutral adjoint Phi,
W = Qt Phi Q.  Fugacities:
  Q_1: z Y a,  Q_2: z Y^-1 a,  Qt_1: z^-1 Y^-1 a,  Qt_2: z^-1 Y a,  Phi: a^-2
  Y = SU(2)_H Cartan fugacity, a = N=4 axial fugacity, topological: w^{-m} with w = W^2.
R-charges: hypers r (=1/2 for the N=4 R-symmetry H+C), Phi 2-2r.
No CS term. Zero-point factors combined analytically (grouping rule):
  prod over the 4 charged chirals (fug x^{r-1})^{-|m|/2} = a^{-2|m|} x^{2(1-r)|m|}   (integer powers, no branch issue).

q-Pochhammers with x^{|m|} = q^{|m|/2}:  (c q^j; q)_inf = (c;q)_inf / (c;q)_j,  c = s b x^{|m| mod 2},
s = (-1)^m, so each chiral needs only two infinite products (per parity) and a cheap recursion in j.
"""
import numpy as np


def qinf(a, q, tol=1e-18):
    """(a;q)_inf for an array a (|a| <= ~1), q real in (0,1)."""
    a = np.asarray(a, dtype=complex)
    out = np.ones_like(a)
    qk = 1.0
    amax = np.max(np.abs(a))
    while qk * amax > tol:
        out *= 1.0 - a * qk
        qk *= q
    return out


def sector_ratios(b_num, b_den, x, M, r):
    """
    For one chiral with Chir ratio  (s b_num x^{2-r+|m|}; q) / (s b_den x^{r+|m|}; q),  s = (-1)^m,
    return dict m -> ratio array for |m| <= M.  b_num, b_den are arrays over the z grid.
    """
    q = x * x
    out = {}
    for p in (0, 1):
        s = (-1.0) ** p
        cn = s * b_num * x ** (2 - r + p)
        cd = s * b_den * x ** (r + p)
        Pn, Pd = qinf(cn, q), qinf(cd, q)
        fn = np.ones_like(Pn); fd = np.ones_like(Pd)      # (c;q)_j
        j = 0
        mabs = p
        while mabs <= M:
            val = (Pn / fn) / (Pd / fd)
            out[mabs] = val
            fn = fn * (1 - cn * q ** j); fd = fd * (1 - cd * q ** j)
            j += 1
            mabs += 2
    return out


def tsu2_sectors(Y, a, x, M, N, r=0.5):
    """
    Return dict m -> integral over |z|=1 of the T[SU(2)] integrand in flux sector m,
    WITHOUT the topological factor w^{-m} (so the w-dependence can be applied afterwards).
    """
    q = x * x
    z = np.exp(2j * np.pi * (np.arange(N) + 0.5) / N)
    # Chir(zeta, nu, R) ratio: (zeta^-1 s x^{2-R+|nu|}; q) / (zeta s x^{R+|nu|}; q)
    fugs = [z * Y * a, z / Y * a, a / (z * Y), a * Y / z]
    R_ = [r, r, r, r]
    tables = [sector_ratios(1 / f, f, x, M, Rr) for f, Rr in zip(fugs, R_)]
    adj = qinf(a ** 2 * x ** (2 - (2 - 2 * r)), q) / qinf(a ** -2 * x ** (2 - 2 * r), q)
    S = {}
    for m in range(-M, M + 1):
        prod = np.ones(N, dtype=complex)
        for t in tables:
            prod *= t[abs(m)]
        S[m] = adj * a ** (-2 * abs(m)) * x ** (2 * (1 - r) * abs(m)) * np.mean(prod)
    return S


def tsu2_index(Y, W, a, x, M=60, N=256, r=0.5):
    S = tsu2_sectors(Y, a, x, M, N, r)
    w = W * W
    return sum(w ** (-m) * S[m] for m in S), S
