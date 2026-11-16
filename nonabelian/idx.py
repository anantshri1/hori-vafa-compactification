"""General S^2xS^1 index engine for the non-abelian checks (Anant's conventions, as in the abelian work):
  Chir[zeta, nu, R] = (zeta x^{R-1})^{-|nu|/2} (zeta^{-1} (-1)^nu x^{2-R+|nu|}; x^2) / (zeta (-1)^nu x^{R+|nu|}; x^2)
Theory = gauge U(1)^n (abelian quiver) or U(N) electric, given by:
  chirals: list of (Q (gauge charge vector), flav (dict flavour->charge), R)
  K: gauge-gauge CS matrix in the convention e^{-i pi z.K.z} on S^3  ->  index factor prod_a z_a^{(K m)_a}
  Kf: gauge-flavour mixed CS (dict (a, f) -> level, same convention, symmetric pair counted once as e^{-i pi 2 Kf z_a f})
  FI: dict a -> list of flavours (FI fugacity omega_f^{-m_a})
All half-integer powers of each fugacity are collected into one exponent BEFORE evaluation (grouping rule).
Signs: integer CS sign rule = +1 (as validated in the abelian CS pair)."""
import numpy as np
from itertools import product
from fractions import Fraction as Fr

def qpoch(a, q, tol=1e-17):
    a = np.asarray(a, dtype=complex); out = np.ones_like(a); t = a.copy()
    while np.max(np.abs(t)) > tol:
        out *= (1 - t); t *= q
    return out

def theory_index(n, chirals, K, Kf, FI, flav_vals, x, M, Ngrid=24, vector=None, weyl=1, signfn=None, bg=None, mshift=None):
    bg = bg or {}
    """n gauge U(1) factors (or Cartan of U(N) with vector=N). flav_vals: dict flavour -> complex fugacity value
    (log-free: we need log for fractional powers -> pass as (value, log) via exp). Returns the index (complex)."""
    q = x*x
    th = 2*np.pi*(np.arange(Ngrid) + 0.3)/Ngrid
    grids = np.meshgrid(*([th]*n), indexing='ij')
    logz = [1j*g for g in grids]
    lf = {f: np.log(v) for f, v in flav_vals.items()}
    total = 0j
    for m in product(range(-M, M+1), repeat=n):
        m = np.array(m) + (np.array(mshift) if mshift is not None else 0)
        # exponents (as Fractions) of each gauge fugacity and flavour, plus constant x-power
        ez = [Fr(0)]*n; ef = {f: Fr(0) for f in flav_vals}; ex = Fr(0); sign = 1
        # CS
        Km = [sum(Fr(K[a][b])*int(m[b]) for b in range(n)) for a in range(n)]
        for a in range(n): ez[a] += Km[a]
        for (a, f), lev in Kf.items():
            ef[f] += Fr(lev)*int(m[a]); ez[a] += Fr(lev)*bg.get(f, 0)
        for a, fl in FI.items():
            for f in fl: ef[f] -= int(m[a]); ez[a] -= bg.get(f, 0)
        val = np.ones(grids[0].shape, dtype=complex)
        for Q, fl, R in chirals:
            nu = int(sum(Q[a]*m[a] for a in range(n)) + sum(c*bg.get(f, 0) for f, c in fl.items()))
            R = Fr(R)
            # zero point (zeta x^{R-1})^{-|nu|/2}
            for a in range(n): ez[a] += Fr(-Q[a]*abs(nu), 2)
            for f, c in fl.items(): ef[f] += Fr(-c*abs(nu), 2)
            ex += (R - 1)*Fr(-abs(nu), 2)
            logzeta = sum(Q[a]*logz[a] for a in range(n)) + sum(c*lf[f] for f, c in fl.items())
            zeta = np.exp(logzeta)
            s = (-1)**nu
            val *= qpoch(s*x**(2 - float(R) + abs(nu))/zeta, q)/qpoch(s*zeta*x**(float(R) + abs(nu)), q)
        if vector:
            N = vector
            for a in range(N):
                for b in range(N):
                    if a == b: continue
                    d = abs(int(m[a] - m[b]))
                    val *= (1 - (-1)**int(m[a] - m[b])*np.exp(logz[a] - logz[b])*x**d)   # (-1)^nu convention of Chir carried to W-bosons
                    ex += Fr(-d, 2)
        for a in range(n):
            assert ez[a].denominator == 1, ("non-integer gauge power", a, ez[a], m)
        mon = np.exp(sum(float(ez[a])*logz[a] for a in range(n)) + sum(float(ef[f])*lf[f] for f in flav_vals)) * x**float(ex)
        if signfn is not None: sign = signfn(m)
        total += sign*np.mean(val*mon)
    return total/weyl
