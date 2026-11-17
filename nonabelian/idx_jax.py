"""JAX port of idx.theory_index (same conventions, same signature + mshift + chunk). Sector bookkeeping (Fractions,
zero-point exponents) in Python; integrand over the gauge grid on device, vmapped over chunks of flux sectors.
Validated against idx.py by jax_check.py. Abelian only (vector=None), which is all the mirror quiver needs."""
import numpy as np
from fractions import Fraction as Fr
from itertools import product
import jax, jax.numpy as jnp
jax.config.update("jax_enable_x64", True)

def theory_index(n, chirals, K, Kf, FI, flav_vals, x, M, Ngrid=24, signfn=None, bg=None, mshift=None, chunk=256, sectors=None, per_sector=False):
    bg = bg or {}
    q = x*x
    Kq = int(np.ceil(17*np.log(10)/(-np.log(q)))) + 2
    th = 2*np.pi*(np.arange(Ngrid) + 0.3)/Ngrid
    logz = jnp.stack([1j*g.ravel() for g in np.meshgrid(*([th]*n), indexing='ij')])      # (n, G)
    lf = {f: np.log(v) for f, v in flav_vals.items()}
    Q = np.array([c[0] for c in chirals], float)                                             # (C, n)
    Rs = np.array([float(c[2]) for c in chirals])
    lflav = np.array([sum(c*lf[f] for f, c in c_[1].items()) for c_ in chirals], complex)   # (C,)
    logzeta = jnp.asarray(Q) @ logz + lflav[:, None]                                         # (C, G)
    kk = jnp.arange(Kq)

    @jax.jit
    def kernel(nus, ez, logc):
        # nus (B, C) ints as floats; ez (B, n); logc (B,)
        s = jnp.where(jnp.mod(nus, 2) == 0, 1.0, -1.0)                                      # (B, C)
        an = jnp.abs(nus)
        a1 = s[:, :, None]*x**(2 - Rs[None, :, None] + an[:, :, None])*jnp.exp(-logzeta)[None]   # (B, C, G)
        a2 = s[:, :, None]*x**(Rs[None, :, None] + an[:, :, None])*jnp.exp(logzeta)[None]
        qk = q**kk
        lv = jnp.sum(jnp.log1p(-a1[..., None]*qk) - jnp.log1p(-a2[..., None]*qk), axis=(1, 3))      # (B, G)
        lm = ez @ logz + logc[:, None]
        return jnp.mean(jnp.exp(lv + lm), axis=1)

    rng = sectors if sectors is not None else product(range(-M, M+1), repeat=n)
    batch_nu, batch_ez, batch_c = [], [], []
    total = 0j; out = []
    def flush():
        nonlocal total
        if batch_nu:
            v = kernel(jnp.array(batch_nu, float), jnp.array(batch_ez, float), jnp.array(batch_c))
            total += complex(jnp.sum(v))
            if per_sector: out.extend(np.asarray(v).tolist())
            batch_nu.clear(); batch_ez.clear(); batch_c.clear()
    for m in rng:
        m = np.array(m) + (np.array(mshift) if mshift is not None else 0)
        sign = signfn(m) if signfn is not None else 1
        if sign == 0: continue
        ez = [Fr(0)]*n; ef = {f: Fr(0) for f in flav_vals}; ex = Fr(0)
        for a in range(n): ez[a] += sum(Fr(K[a][b])*int(m[b]) for b in range(n))
        for (a, f), lev in Kf.items():
            ef[f] += Fr(lev)*int(m[a]); ez[a] += Fr(lev)*bg.get(f, 0)
        for a, fl in FI.items():
            for f in fl: ef[f] -= int(m[a]); ez[a] -= bg.get(f, 0)
        nus = []
        for Qc, fl, R in chirals:
            nu = int(sum(Qc[a]*m[a] for a in range(n)) + sum(c*bg.get(f, 0) for f, c in fl.items()))
            for a in range(n): ez[a] += Fr(-Qc[a]*abs(nu), 2)
            for f, c in fl.items(): ef[f] += Fr(-c*abs(nu), 2)
            ex += (Fr(R) - 1)*Fr(-abs(nu), 2)
            nus.append(nu)
        for a in range(n): assert ez[a].denominator == 1, ("non-integer gauge power", a, ez[a], m)
        logc = sum(float(ef[f])*lf[f] for f in flav_vals) + float(ex)*np.log(x) + np.log(complex(sign))
        batch_nu.append(nus); batch_ez.append([float(e) for e in ez]); batch_c.append(logc)
        if len(batch_nu) >= chunk: flush()
    flush()
    return out if per_sector else total
