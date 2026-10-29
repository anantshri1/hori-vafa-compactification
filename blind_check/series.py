"""Exact truncated series engine for the index.
Monomial key: tuple of Fractions (ex, er, ez, f1, f2, ...) meaning x^ex * x^(r*er) * z^ez * prod f_i^{e_i}.
Weight of a monomial = ex + er*r0 (r0 numeric, only used for truncation)."""
from fractions import Fraction as F
from collections import defaultdict

class Ser:
    def __init__(self, nf, r0, L):
        self.nf = nf; self.r0 = r0; self.L = L
    def w(self, k):
        return float(k[0]) + float(k[1]) * self.r0
    def one(self):
        return {tuple([F(0)] * (3 + self.nf)): 1}
    def mono(self, coeff, ex=0, er=0, ez=0, f=None):
        f = f or [0] * self.nf
        return {tuple([F(ex), F(er), F(ez)] + [F(v) for v in f]): coeff}
    def mul(self, A, B):
        out = defaultdict(int)
        for ka, va in A.items():
            wa = self.w(ka)
            if wa > self.L + 1e-9: continue
            for kb, vb in B.items():
                k = tuple(a + b for a, b in zip(ka, kb))
                if self.w(k) > self.L + 1e-9: continue
                out[k] += va * vb
        return {k: v for k, v in out.items() if v != 0}
    def chir(self, zc, f, nu, R):
        """chiral with zeta = z^zc * prod f_i^{f_i}, flux nu, R=(R0,R1) meaning R0+R1*r."""
        R0, R1 = F(R[0]), F(R[1]); a = abs(nu); s = 1 if nu % 2 == 0 else -1
        # prefactor (zeta x^{R-1})^{-a/2}
        pre = self.mono(1, ex=-F(a, 2) * (R0 - 1), er=-F(a, 2) * R1, ez=-F(a, 2) * zc,
                        f=[-F(a, 2) * v for v in f])
        res = pre
        # numerator prod_j (1 - s zeta^{-1} x^{2-R+a+2j})
        j = 0
        while True:
            ex = 2 - R0 + a + 2 * j; er = -R1
            if float(ex) + float(er) * self.r0 > self.L + 1e-9: break
            fac = self.one()
            fac.update(self.mono(-s, ex=ex, er=er, ez=-zc, f=[-v for v in f]))
            res = self.mul(res, fac); j += 1
        # denominator 1/prod (1 - s zeta x^{R+a+2j}) = prod sum_k (s zeta x^..)^k
        j = 0
        while True:
            ex = R0 + a + 2 * j; er = R1
            wt = float(ex) + float(er) * self.r0
            if wt > self.L + 1e-9: break
            assert wt > 0
            fac = self.one(); kk = 1
            while kk * wt <= self.L + 1e-9:
                fac.update(self.mono((s ** kk), ex=kk * ex, er=kk * er, ez=kk * zc, f=[kk * v for v in f]))
                kk += 1
            res = self.mul(res, fac); j += 1
        return res

def u1(S, chirals, cl, M):
    """chirals: list of (zcharge, fvec, fluxfun(m), R); cl(m) -> dict monomial (classical piece incl z^..)."""
    tot = defaultdict(int)
    for m in range(-M, M + 1):
        t = cl(m)
        for (zc, f, fl, R) in chirals:
            t = S.mul(t, S.chir(zc, f, fl(m), R))
        for k, v in t.items():
            if k[2] == 0:
                tot[k] += v
    return {k: v for k, v in tot.items() if v != 0}

def show(S, ser, names, maxw=None):
    items = sorted(ser.items(), key=lambda kv: S.w(kv[0]))
    for k, v in items:
        if maxw is not None and S.w(k) > maxw: continue
        s = f"x^({k[0]}+{k[1]}r)"
        for nm, e in zip(names, k[3:]):
            if e != 0: s += f" {nm}^{e}"
        print(f"{v:+d}  {s}   [w={S.w(k):.3f}]")
