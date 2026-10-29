"""
Operator-level check of A = U(1)_1[2,0] <-> B = U(1)_{-1}[1,1].
Resolve the index by charges: coefficient of y^j wA^k (A side; B side via wB = y^2, a = wA^{1/2}),
leading power of x and its coefficient, at small x. Then, for the SU(2) triplet, find WHICH flux
sector on each side produces each component.
"""
import numpy as np
from pair_limit import theory_A, theory_B

rA = 0.3; rB = (1 - rA) / 2
NY, NW = 12, 12
ys = np.exp(2j*np.pi*(np.arange(NY)+0.5)/NY)
ws = np.exp(2j*np.pi*(np.arange(NW)+0.5)/NW)

def resolved(xv, side, M=12, N=128):
    """C[j,k](x): Fourier coefficients in y (j) and wA (k), returned per flux sector m too."""
    C = {}
    for iy, y in enumerate(ys):
        if side == "A":
            S = theory_A(y, rA, xv, M, N)
        for iw, w in enumerate(ws):
            if side == "B":
                S = theory_B(w**0.5 if False else None, rB, xv, M, N) if False else None
            pass
    return C

# simpler, explicit version: per-sector Fourier coefficients
def coeffA(xv, M=10, N=128):
    out = {}
    for iy, y in enumerate(ys):
        S = theory_A(y, rA, xv, M, N)
        for m, v in S.items():
            # topological charge of sector m: w^{-m}  -> k = -m
            for j in range(-4, 5):
                out[(j, -m, m)] = out.get((j, -m, m), 0) + v * y**(-j) / NY
    return out           # key (y-power j, wA-power k, flux m)

def coeffB(xv, M=10, N=128):
    """B with a = wA^{1/2}: need integer wA powers -> sample a on the circle, wA = a^2.
       wB = y^2 means B's flux m carries y^{-2m}."""
    out = {}
    NA = 16
    as_ = np.exp(2j*np.pi*(np.arange(NA)+0.5)/NA)
    for a in as_:
        S = theory_B(a, rB, xv, M, N)
        for m, v in S.items():
            for t in range(-8, 9):          # a-power t  -> wA power t/2
                out[(-2*m, t, m)] = out.get((-2*m, t, m), 0) + v * a**(-t) / NA
    return out           # key (y-power, a-power, flux m)

x1, x2 = 0.02, 0.04
A1, A2 = coeffA(x1), coeffA(x2)
B1, B2 = coeffB(x1), coeffB(x2)

def leading(d1, d2, key):
    c1, c2 = d1.get(key, 0), d2.get(key, 0)
    if abs(c1) < 1e-14 or abs(c2) < 1e-14:
        return None
    p = np.log(abs(c2/c1)) / np.log(x2/x1)
    return p, c1 / x1**p

print(f"rA = {rA}, rB = {rB};  2 rA = {2*rA}")
print("\nA side: SU(2) triplet candidates (y^{+-2}, y^0) and the neutral monopole, by flux sector m")
for j in (2, 0, -2):
    for m in range(-3, 4):
        r = leading(A1, A2, (j, -m, m))
        if r and r[0] < 1.5:
            print(f"  y^{j:+d} wA^{-m:+d}  (flux m={m:+d}):  x^{r[0]:.3f}  coeff {r[1].real:+.3f}")
print("\nB side (wA^k = a^{2k}): same charges, by B flux sector m")
for j in (2, 0, -2):
    for m in range(-3, 4):
        for t in range(-8, 9, 2):
            r = leading(B1, B2, (j, t, m))
            if r and r[0] < 1.5 and (-2*m) == j:
                print(f"  y^{j:+d} wA^{t//2:+d}  (B flux m={m:+d}, a^{t:+d}):  x^{r[0]:.3f}  coeff {r[1].real:+.3f}")
