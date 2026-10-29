"""
Pair:  A = U(1)_{+1} with [2,0] (Q1, Q2 charge +1, SU(2) fugacity y, R-charge rA)
       B = U(1)_{-1} with [1,1] (Q charge +1, Qt charge -1, axial fugacity a, R-charge rB)
Index in Anant's convention; CS z^{-k m}; grouped z-powers (all integer here):
  A: zero-points z^{-|m|} x^{(1-rA)|m|}, CS z^{-m}       -> z^{-(m+|m|)}  (neutral for m<0: flat half-Coulomb)
  B: zero-points a^{-|m|} x^{(1-rB)|m|} (z cancels), CS z^{+m}  -> z^{+m}  (no neutral monopole: gapped)
Find the exact map at zero background flux:
  wB = y^p,  a = s * wA^qa,  rB fitted,  optional (-1)^m CS signs on each side.
"""
import sys, numpy as np, itertools
sys.path.insert(0, '.')
from index_lib import chir_ratio

x, N, M = 0.3, 128, 50
z = np.exp(2j*np.pi*(np.arange(N)+0.5)/N)

def sectors_A(y, rA):
    S = {}
    for m in range(-M, M+1):
        f = chir_ratio(z*y, m, rA, x) * chir_ratio(z/y, m, rA, x)
        S[m] = x**((1-rA)*abs(m)) * np.mean(z**(-(m+abs(m))) * f)
    return S

def sectors_B(a, rB):
    S = {}
    for m in range(-M, M+1):
        f = chir_ratio(z*a, m, rB, x) * chir_ratio(a/z, -m, rB, x)
        S[m] = a**(-abs(m)) * x**((1-rB)*abs(m)) * np.mean(z**(m) * f)
    return S

def total(S, w, sgn):
    return sum((sgn**m) * w**(-m) * S[m] for m in S)

rA = 0.4
pts = [(np.exp(0.5j), np.exp(0.8j)), (np.exp(1.4j), np.exp(-0.6j)), (np.exp(2.3j), np.exp(1.9j))]   # (y, wA)
SA = [sectors_A(y, rA) for y, _ in pts]
IA = {sA: [total(SA[i], pts[i][1], sA) for i in range(3)] for sA in (1, -1)}

best = []
rB_grid = np.arange(-0.6, 2.61, 0.05)
for s, qa in itertools.product([1, -1], [1, -1, 0.5, -0.5]):
    for rB in rB_grid:
        SB = [sectors_B(s * pts[i][1]**qa, rB) for i in range(3)]
        for p, sA, sB in itertools.product([1, -1, 2, -2], (1, -1), (1, -1)):
            err = max(abs(IA[sA][i] / total(SB[i], pts[i][0]**p, sB) - 1) for i in range(3))
            best.append((err, s, qa, round(rB, 3), p, sA, sB))
best.sort()
print("rA =", rA)
print("top candidates (err, a-sign, a=wA^qa, rB, wB=y^p, CS sign A, CS sign B):")
for b in best[:8]:
    print("  ", f"{b[0]:.2e}", b[1:])
