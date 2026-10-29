"""Anant's gauge for B: Q fugacity z (no eta), Qt fugacity wA/z (carries eta), R-charges (cQ, 1-r-cQ),
B topological fugacity y^{+-2} (X2 - X1).  Derived from my symmetric gauge by z -> z wA^{1/2} x^{c}:
the CS term z^{m} then produces wA^{m/2} x^{c m}; grouped with the zero-points the wA powers are integer.
Scan the leftover x^{p m} (gauge/topological R-mixing) and the sign of the y-power."""
import sys, numpy as np
sys.path.insert(0, '.')
from index_lib import chir_ratio
x, N, M = 0.3, 256, 50
z = np.exp(2j*np.pi*(np.arange(N)+0.5)/N)
r = 0.3

def IA(y, wA):
    return sum(wA**(-m)*x**((1-r)*abs(m))*np.mean(z**(-(m+abs(m)))*chir_ratio(z*y, m, r, x)*chir_ratio(z/y, m, r, x))
               for m in range(-M, M+1))

def IB_his(y, wA, cQ, p, sy):
    RT = 1 - r - cQ
    tot = 0
    for m in range(-M, M+1):
        zero = wA**(-(m + abs(m))//2) * x**((2 - cQ - RT)*abs(m)/2)      # grouped zero-points + CS-induced wA^{m/2}
        tot += (y**(2*sy))**(-m) * x**(p*m) * zero * np.mean(z**m * chir_ratio(z, m, cQ, x) * chir_ratio(wA/z, -m, RT, x))
    return tot

pts = [(np.exp(0.5j), np.exp(0.8j)), (np.exp(1.4j), np.exp(-0.6j)), (np.exp(2.3j), np.exp(1.9j))]
cQ = 0.1
res = []
for sy in (1, -1):
    for p in np.arange(-1.0, 1.001, 0.05):
        e = max(abs(IA(y, w)/IB_his(y, w, cQ, p, sy) - 1) for y, w in pts)
        res.append((e, sy, round(p, 3)))
res.sort()
print(f"r = {r}, R(Q) = {cQ}, R(Qt) = {1-r-cQ}: best (err, y-power sign, x^(p m)):", res[:3])
