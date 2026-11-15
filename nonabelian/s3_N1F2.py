"""Calibration of the thesis formulas (Unitary_Mirrors_1.tex, eqs u_n_chiral and mirpartfun) at N=1, F=2:
SQCD side U(1) with [2,0] at the induced CS level  vs  planar mirror (single node). This is the abelian CS pair already checked
in the index, so it tests my reading of the thesis conventions before N=2, F=4."""
import mpmath as mp, sys
from sb import sb, Q
mp.mp.dps = 15
iQ2 = 1j*Q/2
def Z_el(X, eta):
    N, F = 1, 2
    pre = mp.e**(-N*mp.pi*1j/2*sum(x**2 for x in X))
    f = lambda u: mp.e**(2j*mp.pi*eta*u + 1j*mp.pi*u**2 + 1j*mp.pi*(F/2-N)*u**2) * sb(X[0]-u)*sb(X[1]-u)
    return pre*mp.quad(f, [-mp.inf, -3, 0, 3, mp.inf])
def Z_mir(X, eta):
    N = 1
    pre = mp.e**(2j*mp.pi*eta*(-(2*N-1)*1j*Q/4 + X[0])) * mp.e**(-1j*mp.pi/2*eta**2)
    f = lambda z: mp.e**(2j*mp.pi*z*(X[1]-X[0]) - 1j*mp.pi*z**2) * sb(iQ2 - eta + z)*sb(iQ2 - z)*mp.e**(1j*mp.pi*z*eta)
    return pre*mp.quad(f, [-mp.inf, -3, 0, 3, mp.inf])
for X, eta in [((0.15+0.3j, -0.15+0.3j), 0.2-0.25j), ((0.3+0.25j, -0.3+0.25j), -0.1-0.2j)]:
    a, c = Z_el(X, eta), Z_mir(X, eta)
    print("X =", X, "eta =", eta, " Z_el =", mp.nstr(a, 10), " Z_mir =", mp.nstr(c, 10), " ratio =", mp.nstr(a/c, 10), flush=True)
