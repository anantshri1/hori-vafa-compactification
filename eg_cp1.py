"""Elliptic genus of the CP^1 GLSM (Higgs side of A = U(1)_1[2,0]): U(1), two chirals of charge +1, R = 0,
SU(2) flavour holonomy +-f. Integrand picks up e^{2 pi i z sum Q} = e^{4 pi i z} under u -> u + tau, so it is elliptic
(and the JK answer eta-independent) only for z in Z/2, i.e. y = +-1. Test: eta>0 vs eta<0 at generic z and at z = 1/2."""
import mpmath as mp
from eg_tcp1 import th1, chiral
mp.mp.dps = 20
tau = 0.31 + 0.9j

def cp1(z, f, side):
    # charge +1 chirals: a = u + f, u - f  (R = 0);  no negatively charged fields -> eta<0 picks poles "at infinity":
    # instead compute eta<0 as minus the sum over all poles in a fundamental cell of the integrand (residue theorem on the torus)
    poles = [-f, f]
    if side > 0:
        return sum(chiral(tau, z, 0, u0 + (f if u0 == -f else -f)) for u0 in poles)
    return None

def total_residue(z, f):
    """Sum of residues of Z(u) du over a fundamental cell = 0 iff Z is elliptic in u. Check it numerically by a contour
    integral around the boundary of a cell (parallelogram), which must vanish for an elliptic integrand."""
    vec = 2*mp.pi*mp.eta(tau)**3 if hasattr(mp, 'eta') else None
    def Z(u):
        return chiral(tau, z, 0, u + f) * chiral(tau, z, 0, u - f) * (th1(tau, -z)/1)  # vector factor up to const
    c = 0.05 + 0.03j     # cell corner, chosen to avoid poles on the boundary
    edges = [(c, c + 1), (c + 1, c + 1 + tau), (c + 1 + tau, c + tau), (c + tau, c)]
    I = 0
    for a, b in edges:
        I += mp.quad(lambda s: Z(a + s*(b - a))*(b - a), [0, 1])
    return I

f = 0.23 + 0.11j
for z in [0.17 + 0.05j, 0.31, 0.5, 0.5 + 0.0j]:
    print(f"z = {z}:  boundary integral of the integrand over a cell = {mp.nstr(total_residue(z, f), 6)}   "
          f"(must vanish for a well-defined, eta-independent elliptic genus)")
print("\nCP^1 elliptic genus at y = -1 (z = 1/2), JK eta>0, as a function of tau (q-series is non-trivial):")
for t in [0.31 + 0.9j, 0.31 + 1.4j, 0.31 + 2.5j]:
    tau = t
    E = sum(chiral(tau, 0.5, 0, u0 + s) for u0, s in [(-f, -f), (f, f)])
    print(f"  tau = {t}:  E(y=-1, f) = {mp.nstr(E, 10)}")
