"""
Elliptic genus of the 2d N=(2,2) T*CP^1 GLSM (the Higgs-limit side of T[SU(2)]), BEHT rank-1 JK residues.
Content: U(1); Q1: u+f+t, Q2: u-f+t (charge +1); Qt1: -u-f+t, Qt2: -u+f+t (charge -1); Phi: -2t.
R-charges: hypers 0, Phi 2 (W = Qt Phi Q). f = SU(2)_H, t = U(1)_t (N=4 axial-type), z = left R.
One-loop: vector 2pi eta^3/theta1(-z); chiral (R, a) -> theta1((R/2-1) z + a)/theta1((R/2) z + a).
At a pole of a charge +1 chiral, vector x pole-chiral gives 1/(u-u0): residue = product of the other factors.
JK with eta > 0 sums the poles of positively charged fields; eta < 0 the negative ones (with a sign).
Checks: (i) eta-independence (needs sum of gauge charges = 0, true here); (ii) z -> 0 limit."""
import mpmath as mp
mp.mp.dps = 20

def th1(tau, z):
    return mp.jtheta(1, mp.pi*z, mp.e**(1j*mp.pi*tau))

def chiral(tau, z, R, a):
    return th1(tau, (R/2 - 1)*z + a) / th1(tau, (R/2)*z + a)

def fields(u, f, t):
    return [(+1, 0, u + f + t), (+1, 0, u - f + t), (-1, 0, -u - f + t), (-1, 0, -u + f + t), (0, 2, -2*t)]

def eg(tau, z, f, t, side=+1):
    # poles: for each charge Q = side chiral, solve (R/2) z + a(u) = 0 -> u0
    tot = 0
    flist = fields(0, f, t)
    for i, (Q, R, a0) in enumerate(flist):
        if Q != side:
            continue
        u0 = -(a0 + (R/2)*z)/Q
        prod = 1
        for j, (Qj, Rj, aj) in enumerate(fields(u0, f, t)):
            if j == i:
                continue
            prod *= chiral(tau, z, Rj, aj)
        tot += prod if side == +1 else prod      # orientation sign: Q=-1 pole gives -1, eta<0 gives -1
    return tot

tau = 0.31 + 0.9j
for z, f, t in [(0.17 + 0.05j, 0.23 + 0.11j, 0.07 - 0.04j), (0.33 - 0.02j, 0.12 + 0.2j, 0.19 + 0.03j)]:
    Ep, Em = eg(tau, z, f, t, +1), eg(tau, z, f, t, -1)
    print(f"z={z}, f={f}, t={t}:  eta>0: {mp.nstr(Ep, 12)}   eta<0: {mp.nstr(Em, 12)}   |diff| = {mp.nstr(abs(Ep-Em), 3)}")
# z -> 0
for z in [1e-3, 1e-5]:
    print(f"z={z}: EG = {mp.nstr(eg(tau, z, 0.23+0.11j, 0.07-0.04j, +1), 10)}  (Euler characteristic of T*CP^1 is 2)")
# double periodicity in u is what makes the JK answer eta-independent; check flavour periodicity f -> f+1, f+tau
z, f, t = 0.17+0.05j, 0.23+0.11j, 0.07-0.04j
E0 = eg(tau, z, f, t)
print("f -> f+1:", mp.nstr(abs(eg(tau, z, f+1, t)/E0 - 1), 3), "  f -> f+tau:", mp.nstr(abs(eg(tau, z, f+tau, t)/E0 - 1), 3))
