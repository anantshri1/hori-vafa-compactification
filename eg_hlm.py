"""(1) My T*CP^1 JK elliptic genus vs HLM eq. (3.43) (GLSM II for A1).
    Map read off the one-loop factors: mine(z, f, t) = HLM(z_H = z + xi2, xi1 = f, xi2 = -t), up to sign/normalisation.
(2) HLM eq. (3.17), N = 2 Taub-NUT (ALF) elliptic genus at finite radius g^2, unfolded to the plane:
    E_2 = (g^2/2) sum_{a,b=0,1} int_C dudubar/tau2  phi(u + xi1) exp(-(g^2 pi/tau2) |u + (a+b tau)/2|^2),
    phi(u) = th1(u+z) th1(u-z) / (th1(u+xi2) th1(u-xi2)).  Normalisation: g^2 int dudubar/tau2 e^{-g^2 pi|u|^2/tau2} = 1.
    g^2 -> infinity should give HLM (4.17) = (1/2) sum_{a,b} phi(-(a+b tau)/2 + xi1)."""
import numpy as np, mpmath as mp
from eg_tcp1 import eg, th1

tau = 0.31 + 0.9j
def hlm343(zH, x1, x2):
    t1 = lambda w: th1(tau, w)
    return -(t1(-2*x1 + x2 - zH)*t1(2*x1 - x2 - zH))/(t1(-2*x1)*t1(2*x1 - 2*x2)) \
           -(t1(2*x1 + x2 - zH)*t1(-2*x1 - x2 - zH))/(t1(2*x1)*t1(-2*x1 - 2*x2))
for z, f, t in [(0.17+0.05j, 0.23+0.11j, 0.07-0.04j), (0.33-0.02j, 0.12+0.2j, 0.19+0.03j)]:
    mine = eg(tau, z, f, t); theirs = hlm343(z - t, f, -t)
    print(f"(1) mine = {mp.nstr(mine, 12)}   HLM (3.43) = {mp.nstr(theirs, 12)}   ratio = {mp.nstr(mine/theirs, 12)}")

# (2) Taub-NUT N=2 at finite g^2 (numerical 2d integral, polar grid around each Gaussian centre)
z, x1, x2 = 0.21 + 0.03j, 0.13 + 0.07j, 0.29 + 0.05j
t2 = tau.imag
def phi(u):
    return complex(th1(tau, u + x1 + z)*th1(tau, u + x1 - z)/(th1(tau, u + x1 + x2)*th1(tau, u + x1 - x2)))
ale = 0.5*sum(phi(-(a + b*tau)/2) for a in (0, 1) for b in (0, 1))
print(f"\n(2) ALE limit, HLM (4.17) for N=2: {ale:.10f}")
phiv = np.vectorize(phi)
for g2 in [3, 6, 12, 24, 48]:
    w = np.sqrt(t2/(g2*np.pi))                    # Gaussian width
    rs = (np.arange(160) + 0.5)/160 * 6*w         # radial grid up to 6 widths
    ths = 2*np.pi*(np.arange(96) + 0.5)/96
    R, T = np.meshgrid(rs, ths, indexing='ij')
    tot = 0
    for a in (0, 1):
        for b in (0, 1):
            c = -(a + b*tau)/2
            U = c + R*np.exp(1j*T)
            vals = phiv(U) * np.exp(-(g2*np.pi/t2)*R**2)
            # measure: dudubar -> normalised so that g2/t2 * int e^{-g2 pi r^2/t2} (2 pi r dr) * K = 1  =>  K = 1
            tot += np.sum(vals * R) * (rs[1]-rs[0]) * (ths[1]-ths[0])
    E = (g2/2) * tot / t2
    print(f"    g^2 = {g2:>3}: E_TN = {E:.8f}   |E_TN - E_ALE| = {abs(E - ale):.2e}")
