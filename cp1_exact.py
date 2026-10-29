"""Exponentially convergent evaluation of both sides of the CP^1 identity.
GLSM: close the sigma contour in the lower half plane (kappa < 0); residues at the poles of the two numerators:
  Z_m = sum_{s=+-} sum_k (-1)^k / (k! (k+|m|)!) e^{kappa(r/2+|m|/2+k)} e^{-i s mu kappa} G(-k+2i s mu)/G(1+|m|+k-2i s mu)
J0 side: u = sinh(rho/2) makes the oscillation uniform (period 2pi/L); oscillatory quadrature with quadosc."""
import mpmath as mp
mp.mp.dps = 25

def glsm(r, mu, kap, th, MM=80, K=80):
    def Zm(m):
        tot = 0
        for s in (1, -1):
            for k in range(K):
                tot += (-1)**k/(mp.factorial(k)*mp.factorial(k+abs(m))) * mp.e**(kap*(r/2+abs(m)/2+k)) \
                       * mp.e**(-1j*s*mu*kap) * mp.gamma(-k+2j*s*mu)/mp.gamma(1+abs(m)+k-2j*s*mu)
        return tot
    return Zm(0) + mp.fsum(2*mp.cos(th*m)*Zm(m) for m in range(1, MM+1))

def j0side(r, mu, kap, th):
    L = 4*mp.e**(kap/4); s2 = mp.sin(th/2)**2
    f = lambda u: 2*mp.cos(4*mu*mp.asinh(u))*mp.besselj(0, L*mp.sqrt(u*u+s2))*2/mp.sqrt(1+u*u)
    return mp.e**(r*kap/2)*mp.quadosc(f, [0, mp.inf], omega=L)

if __name__ == "__main__":
 for r, mu, kap, th in [(0.3, 0.25, -2.0, 0.8), (0.4, 0.6, -1.0, 2.0), (0.2, 0.1, -3.0, 1.3), (0.45, 0.8, -0.5, 2.7)]:
    G = glsm(r, mu, kap, th); J = j0side(r, mu, kap, th)
    print(f"r={r} mu={mu} kappa={kap} theta={th}:  GLSM(residues) = {mp.nstr(G, 15)}   J0 side = {mp.nstr(J, 15)}   rel diff = {mp.nstr(abs(G/J-1), 3)}")
