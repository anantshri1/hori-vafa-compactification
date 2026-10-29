"""CS pair, vacuum by vacuum. Z_S2[CP^1] = Z_+ + Z_-, Z_s = Higgs-residue family at sigma = -s mu (s = +-1):
  Z_s = sum_m e^{-i theta m} sum_k (-1)^k/(k!(k+|m|)!) e^{kappa(r/2+|m|/2+k)} e^{-i s mu kappa} G(-k+2i s mu)/G(1+|m|+k-2i s mu)
Coulomb side: J0 = (H0^(1) + H0^(2))/2  ->  F_j = e^{r kappa/2} int drho e^{-2i mu rho} H0^(j)(L sqrt(sinh^2(rho/2)+sin^2(theta/2))),
L = 4 e^{kappa/4};  (F_1 + F_2)/2 = Z_+ + Z_-  (verified earlier to 1e-16).
Question: is {F_1/2, F_2/2} = {Z_+, Z_-} (vacuum by vacuum), or related by a non-trivial (Stokes-type) matrix?"""
import mpmath as mp
mp.mp.dps = 20

def Zs(s, r, mu, kap, th, MM=60, K=60):
    tot = 0
    for m in range(-MM, MM+1):
        zm = 0
        for k in range(K):
            zm += (-1)**k/(mp.factorial(k)*mp.factorial(k+abs(m))) * mp.e**(kap*(r/2+abs(m)/2+k)) \
                  * mp.e**(-1j*s*mu*kap) * mp.gamma(-k+2j*s*mu)/mp.gamma(1+abs(m)+k-2j*s*mu)
        tot += mp.e**(-1j*th*m)*zm
    return tot

def F(j, r, mu, kap, th):
    L = 4*mp.e**(kap/4); s2 = mp.sin(th/2)**2
    Y = (lambda X: mp.bessely(0, X))
    J = (lambda X: mp.besselj(0, X))
    H = (lambda X: J(X) + 1j*Y(X)) if j == 1 else (lambda X: J(X) - 1j*Y(X))
    f = lambda u: 2*mp.cos(4*mu*mp.asinh(u))*H(L*mp.sqrt(u*u+s2))*2/mp.sqrt(1+u*u)
    return mp.e**(r*kap/2)*mp.quadosc(f, [0, mp.inf], omega=L)

if __name__ == '__main__':
  for r, mu, kap, th in [(0.3, 0.25, -2.0, 0.8), (0.4, 0.6, -1.0, 2.0), (0.2, 0.1, -3.0, 1.3)]:
    zp, zm = Zs(+1, r, mu, kap, th), Zs(-1, r, mu, kap, th)
    f1, f2 = F(1, r, mu, kap, th)/2, F(2, r, mu, kap, th)/2
    print(f"r={r} mu={mu} kappa={kap} theta={th}")
    print(f"   Z_+ = {mp.nstr(zp, 10)}    Z_- = {mp.nstr(zm, 10)}    sum = {mp.nstr(zp+zm, 10)}")
    print(f"   F1/2 = {mp.nstr(f1, 10)}   F2/2 = {mp.nstr(f2, 10)}   sum = {mp.nstr(f1+f2, 10)}")
