"""CS pair, vacuum by vacuum, done via factorisation (blocks) instead of the Hankel split.
Claim A (Higgs/GLSM side): each Higgs-residue family factorises exactly,
   Z_s = e^{r kappa/2} * pi/sin(pi a_s) * I_{-a_s}(x) I_{-a_s}(xbar),   a_s = 2 i s mu,  x = 2 e^{kappa/4 - i theta/2}
   (vortex block f_a(z) = sum_n z^n/(n! Gamma(1+n-a)) = z^{a/2} I_{-a}(2 sqrt z), z = e^{kappa/2 - i theta}).
Claim B (Coulomb/LG side, W = 2 Lam cosh Y with twisted mass): thimble periods are K_a(x) (real-line thimble) and I_a(x);
   Z = e^{r kappa/2} [ 2 (K_a I_a' + I_a K_a') + (4 sin(pi a)/pi) K_a K_a' ],  ' = argument xbar, a = 2 i mu.
   This equals the Coulomb integral e^{r kappa/2} int drho e^{-2i mu rho} J0(4 e^{kappa/4} sqrt(sinh^2(rho/2)+sin^2(theta/2))).
So: residue (vacuum) basis {I_{-a}, I_{a}} and thimble basis {K_a, I_a} are related by the connection matrix
   I_{-a} = I_a + (2 sin(pi a)/pi) K_a  -- depends on mu only (not kappa, theta). The Hankel split of J0 is neither basis."""
import mpmath as mp
from cs_vacua import Zs
mp.mp.dps = 20

def coulomb(r, mu, kap, th):   # same quadrature as cp1_exact (u = sinh(rho/2))
    L = 4*mp.e**(kap/4); s2 = mp.sin(th/2)**2
    f = lambda u: 2*mp.cos(4*mu*mp.asinh(u))*mp.besselj(0, L*mp.sqrt(u*u+s2))*2/mp.sqrt(1+u*u)
    return mp.e**(r*kap/2)*mp.quadosc(f, [0, mp.inf], omega=L)

for r, mu, kap, th in [(0.3, 0.25, -2.0, 0.8), (0.4, 0.6, -1.0, 2.0), (0.2, 0.1, -3.0, 1.3), (0.35, 0.05, -0.5, 2.7)]:
    x = 2*mp.e**(kap/4 - 1j*th/2); xb = mp.conj(x); pre = mp.e**(r*kap/2)
    print(f"(r,mu,kappa,theta) = ({r},{mu},{kap},{th})")
    for s in (+1, -1):
        a = 2j*s*mu
        closed = pre*mp.pi/mp.sin(mp.pi*a)*mp.besseli(-a, x)*mp.besseli(-a, xb)
        num = Zs(s, r, mu, kap, th)
        print(f"   A  s={s:+d}: residues {mp.nstr(num, 12)}   closed form {mp.nstr(closed, 12)}   rel diff {mp.nstr(abs(num/closed-1), 3)}")
    a = 2j*mu
    K, I, Kb, Ib = mp.besselk(a, x), mp.besseli(a, x), mp.besselk(a, xb), mp.besseli(a, xb)
    thimble = pre*(2*(K*Ib + I*Kb) + 4*mp.sin(mp.pi*a)/mp.pi*K*Kb)
    C = coulomb(r, mu, kap, th)
    print(f"   B  thimble form {mp.nstr(thimble, 12)}   Coulomb integral {mp.nstr(C, 12)}   rel diff {mp.nstr(abs(thimble/C-1), 3)}")
    print(f"      K K' coefficient 4 sin(pi a)/pi = {mp.nstr(4*mp.sin(mp.pi*a)/mp.pi, 8)}  (= 4i sinh(2 pi mu)/pi)")
# mu -> 0: residue families blow up as 1/mu, the sum stays finite -> residue pieces are not observables of a single vacuum at mu=0
for mu in [1e-2, 1e-4]:
    zp, zm = Zs(+1, 0.3, mu, -2.0, 0.8), Zs(-1, 0.3, mu, -2.0, 0.8)
    print(f"mu={mu}: |Z_+| = {mp.nstr(abs(zp), 6)}   Z_+ + Z_- = {mp.nstr(zp+zm, 10)}")
