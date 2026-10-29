"""Check that the residue sum in cp1_exact.py IS the Benini-Cremonesi / Doroud et al. S^2 contour integral
  Z = sum_m e^{-i theta m} int_R dsigma/(2pi) e^{i kappa sigma} prod_{s=+-} G(r/2 - i(sigma+s mu) - m/2)/G(1 - r/2 + i(sigma+s mu) - m/2)
(kappa = -4 pi xi, twisted masses +-mu, R-charge r). Direct quadrature in sigma (quadosc: plain quad fails, 0.2-30% error), flux sum truncated."""
import mpmath as mp
from cp1_exact import glsm
mp.mp.dps = 15
def contour(r, mu, kap, th, MM=None):
    MM = MM or int(30/abs(kap)) + 2   # flux tail ~ e^{kappa |m|/2}
    def Zm(m):
        f = lambda sg: mp.e**(1j*kap*sg)*mp.fprod(mp.gamma(r/2-1j*(sg+s*mu)-m/2)/mp.gamma(1-r/2+1j*(sg+s*mu)-m/2) for s in (1,-1))
        w = abs(kap)  # integrand decays only like |sigma|^{2(r-1)} and oscillates: use oscillatory quadrature
        return (mp.quad(f, [-3, 0, 3]) + mp.quadosc(f, [3, mp.inf], omega=w) + mp.quadosc(f, [-mp.inf, -3], omega=w))/(2*mp.pi)
    return Zm(0) + mp.fsum(2*mp.cos(th*m)*Zm(m) for m in range(1, MM+1))
for p in [(0.3, 0.25, -2.0, 0.8), (0.4, 0.6, -1.0, 2.0)]:
    C, G = contour(*p), glsm(*p)
    print(p, " contour:", mp.nstr(C, 12), "  residues:", mp.nstr(G, 12), "  rel diff:", mp.nstr(abs(C/G-1), 3), flush=True)
