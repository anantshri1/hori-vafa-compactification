"""Direct S^2 contour integral (chirp-aware tail, see cp1_contour_chirp.py) vs residue sum, all 4 test points, full flux sum."""
import mpmath as mp
from cp1_exact import glsm
mp.mp.dps = 20
def contour(r, mu, kap, th, cut=6):
    MM = int(40/abs(kap)) + 2
    g = lambda s: kap*s - 4*s*mp.log(s) + 4*s
    g0 = g(cut)
    zero = lambda n: mp.findroot(lambda s: g(s) - (g0 - n*mp.pi), cut + n*mp.pi/(4*mp.log(cut + n + 1)))
    def Zm(m):
        f = lambda sg: mp.e**(1j*kap*sg)*mp.fprod(mp.gamma(r/2-1j*(sg+s*mu)-m/2)/mp.gamma(1-r/2+1j*(sg+s*mu)-m/2) for s in (1,-1))
        return 2*mp.re(mp.quad(f, mp.linspace(0, cut, 13)) + mp.quadosc(f, [cut, mp.inf], zeros=zero))/(2*mp.pi)
    return Zm(0) + mp.fsum(2*mp.cos(th*m)*Zm(m) for m in range(1, MM+1))
for p in [(0.3, 0.25, -2.0, 0.8), (0.4, 0.6, -1.0, 2.0), (0.2, 0.1, -3.0, 1.3), (0.45, 0.8, -0.5, 2.7)]:
    C, G = contour(*p), glsm(*p)
    print(p, " contour:", mp.nstr(C, 15), "  residues:", mp.nstr(mp.re(G), 15), "  rel diff:", mp.nstr(abs(C/G-1), 3), flush=True)
