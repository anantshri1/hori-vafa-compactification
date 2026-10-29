"""Diagnosis of the cut-independent ~5e-7 per-sector error in cp1_contour_sector.py.
The sigma-integrand is NOT a fixed-frequency oscillation: by Stirling the Gamma ratio carries the phase ~ -4 sigma log sigma
(two chirals), so the total phase is g(sigma) = kappa sigma - 4 sigma log sigma + 4 sigma (a chirp). quadosc with omega = |kappa|
uses the wrong period. Here: f(-sigma) = conj f(sigma) (mu -> -mu symmetry of the pair), so int_R = 2 Re int_0^inf, and the tail
uses quadosc with zeros g(sigma_n) = g(cut) - n pi."""
import mpmath as mp, sys, time
from cp1_contour_sector import resid
mp.mp.dps = 20
r, mu, kap = 0.4, 0.6, -1.0
m = int(sys.argv[1]); cut = 6
f = lambda sg: mp.e**(1j*kap*sg)*mp.fprod(mp.gamma(r/2-1j*(sg+s*mu)-m/2)/mp.gamma(1-r/2+1j*(sg+s*mu)-m/2) for s in (1,-1))
g = lambda s: kap*s - 4*s*mp.log(s) + 4*s
g0 = g(cut)
def zero(n):
    return mp.findroot(lambda s: g(s) - (g0 - n*mp.pi), cut + n*mp.pi/(4*mp.log(cut + n + 1)))
t = time.time()
I = 2*mp.re(mp.quad(f, mp.linspace(0, cut, 13)) + mp.quadosc(f, [cut, mp.inf], zeros=zero))/(2*mp.pi)
R = resid(m)
print(f"m={m}: chirp-aware direct {mp.nstr(I, 15)}  residues {mp.nstr(mp.re(R), 15)}  |diff| {mp.nstr(abs(I-R), 3)}  ({time.time()-t:.0f}s)", flush=True)
