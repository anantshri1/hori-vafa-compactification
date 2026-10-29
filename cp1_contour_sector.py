"""Sector-by-sector version of the contour check at the hard point (r,mu,kappa,theta) = (0.4,0.6,-1,2).
The full check (cp1_contour_hp.py) hit the 2 h limit. Here: compare the direct sigma-integral and the residue sum for single
fluxes m = 0, 1, 5, at two quadrature settings, to separate quadrature error from a real discrepancy."""
import mpmath as mp, sys, time
mp.mp.dps = 20
r, mu, kap = 0.4, 0.6, -1.0
def resid(m, K=80):
    return mp.fsum((-1)**k/(mp.factorial(k)*mp.factorial(k+abs(m))) * mp.e**(kap*(r/2+abs(m)/2+k)) * mp.e**(-1j*s*mu*kap)
                   * mp.gamma(-k+2j*s*mu)/mp.gamma(1+abs(m)+k-2j*s*mu) for s in (1,-1) for k in range(K))
def direct(m, cut):
    f = lambda sg: mp.e**(1j*kap*sg)*mp.fprod(mp.gamma(r/2-1j*(sg+s*mu)-m/2)/mp.gamma(1-r/2+1j*(sg+s*mu)-m/2) for s in (1,-1))
    return (mp.quad(f, mp.linspace(-cut, cut, 2*cut+1)) + mp.quadosc(f, [cut, mp.inf], omega=abs(kap))
            + mp.quadosc(f, [-mp.inf, -cut], omega=abs(kap)))/(2*mp.pi)
if __name__ == "__main__":
  m = int(sys.argv[1]); cut = int(sys.argv[2])
  t = time.time(); R = resid(m); D = direct(m, cut)
  print(f"m={m} cut={cut}: direct {mp.nstr(D, 15)}  residues {mp.nstr(R, 15)}  |diff| {mp.nstr(abs(D-R), 3)}  ({time.time()-t:.0f}s)", flush=True)
