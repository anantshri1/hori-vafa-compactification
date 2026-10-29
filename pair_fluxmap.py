"""
Flux-dependent map: A with U(1)_T background flux nA  vs  B with axial background flux nax.
A: Z_cl = wA^{-m} z^{-m - nA}  (k = 1), zero-points z^{-|m|} x^{(1-rA)|m|}.
B: Q (fug z a, flux m + nax), Qt (fug a/z, flux -m + nax), CS z^{+m}, zero-points principal branch
   (z a x^{rB-1})^{-|m+nax|/2} (a z^{-1} x^{rB-1})^{-|m-nax|/2}  -> grouped analytically (integer powers).
Fit: I_A(y, wA, nA) / I_B(a = wA^{1/2}, wB = y^2, nax) = monomial c * y^alpha * wA^beta * x^gamma ?
"""
import sys, numpy as np
sys.path.insert(0, '.')
from index_lib import chir_ratio

x, N, M = 0.3, 256, 40
z = np.exp(2j*np.pi*(np.arange(N)+0.5)/N)
rA = 0.3; rB = (1 - rA)/2

def IA(y, wA, nA):
    tot = 0
    for m in range(-M, M+1):
        f = chir_ratio(z*y, m, rA, x) * chir_ratio(z/y, m, rA, x)
        tot += wA**(-m) * x**((1-rA)*abs(m)) * np.mean(z**(-(m + abs(m)) - nA) * f)
    return tot

def IB(a, wB, nax):
    tot = 0
    for m in range(-M, M+1):
        f1, f2 = m + nax, -m + nax
        zp = (-abs(f1) + abs(f2)) // 2          # z-power from zero-points (integer: |f1|,|f2| same parity)
        ap = -(abs(f1) + abs(f2)) // 2
        f = chir_ratio(z*a, f1, rB, x) * chir_ratio(a/z, f2, rB, x)
        tot += wB**(-m) * a**ap * x**((1-rB)*(abs(f1)+abs(f2))/2) * np.mean(z**(zp + m) * f)
    return tot

pts = [(np.exp(0.5j), np.exp(0.8j)), (np.exp(1.4j), np.exp(-0.6j)), (np.exp(2.3j), np.exp(1.9j)), (np.exp(-1.1j), np.exp(2.6j))]
for nA in [2, -2, 4, 1]:
    for nax_rule, nax in [("nA/2", nA/2), ("nA", nA), ("-nA/2", -nA/2)]:
        if nax != int(nax):
            continue
        nax = int(nax)
        R = np.array([IA(y, w, nA) / IB(w**0.5, y**2, nax) for y, w in pts])
        th_y = np.array([np.angle(y) for y, _ in pts]); th_w = np.array([np.angle(w) for _, w in pts])
        # fit arg R = alpha th_y + beta th_w + c  (least squares), check |R| constant
        Amat = np.c_[th_y, th_w, np.ones(len(pts))]
        ph = np.angle(R)
        best = None
        for k in np.array(np.meshgrid(*[range(-2, 3)]*4)).T.reshape(-1, 4):   # unwrap ambiguity
            sol, res, *_ = np.linalg.lstsq(Amat, ph + 2*np.pi*k, rcond=None)
            r = np.linalg.norm(Amat @ sol - ph - 2*np.pi*k)
            if best is None or r < best[0]:
                best = (r, sol)
        print(f"nA={nA:+d} nax={nax:+d} ({nax_rule}): |R| = {np.round(abs(R), 6)}  "
              f"phase fit: y^{best[1][0]:+.3f} wA^{best[1][1]:+.3f} (resid {best[0]:.1e})")
