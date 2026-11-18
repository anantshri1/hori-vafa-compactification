"""LOG 34: single chiral in the Bessel regime x^{nu} ~ eps. Prediction:
Chir(zeta, nu, R) ~ zeropoint * exp( s (zeta x^{R+nu} - x^{2-R+nu}/zeta)/(1-q) ),  s = (-1)^nu,
i.e. W - Wbar with w = x^nu zeta: X+ = x^R w, X- = x^{2-R} wbar."""
import numpy as np
from idx import qpoch
R = 0.3
for eps in (0.2, 0.1, 0.05, 0.025):
    x = np.exp(-eps/2); q = x*x
    nu = int(round(2*np.log(1/eps)/eps))           # x^nu ~ eps
    errs = []
    for th in np.linspace(0, 2*np.pi, 7)[:-1]:
        s = (-1)**nu; zeta = s*np.exp(1j*th)
        ex = qpoch(x**(2 - R + nu)/zeta, q)/qpoch(zeta*x**(R + nu), q)
        pr = np.exp((zeta*x**(R + nu) - x**(2 - R + nu)/zeta)/(1 - q))
        errs.append(abs(ex/pr - 1))
    print(f"eps={eps} nu={nu} x^nu/eps={x**nu/eps:.3f}: max rel err {max(errs):.2e}", flush=True)
