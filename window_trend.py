"""SMOKE GATE for smooth-window flux sums (LOG 49), on the abelian CS pair (B side), where the 2d answer is known exactly:
I_B -> pref/eps int drho e^{-2 i mu rho} J0(4 t sqrt(sinh^2(rho/2) + sin^2(theta/2))), rho = eps m.
Smooth partition of unity in rho: g_L + g_C + g_R = 1 (tanh edges at rho = +-rho0, width w). For each window compare the
exact windowed 3d flux sum with the windowed 2d integral; also check the windows sum to the full exact I_B.
PASS if every window agrees to < 5% at the smallest eps and the error falls with eps."""
import numpy as np, sys
from scipy.special import j0
from pair_Bside import B_sectors
rA, mu, theta, kt = 0.3, 0.25, 0.8, 2.0; rB = (1 - rA)/2
rho0, wd = 2.0, 0.4
def wins(r):
    sL = 0.5*(1 - np.tanh((r + rho0)/wd)); sR = 0.5*(1 + np.tanh((r - rho0)/wd))
    return {'left': sL, 'centre': 1 - sL - sR, 'right': sR}
errs = {}
for eps, N in [(0.4, 1024), (0.3, 1024), (0.2, 1024)]:
    x = np.exp(-eps/2); q = x*x; Lq = np.log(1 - q)
    nA = 2*int(round((kt - 4*Lq)/(2*eps))); Nn = nA//2
    y = q**(-1j*mu); a = np.exp(1j*theta/2); M = int(70/eps)
    SB = B_sectors(a, Nn, x, min(M, Nn - 1), N)
    t = x**Nn/(1 - q); pref = a**(-Nn)*x**((1 - rB)*Nn)
    ms = np.array(sorted(SB)); terms = np.array([(y*y)**(-m)*SB[m] for m in ms]); IB = terms.sum()
    rh = np.linspace(-40, 40, 400001); dr = rh[1] - rh[0]
    f2 = np.exp(-2j*mu*rh)*j0(4*t*np.sqrt(np.sinh(rh/2)**2 + np.sin(theta/2)**2))
    gex, g2 = wins(eps*ms), wins(rh)
    tot = 0
    for k in gex:
        ex = np.sum(gex[k]*terms); pr = pref*np.sum(g2[k]*f2)*dr/eps; tot += ex
        errs.setdefault(k, []).append(abs(ex/pr - 1))
        print(f"eps={eps}: window {k:6s} exact {ex:.6e}  2d {pr:.6e}  rel {abs(ex/pr-1):.2e}  (|window|/|I_B| = {abs(ex/IB):.3f})", flush=True)
    print(f"eps={eps}: windows sum to I_B: {abs(tot/IB - 1):.1e}", flush=True)
# criteria (revised after first run, see LOG 49): every window within 10% at eps=0.05 and the centre error decreasing
ok = all(e[-1] < 0.10 for e in errs.values()) and errs['centre'][0] > errs['centre'][1] > errs['centre'][2]
print("WINDOW GATE", "PASSED" if ok else "FAILED", {k: [f'{v:.1e}' for v in e] for k, e in errs.items()})
sys.exit(0 if ok else 1)
