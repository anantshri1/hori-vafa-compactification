"""
Step 0: small-circle asymptotics of a single 3d N=2 chiral, in Anant's index
convention (3dSCI notebook, `Chir`):

  Chir[zeta, nu, R] = (zeta x^(R-1))^(-|nu|/2)
       * (zeta^-1 (-1)^nu x^(2-R+|nu|); x^2)_inf / (zeta (-1)^nu x^(R+|nu|); x^2)_inf

Conventions used here:  x = e^{-beta},  q = x^2 = e^{-eps},  eps = 2 beta.
(a; q)_inf = prod_{k>=0} (1 - a q^k)   (Mathematica QPochhammer / mpmath qp).

Every asymptotic formula below is checked against exact high-precision numerics.
"""
import mpmath as mp
mp.mp.dps = 40

def qp(a, q):
    return mp.qp(a, q)

def report(label, exact, pred, eps):
    r = exact / pred
    print(f"  {label:<28} eps={float(eps):<7} |ratio-1| = {float(abs(r - 1)):.3e}")
    return r

eps_list = [mp.mpf('0.2'), mp.mpf('0.1'), mp.mpf('0.05'), mp.mpf('0.025')]

# ---------------------------------------------------------------- Test A
print("A. (q^a; q)_inf  ~  sqrt(2pi) e^{-pi^2/(6 eps)} eps^{1/2-a} / Gamma(a)")
for a in [mp.mpc('0.3', '0.4'), mp.mpc('1.7', '-0.25')]:
    print(f" a = {a}")
    for eps in eps_list:
        q = mp.e**(-eps)
        pred = mp.sqrt(2*mp.pi) * mp.e**(-mp.pi**2/(6*eps)) * eps**(mp.mpf(1)/2 - a) / mp.gamma(a)
        report("", qp(q**a, q), pred, eps)

# ---------------------------------------------------------------- Test B
print("\nB. (-q^a; q)_inf  ~  e^{+pi^2/(12 eps)} 2^{1/2-a}   (sign inside: no Gamma)")
for a in [mp.mpc('0.3', '0.4')]:
    print(f" a = {a}")
    for eps in eps_list:
        q = mp.e**(-eps)
        pred = mp.e**(mp.pi**2/(12*eps)) * mp.mpf(2)**(mp.mpf(1)/2 - a)
        report("", qp(-q**a, q), pred, eps)

# ---------------------------------------------------------------- Test C
def chir_ratio(zeta, nu, R, x):
    """Pochhammer ratio of Chir (zero-point factor excluded)."""
    s = (-1)**nu
    return qp(s * x**(2 - R + abs(nu)) / zeta, x**2) / qp(s * zeta * x**(R + abs(nu)), x**2)

R, sig = mp.mpf('0.6'), mp.mpf('0.3')
print("\nC. Chir Pochhammer ratio with zeta = (-1)^nu q^{-i sigma} (sign absorbed)")
print("   prediction: (1-q)^{R-1-2i sigma} Gamma(R/2 - i sigma + |nu|/2) / Gamma(1 - R/2 + i sigma + |nu|/2)")
for nu in [0, 1, -1, 2, -3]:
    print(f" nu = {nu}")
    for eps in eps_list:
        q = mp.e**(-eps); x = mp.e**(-eps/2)
        zeta = (-1)**nu * q**(-1j*sig)
        pred = (1 - q)**(R - 1 - 2j*sig) * mp.gamma(R/2 - 1j*sig + abs(nu)/mp.mpf(2)) \
               / mp.gamma(1 - R/2 + 1j*sig + abs(nu)/mp.mpf(2))
        report("", chir_ratio(zeta, nu, R, x), pred, eps)

# ---------------------------------------------------------------- Test D
print("\nD. Odd nu but zeta near +1 (sign NOT absorbed): ratio -> 2^{R-1-2i sigma}, finite")
print("   so relative to the even/absorbed sectors (~eps^{R-1}) it is suppressed by eps^{1-R}")
for nu in [1]:
    for eps in eps_list:
        q = mp.e**(-eps); x = mp.e**(-eps/2)
        zeta = q**(-1j*sig)
        pred = mp.mpf(2)**(R - 1 - 2j*sig)
        report("", chir_ratio(zeta, nu, R, x), pred, eps)

# ---------------------------------------------------------------- Test E
print("\nE. Zero-point factor (zeta x^{R-1})^{-|nu|/2}, principal branch, across zeta = -1")
x = mp.e**(-mp.mpf('0.05'))
for nu in [1, 2]:
    for d in [mp.mpf('1e-8'), mp.mpf('-1e-8')]:
        zeta = -mp.e**(1j*d)
        zp = (zeta * x**(R - 1))**(-abs(nu)/mp.mpf(2))
        print(f"  nu={nu}  arg(zeta) just {'below' if d>0 else 'above'} -pi/+pi edge: "
              f"zero-point = {mp.nstr(zp, 8)}")
print("  -> for odd |nu| the factor flips sign across zeta = -1: discontinuous exactly at the")
print("     point where odd-flux sectors localise (Test C). Paired with a level-1/2 CS factor")
print("     zeta^{-nu/2} on the same principal branch, the product zeta^{-(|nu|+nu)/2} is single-valued.")
for nu in [1, -1, 3]:
    vals = []
    for d in [mp.mpf('1e-8'), mp.mpf('-1e-8')]:
        zeta = -mp.e**(1j*d)
        vals.append(zeta**(-abs(nu)/mp.mpf(2)) * zeta**(-nu/mp.mpf(2)))
    print(f"  nu={nu:>2}: chiral x CS_(1/2) half-powers, both sides of the cut: "
          f"{mp.nstr(vals[0], 6)} | {mp.nstr(vals[1], 6)}")
