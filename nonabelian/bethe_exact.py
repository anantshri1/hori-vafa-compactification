"""Electric Bethe count by elimination. sigma=+1: z1 + w z2 P(z1) = 0, z2 + w z1 P(z2) = 0 (after dividing by z_a);
sigma=-1: z_a^{-2} -> multiply: 1 + w z1 z2 z_a^2 P(z_a) = 0. Count solutions with z1 != z2, z != 0, divide by 2."""
import sympy as sp, numpy as np
z1, z2 = sp.symbols('z1 z2')
rng = np.random.default_rng(7)
ys = [sp.Rational(int(a), 97) + sp.I*sp.Rational(int(b), 89) for a, b in rng.integers(-90, 90, (4, 2))]
w = sp.Rational(5, 7) + sp.I*sp.Rational(3, 11)
P = lambda z: sp.prod([1 - yj*z for yj in ys])
for sig, (E1, E2) in {+1: (z1 + w*z2*P(z1), z2 + w*z1*P(z2)), -1: (1 + w*z1*z2*z1**2*P(z1), 1 + w*z1*z2*z2**2*P(z2))}.items():
    R = sp.Poly(sp.resultant(sp.expand(E1), sp.expand(E2), z2), z1)
    roots = [complex(r) for r in sp.Poly(R, z1).nroots(n=30, maxsteps=200)]
    good = 0
    for r in roots:
        if abs(r) < 1e-8: continue
        c2 = [complex(c) for c in sp.Poly(sp.expand(E1.subs(z1, r)), z2).nroots(n=30, maxsteps=200)] if sp.Poly(sp.expand(E1.subs(z1, r)), z2).degree() > 0 else []
        for s in c2:
            if abs(complex(E2.subs({z1: r, z2: s}).evalf())) < 1e-6*(1 + abs(s))**6 and abs(s - r) > 1e-6*(1 + abs(r)) and abs(s) > 1e-8:
                good += 1
    print(f"sigma={sig:+d}: resultant degree {R.degree()}, ordered solutions z1!=z2: {good}, Bethe vacua {good/2}", flush=True)
