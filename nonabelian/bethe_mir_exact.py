"""LOG 43: exact mirror Bethe-vacuum count. Node equations prod_b z_b^{K'_ab} prod_i (1 - X_i)^{-Q_ia} = t_a,
K' = K + sig/2 sum_i Q_i Q_i^T, X_i = f_i z^{Q_i}. Polynomialise, saturate by z_a and (1 - X_i), count solutions
(degree of the zero-dimensional ideal) with a Groebner basis over QQ (random rational parameters)."""
import sympy as sp, numpy as np, sys
import optc_defs as O
sig = int(sys.argv[1])
z = sp.symbols('zL zT zB zR'); s = sp.Symbol('s')
Qs = [c[0] for c in [(O.one(O.B, 1), {}, 0)] + O.rest]
K = np.array(O.Kmat(1), float); Kp = K + sig*0.5*np.array(Qs, float).T @ np.array(Qs, float)
Kp = np.round(Kp).astype(int)
rng = np.random.default_rng(3)
R = lambda: sp.Rational(int(rng.integers(2, 19)), int(rng.integers(2, 19)))*(1 if rng.random() < .5 else -1)
f = [R() for _ in Qs]; t = [R() for _ in range(4)]
X = [f[i]*sp.Mul(*[z[a]**Qs[i][a] for a in range(4)]) for i in range(len(Qs))]
eqs = []
for a in range(4):
    lhs = sp.Mul(*[z[b]**int(Kp[a][b]) for b in range(4)])*sp.Mul(*[(1 - X[i])**(-Qs[i][a]) for i in range(len(Qs))])
    num, den = sp.fraction(sp.together(lhs - t[a]))
    eqs.append(sp.expand(num))
sat = sp.expand(s*sp.Mul(*z)*sp.Mul(*[sp.fraction(sp.together(1 - x))[0] for x in X]) - 1)
G = sp.groebner(eqs + [sat], s, *z, order='grevlex')
# count solutions: number of standard monomials
def count(G, gens):
    import itertools
    lead = [sp.Poly(g, *gens).monoms(order='grevlex')[0] for g in G.exprs]
    # enumerate standard monomials up to a bound
    bound = 30; cnt = 0
    for e in itertools.product(range(bound), repeat=len(gens)):
        if not any(all(e[k] >= l[k] for k in range(len(gens))) for l in lead): cnt += 1
    return cnt
print(f"sigma={sig:+d}: K' =", Kp.tolist(), " #Bethe vacua =", count(G, (s,)+z), flush=True)
