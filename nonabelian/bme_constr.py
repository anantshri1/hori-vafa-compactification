"""LOG 44: mirror Bethe count with the superpotential constraints imposed on the parameters:
triangles: f_LB f_TL f_BT = 1, f_RB f_TR f_BT = 1 (W = cubic terms, R = 2 each => fugacity product 1);
monopole / column structure: t_T = s * t_B (same column FI y3/y2; sign s to be tested). Symmetric quantisation (sigma=+1)."""
import sympy as sp, numpy as np, subprocess, itertools, sys
import optc_defs as O
z = sp.symbols('zL zT zB zR'); s_ = sp.Symbol('s')
Qs = [c[0] for c in [(O.one(O.B, 1), {}, 0)] + O.rest]     # B, L-B, T-L, -T, R-B, T-R, B-T
Qa = np.array(Qs, float); K = np.array(O.Kmat(1), float)
def count(f, t, sgv):
    Kp = np.round(K + 0.5*sum(sgv[i]*np.outer(Qa[i], Qa[i]) for i in range(7))).astype(int)
    X = [f[i]*sp.Mul(*[z[a]**Qs[i][a] for a in range(4)]) for i in range(7)]
    eqs = [sp.expand(sp.fraction(sp.together(sp.Mul(*[z[b]**int(Kp[a][b]) for b in range(4)])*sp.Mul(*[(1 - X[i])**(-Qs[i][a]) for i in range(7)]) - t[a]))[0]) for a in range(4)]
    sat = sp.expand(s_*sp.Mul(*z)*sp.Mul(*[sp.fraction(sp.together(1 - x))[0] for x in X]) - 1)
    open('bmc.sing', 'w').write("ring r = 0,(s,zL,zT,zB,zR),dp;\nideal I = " + ",\n".join(str(e).replace('**', '^') for e in eqs + [sat]) + ";\nideal G = std(I);\nprint(vdim(G));\nquit;\n")
    return subprocess.run(['Singular', '-q', 'bmc.sing'], capture_output=True, text=True, timeout=500).stdout.strip()
rng = np.random.default_rng(int(sys.argv[1]) if len(sys.argv) > 1 else 3)
R = lambda: sp.Rational(int(rng.integers(2, 19)), int(rng.integers(2, 19)))*(1 if rng.random() < .5 else -1)
fB, fLB, fTL, fRB, f_T = R(), R(), R(), R(), R()
fBT = 1/(fLB*fTL); fTR = 1/(fRB*fBT)
f = [fB, fLB, fTL, f_T, fRB, fTR, fBT]
tL, tB, tR = R(), R(), R()
for sgn in (1, -1):
    t = [tL, sgn*tB, tB, tR]   # order (L,T,B,R)
    print(f"constraints on, t_T = {sgn:+d} t_B, uniform sigma=+1: #vacua =", count(f, t, [1]*7), flush=True)
print("constraints on, t_T free, uniform sigma=+1: #vacua =", count(f, [tL, R(), tB, tR], [1]*7), flush=True)
print("constraints on, uniform sigma=-1: #vacua =", count(f, [tL, tB, tB, tR], [-1]*7), flush=True)
print("constraints OFF (generic f), uniform sigma=+1: #vacua =", count([R() for _ in range(7)], [tL, tB, tB, tR], [1]*7), flush=True)
