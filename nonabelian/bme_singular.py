import sympy as sp, numpy as np, sys, subprocess
import optc_defs as O
z = sp.symbols('zL zT zB zR'); s = sp.Symbol('s')
Qs = [c[0] for c in [(O.one(O.B, 1), {}, 0)] + O.rest]
for sig in (+1, -1):
    for seed in (3, 8):
        K = np.array(O.Kmat(1), float); Kp = np.round(K + sig*0.5*np.array(Qs, float).T @ np.array(Qs, float)).astype(int)
        rng = np.random.default_rng(seed)
        R = lambda: sp.Rational(int(rng.integers(2, 19)), int(rng.integers(2, 19)))*(1 if rng.random() < .5 else -1)
        f = [R() for _ in Qs]; t = [R() for _ in range(4)]
        X = [f[i]*sp.Mul(*[z[a]**Qs[i][a] for a in range(4)]) for i in range(len(Qs))]
        eqs = []
        for a in range(4):
            lhs = sp.Mul(*[z[b]**int(Kp[a][b]) for b in range(4)])*sp.Mul(*[(1 - X[i])**(-Qs[i][a]) for i in range(len(Qs))])
            eqs.append(sp.expand(sp.fraction(sp.together(lhs - t[a]))[0]))
        sat = sp.expand(s*sp.Mul(*z)*sp.Mul(*[sp.fraction(sp.together(1 - x))[0] for x in X]) - 1)
        txt = "ring r = 0,(s,zL,zT,zB,zR),dp;\nideal I = " + ",\n".join(str(e).replace('**', '^') for e in eqs + [sat]) + ";\noption(redSB);\nideal G = std(I);\nprint(vdim(G));\nquit;\n"
        open('bme.sing', 'w').write(txt)
        out = subprocess.run(['Singular', '-q', 'bme.sing'], capture_output=True, text=True, timeout=1200).stdout.strip()
        print(f"sigma={sig:+d} seed={seed}: K'={Kp.tolist()}  #mirror Bethe vacua = {out}", flush=True)
