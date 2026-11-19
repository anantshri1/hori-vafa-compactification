import sympy as sp, numpy as np, sys, subprocess, itertools, time, os
from multiprocessing import Pool
import optc_defs as O
z = sp.symbols('zL zT zB zR'); s = sp.Symbol('s')
Qs = [c[0] for c in [(O.one(O.B, 1), {}, 0)] + O.rest]
nm = ['B', 'L-B', 'T-L', '-T', 'R-B', 'T-R', 'B-T']
rng = np.random.default_rng(3)
R = lambda: sp.Rational(int(rng.integers(2, 19)), int(rng.integers(2, 19)))*(1 if rng.random() < .5 else -1)
f = [R() for _ in Qs]; t = [R() for _ in range(4)]
X = [f[i]*sp.Mul(*[z[a]**Qs[i][a] for a in range(4)]) for i in range(len(Qs))]
K = np.array(O.Kmat(1), float); Qa = np.array(Qs, float)
def run(sg):
    Kp = K + 0.5*sum(sg[i]*np.outer(Qa[i], Qa[i]) for i in range(7))
    if not np.allclose(Kp, np.round(Kp)): return sg, None
    Kp = np.round(Kp).astype(int); eqs = []
    for a in range(4):
        lhs = sp.Mul(*[z[b]**int(Kp[a][b]) for b in range(4)])*sp.Mul(*[(1 - X[i])**(-Qs[i][a]) for i in range(7)])
        eqs.append(sp.expand(sp.fraction(sp.together(lhs - t[a]))[0]))
    sat = sp.expand(s*sp.Mul(*z)*sp.Mul(*[sp.fraction(sp.together(1 - x))[0] for x in X]) - 1)
    fn = f"bme_{os.getpid()}.sing"
    open(fn, 'w').write("ring r = 0,(s,zL,zT,zB,zR),dp;\nideal I = " + ",\n".join(str(e).replace('**', '^') for e in eqs + [sat]) + ";\nideal G = std(I);\nprint(vdim(G));\nquit;\n")
    try: out = subprocess.run(['Singular', '-q', fn], capture_output=True, text=True, timeout=300).stdout.strip()
    except subprocess.TimeoutExpired: out = 'timeout'
    return sg, out
if __name__ == "__main__":
    sgs = list(itertools.product((1, -1), repeat=7))
    sgs = [g for g in sgs if (g[1], g[2]) != (g[4], g[5])]     # non-symmetric remainder
    t0 = time.time()
    with Pool(4) as p:
        for sg, out in p.imap_unordered(run, sgs):
            print(dict(zip(nm, sg)), "->", out, f"[{time.time()-t0:.0f}s]", flush=True)
