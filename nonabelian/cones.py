"""Power counting of the 2d Coulomb integrand of the reduced mirror quiver (B chiral removed, option c).
For every subspace V cut out by setting a subset of chiral weights to zero, generic |f| ~ |s|^{sum_{charged}(R_i-1)};
convergence needs exponent < -dim V. Print all divergent subspaces as a function of r. Also with B chiral flat (same thing)."""
import numpy as np, itertools, sympy as sp
r = sp.symbols('r')
# charges on (L,T,B,R), R-charge
CH = {'L-B': ((1,0,-1,0), r), 'T-L': ((-1,1,0,0), r), '-T': ((0,-1,0,0), r),
      'R-B': ((0,0,-1,1), r), 'T-R': ((0,1,0,-1), r), 'B-T': ((0,-1,1,0), 2-2*r)}
names = list(CH)
seen = {}
for k in range(len(names)+1):
    for zero in itertools.combinations(names, k):
        A = sp.Matrix([CH[n][0] for n in zero]) if zero else sp.zeros(0, 4)
        V = A.nullspace() if zero else [sp.eye(4)[:, i] for i in range(4)]
        d = len(V)
        if d == 0: continue
        Vm = sp.Matrix.hstack(*V)
        charged = [n for n in names if any(x != 0 for x in (sp.Matrix([CH[n][0]])*Vm))]
        key = tuple(sorted(charged)) + (d,)
        if key in seen: continue
        p = sp.simplify(sum(CH[n][1]-1 for n in charged))
        margin = sp.simplify(p + d)   # need margin < 0
        seen[key] = (d, charged, p, margin, [list(v) for v in V])
for key, (d, ch, p, mg, V) in sorted(seen.items(), key=lambda t: t[1][0]):
    bad = [rv for rv in (0.3, 0.0, -0.3) if float(mg.subs(r, rv)) >= 0]
    print(f"dim {d}: span {V}  charged {ch}  exponent {p}  margin {mg}  divergent at r in {bad}")
