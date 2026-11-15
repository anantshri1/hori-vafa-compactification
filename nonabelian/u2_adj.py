"""Diagnostic 3: the mismatch scales like x^{~1.3}, the R-charge of the mirror adjoint (1+r). Scan the adjoint's R-charge."""
import u2_defs as U
U.x = 0.07; U.M_mir = 2; U.Ng_mir = 12
r, h = U.r, U.h
yv, w = U.pts[0]
e = U.electric(yv, w)
for name, Ra in [('1+r', 1+r), ('2r', 2*r), ('1', 1.0), ('2', 2.0), ('2-2r', 2-2*r), ('1+r-h', 1+r-h), ('1+r+h', 1+r+h), ('r', r)]:
    U.Radj = Ra
    m = U.mirror(yv, w, 'A', 1, 1, 1)
    print(f"R_adj={name:7s}: |e/m-1| = {abs(e/m-1):.3e}", flush=True)
