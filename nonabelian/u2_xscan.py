"""Diagnostic: ratio electric/mirror (variant A, s=s1=s3=+1) as x varies; (ratio-1) ~ x^alpha identifies the R-charge of the
first mismatched operator."""
import numpy as np, u2_defs as U
yv, w = U.pts[0]
for x in [0.04, 0.07, 0.1]:
    U.x = x
    U.M_mir = 2; U.Ng_mir = 12
    import idx
    e = U.electric(yv, w); m = U.mirror(yv, w, 'A', 1, 1, 1)
    print(x, e, m, e/m, abs(e/m-1), flush=True)
