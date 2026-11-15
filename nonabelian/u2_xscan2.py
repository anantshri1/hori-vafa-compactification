"""Diagnostic 2: flip the eta-fugacity (w^{-1} -> w), and the FI orientation, at x = 0.07."""
import numpy as np, u2_defs as U
U.x = 0.07; U.M_mir = 2; U.Ng_mir = 12
for yv, w in U.pts:
    e = U.electric(yv, w)
    for eta in (-1, 1):
        U.ETA = eta
        for s in (1, -1):
            m = U.mirror(yv, w, 'A', s, 1, 1)
            print(f"w={w:.3f} eta^{eta} s={s}: e/m = {e/m:.6f}", flush=True)
