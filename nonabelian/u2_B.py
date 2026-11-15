"""Diagnostic 5: the w-charge analysis shows variant A has two gauge-invariant R=4h paths F1-L-L-F2 and F1-R-R-F2 where the
electric side has one monopole of R = 2(1-r) = 4h. Variant B closes only one path. Scan variant B's sign choices."""
import itertools, u2_defs as U
U.x = 0.07; U.M_mir = 2; U.Ng_mir = 12
for (yv, w) in U.pts:
    e = U.electric(yv, w)
    for s, s1, s3 in itertools.product((1, -1), (1, -1), (1, -1)):
        m = U.mirror(yv, w, 'B', s, s1, s3)
        print(f"w={w:.3f} B s={s} s1={s1} s3={s3}: e/m = {e/m:.6f}", flush=True)
