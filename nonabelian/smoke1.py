"""SMOKE 1: Fourier-convolution sectors vs the validated grid engine (idx_jax, per_sector) at small flux."""
import numpy as np, optc_defs as O, idx_jax, fourier_sector as FS
x, n = 0.3, 4; fl = dict(O.fl); fl['eta'] = np.exp(-0.9j)
ch = [(O.one(O.B, 1), {'eta': -1}, O.r)] + O.rest
rng = np.random.default_rng(0); secs = [tuple(int(v) for v in rng.integers(-3, 3, 4)) for _ in range(6)] + [(0, 0, -4, 0), (-1, 2, -4, -1)]
g = idx_jax.theory_index(4, ch, O.Kmat(1), {(O.B, 'eta'): -.5}, O.FIfull, fl, x, M=0, Ngrid=int(__import__("sys").argv[1]), bg={'eta': -n}, sectors=secs, per_sector=True)
for s, gv in zip(secs, g):
    fv = FS.sector(s, x, fl, n)
    print(f"sector {s}: grid {gv:.10e}  fourier {fv:.10e}  rel {abs(fv/gv - 1) if abs(gv) > 1e-14 else abs(fv-gv):.1e}", flush=True)
