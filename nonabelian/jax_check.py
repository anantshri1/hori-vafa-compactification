import numpy as np, time, optc_defs as O, idx, idx_jax
ch = [(O.one(O.B, 1), {'eta': -1}, O.r)] + O.rest
fl = dict(O.fl); fl['eta'] = np.exp(-0.9j)
for x, bgn, sh in [(0.3, -4, None), (0.3, -4, [0, 0, -4, 0])]:
    a = (4, ch, O.Kmat(1), {(O.B, 'eta'): -.5}, O.FIfull, fl, x)
    t = time.time(); I1 = idx.theory_index(*a, M=1, Ngrid=12, bg={'eta': bgn}, mshift=sh); t1 = time.time() - t
    t = time.time(); I2 = idx_jax.theory_index(*a, M=1, Ngrid=12, bg={'eta': bgn}, mshift=sh); t2 = time.time() - t
    print(f"x={x} shift={sh}: numpy {I1:.10e} ({t1:.1f}s)  jax {I2:.10e} ({t2:.1f}s)  rel {abs(I1/I2-1):.1e}", flush=True)
