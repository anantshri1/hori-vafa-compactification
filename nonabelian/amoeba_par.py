import sys, numpy as np
from multiprocessing import Pool
import amoeba as A, optc_defs as O
n, x, mB, lo, hi = int(sys.argv[1]), float(sys.argv[2]), int(sys.argv[3]), int(sys.argv[4]), int(sys.argv[5])
fl = dict(O.fl); fl['eta'] = np.exp(-0.9j)
def job(mL):
    rng = np.random.default_rng(mL + 100)
    ok = sum(A.residual((mL, mT, mB, mR), x, n, fl, starts=8, rng=rng) < 1e-3 for mT in (0, 2, 4) for mR in range(lo, hi + 1))
    return mL, ok/(3*(hi - lo + 1))
with Pool() as p: out = p.map(job, range(lo, hi + 1))
print(f"n={n} x={x} m_B={mB}: solvable fraction vs m_L: " + " ".join(f"{l}:{f:.2f}" for l, f in out), flush=True)
