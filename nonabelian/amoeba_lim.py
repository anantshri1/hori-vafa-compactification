"""LOG 27 (i): predicted L-support in u = eps m_L at fixed n*eps = 4 as eps -> 0. m_T, m_R sampled on a u-grid."""
import sys, numpy as np
from multiprocessing import Pool
import amoeba as A, optc_defs as O
fl = dict(O.fl); fl['eta'] = np.exp(-0.9j)
def edges(args):
    n, frac = args; eps = 4.0/n; x = np.exp(-eps/2); mB = int(round(frac*n))
    uT = np.arange(-0.5, 1.51, 0.5); uR = np.arange(-4.0, 3.01, 1.0); uL = np.arange(-6.0, 3.51, 0.5)
    ok = []
    for ul in uL:
        mL = int(round(ul/eps)); hit = False
        for ut in uT:
            for ur in uR:
                if A.residual((mL, int(round(ut/eps)), mB, int(round(ur/eps))), x, n, fl, starts=4, rng=np.random.default_rng(1)) < 1e-3:
                    hit = True; break
            if hit: break
        ok.append(hit)
    ok = np.array(ok); sup = uL[ok]
    return n, frac, (sup.min(), sup.max()) if ok.any() else None
if __name__ == "__main__":
    ns = [int(a) for a in sys.argv[1].split(',')]
    with Pool() as p:
        for n, frac, e in p.map(edges, [(n, f) for n in ns for f in (0.0, -1.0)]):
            print(f"n={n} eps={4/n:.4f} m_B/n={frac}: predicted u_L support {e}", flush=True)
