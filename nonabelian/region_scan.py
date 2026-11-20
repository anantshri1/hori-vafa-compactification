"""LOG 47: weight map of the exact mirror index in the Bessel scaling, by branch pattern (signs of nu_i).
Fluxes m = nu* (c + k/S), c = ladder centre (-2,-1,-3,-2), k in [-2S, 2S]^4 (u-box +-2 nu* around the ladder point),
eta-flux 4 nu*. Aggregates sum|S| per sign pattern. Usage: python region_scan.py eps S"""
import sys, numpy as np, itertools, time
from multiprocessing import Pool
import optc_defs as O, fourier_sector as FS
eps, S = float(sys.argv[1]), int(sys.argv[2]); Bx = float(sys.argv[3]) if len(sys.argv) > 3 else 2.0
x = np.exp(-eps/2); ns = int(round(2*np.log(1/eps)/eps)); n = 4*ns
fl = dict(O.fl); fl['eta'] = np.exp(-0.9j); fl['fT'] = fl['fT']/O.x**(O.r/2)*x**(O.r/2); fl['fB'] = fl['fB']/O.x**(-O.r/2)*x**(-O.r/2)
c = np.array([-2, -1, -3, -2])
ks = list(itertools.product(range(-int(Bx*S), int(Bx*S) + 1), repeat=4))
def job(k):
    m = tuple(int(round(v)) for v in ns*(c + np.array(k)/S))
    v = FS.sector(m, x, fl, n, Kc=16, Nf=256)
    _, nus, _, _ = FS.bookkeeping(m, x, fl, n)
    return ''.join('+' if s > 0 else ('-' if s < 0 else '0') for s in np.sign(nus)), abs(v), m
if __name__ == "__main__":
    t0 = time.time(); agg = {}; best = {}
    with Pool() as p:
        for pat, a, m in p.imap_unordered(job, ks, chunksize=64):
            agg[pat] = agg.get(pat, 0) + a
            if a > best.get(pat, (0,))[0]: best[pat] = (a, m)
    tot = sum(agg.values())
    print(f"eps={eps} nu*={ns} S={S}: {len(ks)} sectors in {time.time()-t0:.0f}s; total sum|S| = {tot:.4e}")
    for pat, a in sorted(agg.items(), key=lambda t: -t[1])[:12]:
        print(f"  pattern {pat} (B,L-B,T-L,-T,R-B,T-R,B-T): share {a/tot:.4f}  max sector {best[pat][0]:.3e} at m={best[pat][1]}")
