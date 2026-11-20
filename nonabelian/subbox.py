"""LOG 48: COHERENT test. Exact contiguous sums of mirror sectors over sub-boxes (stride 1) centred in each region,
vs the same sums with the ladder-LG integrand (sector_lg). Usage: python subbox.py eps half_width"""
import sys, numpy as np, itertools, time
from multiprocessing import Pool
import optc_defs as O, fourier_sector as FS
eps, hw = float(sys.argv[1]), int(sys.argv[2]); only = sys.argv[3].split(',') if len(sys.argv) > 3 else None
x = np.exp(-eps/2); ns = int(round(2*np.log(1/eps)/eps)); n = 4*ns
fl = dict(O.fl); fl['eta'] = np.exp(-0.9j); fl['fT'] = fl['fT']/O.x**(O.r/2)*x**(O.r/2); fl['fB'] = fl['fB']/O.x**(-O.r/2)*x**(-O.r/2)
# region centres in units of nu* (from LOG 47 maxima at eps=0.2: ladder (-32,-20,-44,-32), -T flip (-28,4,-60,-28), B flip (-36,-4,-68,-36))
centres = {'ladder ++++++-': (-2, -1.25, -2.75, -2), '-T flip +++-++-': (-1.75, 0.25, -3.75, -1.75),
           'B flip -+++++-': (-2.25, -0.25, -4.25, -2.25), 'far (+3nu* off)': (1, 2, 0, 1)}
def job(args):
    m = args; return FS.sector(m, x, fl, n), FS.sector_lg(m, x, fl, n)
if __name__ == "__main__":
    for name, c in centres.items():
        if only and not any(o in name for o in only): continue
        c0 = np.round(np.array(c)*ns).astype(int)
        secs = [tuple(int(v) for v in c0 + np.array(d)) for d in itertools.product(range(-hw, hw + 1), repeat=4)]
        t0 = time.time()
        with Pool() as p: res = p.map(job, secs, chunksize=64)
        ex = sum(r[0] for r in res); lg = sum(r[1] for r in res); ab = sum(abs(r[0]) for r in res)
        print(f"eps={eps} hw={hw} {name:18s} centre {tuple(c0)}: exact sum {ex:.6e} (sum|.| {ab:.3e}, coherence {abs(ex)/ab:.2e}) "
              f"| LG sum {lg:.6e} | |exact-LG|/sum|.| {abs(ex-lg)/ab:.2e}  [{time.time()-t0:.0f}s]", flush=True)
