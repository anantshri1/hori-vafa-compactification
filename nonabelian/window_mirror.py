"""LOG 49: Gaussian-windowed flux sums of the planar mirror (exact vs symmetric ladder LG), eps fixed.
Window g(m) = exp(-|eps(m - m_c)|^2 / (2 s^2)) in u = eps m, summed over |m - m_c| <= 3 s/eps (stride 1).
Usage: python window_mirror.py eps s region[,region]"""
import sys, numpy as np, itertools, time
from multiprocessing import Pool
import optc_defs as O, fourier_sector as FS
eps, s = float(sys.argv[1]), float(sys.argv[2]); only = sys.argv[3].split(',')
x = np.exp(-eps/2); ns = int(round(2*np.log(1/eps)/eps)); n = 4*ns
fl = dict(O.fl); fl['eta'] = np.exp(-0.9j); fl['fT'] = fl['fT']/O.x**(O.r/2)*x**(O.r/2); fl['fB'] = fl['fB']/O.x**(-O.r/2)*x**(-O.r/2)
centres = {'ladder': (-2, -1.25, -2.75, -2), 'Tflip': (-1.75, 0.25, -3.75, -1.75), 'Bflip': (-2.25, -0.25, -4.25, -2.25), 'far': (1, 2, 0, 1)}
def job(m): return FS.sector(m, x, fl, n), FS.sector_lg(m, x, fl, n)
if __name__ == "__main__":
    h = int(np.ceil(3*s/eps))
    for name in only:
        c0 = np.round(np.array(centres[name])*ns).astype(int)
        D = list(itertools.product(range(-h, h + 1), repeat=4))
        secs = [tuple(int(v) for v in c0 + np.array(d)) for d in D]
        g = np.array([np.exp(-eps**2*np.dot(d, d)/(2*s*s)) for d in D])
        t0 = time.time()
        with Pool() as p: res = np.array(p.map(job, secs, chunksize=64))
        ex, lg = np.sum(g*res[:, 0]), np.sum(g*res[:, 1]); ab = np.sum(g*np.abs(res[:, 0]))
        print(f"eps={eps} s={s} window {name:6s} centre {tuple(int(v) for v in c0)} ({len(D)} sectors): exact {ex:.6e} | LG {lg:.6e} | "
              f"|exact/LG - 1| {abs(ex/lg - 1) if abs(lg) > 0 else np.inf:.2e} | coherence {abs(ex)/ab:.2e}  [{time.time()-t0:.0f}s]", flush=True)
