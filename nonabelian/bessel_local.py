"""Local look (eps = 0.3, Gaussian windows sigma_m = 3, 11^4 sectors): branch model vs Bessel-switch model."""
import numpy as np, itertools, sys
from multiprocessing import Pool
import fourier_sector as FS
FS.set_r(0.3); x = np.exp(-0.15); n = 32
y3 = np.exp(1j*np.array([0.4, 1.1, -0.6])); y = list(y3) + [1/np.prod(y3)]; e = np.exp(0.9j); rat = [y[i+1]/y[i] for i in range(3)]
fl = {'eta': 1/e, 'f1': rat[0], 'fT': rat[1]*x**0.15, 'fB': rat[1]*x**-0.15, 'f3': rat[2]}
def job(m): return FS.sector(m, x, fl, n, 12, 128), FS.sector_branch(m, x, fl, n, 12, 128), FS.sector_branch(m, x, fl, n, 12, 128, True), FS.sector_bessel(m, x, fl, n, 12, 128, 1.0)
if __name__ == "__main__":
    for name, c in {'ladder': (-16, -8, -24, -16), 'Bflip': (-18, -2, -34, -18), 'Tflip': (-14, 2, -30, -14), 'both': (-14, 6, -36, -14), 'diagflip': (-30, -8, -24, -16)}.items():
        D = list(itertools.product(range(-4, 5), repeat=4)); g = np.array([np.exp(-np.dot(d, d)/12) for d in D])
        with Pool() as p: R = np.array(p.map(job, [tuple(np.array(c) + np.array(d)) for d in D], chunksize=64))
        s = (g[:, None]*R).sum(0)
        print(f"{name:8s} exact {s[0]:.4e} | branch {abs(s[1]/s[0]-1):.3f} | vertical exact {abs(s[2]/s[0]-1):.3f} | switch nu* {abs(s[3]/s[0]-1):.3f}", flush=True)
