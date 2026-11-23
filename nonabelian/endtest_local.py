"""Quick local look: smooth-windowed (Gaussian, sigma_m = 3) sums at eps = 0.3 around the flip maxima: exact vs ladder LG
vs end-corrected LG (|nu| Bessel ends) vs ends exact."""
import numpy as np, itertools
from multiprocessing import Pool
import fourier_sector as FS
FS.set_r(0.3); x = np.exp(-0.15); n = 32
y3 = np.exp(1j*np.array([0.4, 1.1, -0.6])); y = list(y3) + [1/np.prod(y3)]; e = np.exp(0.9j); rat = [y[i+1]/y[i] for i in range(3)]
fl = {'eta': 1/e, 'f1': rat[0], 'fT': rat[1]*x**0.15, 'fB': rat[1]*x**-0.15, 'f3': rat[2]}
def job(m): return FS.sector(m, x, fl, n, 12, 128), FS.sector_lg(m, x, fl, n, 12, 128), FS.sector_lg_end(m, x, fl, n, 12, 128), FS.sector_lg_end(m, x, fl, n, 12, 128, True)
if __name__ == "__main__":
    for name, c in {'ladder': (-16, -8, -24, -16), 'Bflip': (-18, -2, -34, -18), 'Tflip': (-14, 2, -30, -14), 'both': (-14, 6, -36, -14)}.items():
        D = list(itertools.product(range(-5, 6), repeat=4)); g = np.array([np.exp(-np.dot(d, d)/18) for d in D])
        with Pool() as p: R = np.array(p.map(job, [tuple(np.array(c) + np.array(d)) for d in D], chunksize=64))
        s = (g[:, None]*R).sum(0)
        print(f"{name:6s} exact {s[0]:.4e} | ladderLG {s[1]:.4e} ({abs(s[1]/s[0]-1):.2f}) | end-Bessel {s[2]:.4e} ({abs(s[2]/s[0]-1):.2f}) | ends exact {s[3]:.4e} ({abs(s[3]/s[0]-1):.2f})", flush=True)
