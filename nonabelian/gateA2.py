"""Gate A' (LOG 52): n = 4 flux map with the mirror from the Fourier engine over a box containing m_B in [-n-3, 3]."""
import numpy as np, itertools, sys, time
from multiprocessing import Pool
import fourier_sector as FS
x, r, n = 0.03, 0.45, int(sys.argv[1]) if len(sys.argv) > 1 else 4
FS.set_r(r)
def fl_of(yang, wang):
    y3 = np.exp(1j*np.array(yang)); y = list(y3) + [1/np.prod(y3)]; e = np.exp(1j*wang)
    rat = [y[i+1]/y[i] for i in range(3)]
    return y, e, {'eta': 1/e, 'f1': rat[0], 'fT': rat[1]*x**(r/2), 'fB': rat[1]*x**(-r/2), 'f3': rat[2]}
def mirror_total(fl, box):
    secs = list(itertools.product(*box))
    with Pool() as p: vals = p.starmap(FS.sector, [(s, x, fl, n, 12, 128) for s in secs], chunksize=64)
    return sum(vals)
if __name__ == "__main__":
    for yang, wang in [([0.4, 1.1, -0.6], 0.9)]:
        y, e, fl = fl_of(yang, wang)
        for h in (4, 5):
            box = [range(-h, h + 1), range(-h, h + 1), range(-n - h, h + 1), range(-h, h + 1)]
            t0 = time.time(); Mt = mirror_total(fl, box)
            print(f"n={n} mirror (Fourier, box h={h}, m_B in [{-n-h},{h}]) = {Mt:.10e}  [{time.time()-t0:.0f}s]", flush=True)
