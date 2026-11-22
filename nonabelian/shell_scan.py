"""LOG 55: where is the missing ~20%? Coarse |S| (stride 4) over a bigger box, outside the LOG 54 box, by face."""
import numpy as np, itertools, time
from multiprocessing import Pool
import fourier_sector as FS
eps, r = 0.3, 0.3; x = np.exp(-eps/2); n = 32
FS.set_r(r)
y3 = np.exp(1j*np.array([0.4, 1.1, -0.6])); y = list(y3) + [1/np.prod(y3)]; e = np.exp(0.9j); rat = [y[i+1]/y[i] for i in range(3)]
fl = {'eta': 1/e, 'f1': rat[0], 'fT': rat[1]*x**(r/2), 'fB': rat[1]*x**(-r/2), 'f3': rat[2]}
old = [(-40, 8), (-12, 12), (-40, 8), (-40, 8)]        # (L,T,B,R)
rng = [range(-72, 33, 6), range(-36, 37, 6), range(-72, 33, 6), range(-72, 33, 6)]
pts = [m for m in itertools.product(*rng) if not all(old[a][0] <= m[a] <= old[a][1] for a in range(4))]
inside = [m for m in itertools.product(*rng) if all(old[a][0] <= m[a] <= old[a][1] for a in range(4))]
def job(m): return m, abs(FS.sector(m, x, fl, n, Kc=12, Nf=128))
if __name__ == "__main__":
    t0 = time.time()
    with Pool() as p: out = p.map(job, pts + inside, chunksize=64)
    nm = 'LTBR'; faces = {}; tin = 0
    for m, a in out:
        if all(old[k][0] <= m[k] <= old[k][1] for k in range(4)): tin += a; continue
        key = ' '.join(f"{nm[k]}{'<' if m[k] < old[k][0] else '>'}" for k in range(4) if not (old[k][0] <= m[k] <= old[k][1]))
        faces[key] = faces.get(key, 0) + a
    tout = sum(faces.values())
    print(f"{len(pts)} outside + {len(inside)} inside sampled [{time.time()-t0:.0f}s]; sum|S| inside {tin:.3e}, outside {tout:.3e} ({tout/(tin+tout):.3f})")
    for k, v in sorted(faces.items(), key=lambda t: -t[1])[:10]: print(f"  outside via {k:12s}: {v/(tin+tout):.4f}")
