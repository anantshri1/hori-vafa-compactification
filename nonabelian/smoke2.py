"""SMOKE 2: Bessel-regime sectors: cutoff convergence (Kc, Nf) and timing, eps = 0.2, 0.1."""
import numpy as np, time, optc_defs as O, fourier_sector as FS
fl = dict(O.fl); fl['eta'] = np.exp(-0.9j)
for eps in (0.2, 0.1):
    x = np.exp(-eps/2); ns = int(round(2*np.log(1/eps)/eps)); n = 4*ns
    fl2 = dict(fl); fl2['fT'] = fl['fT']/O.x**(O.r/2)*x**(O.r/2); fl2['fB'] = fl['fB']/O.x**(-O.r/2)*x**(-O.r/2)   # R-mixing at this x
    m = tuple(int(v) for v in np.array([-2, -1, -3, -2])*ns + np.array([1, 0, -1, 2]))
    for Kc, Nf in [(16, 256), (24, 256), (32, 512)]:
        t0 = time.time(); v = FS.sector(m, x, fl2, n, Kc=Kc, Nf=Nf)
        print(f"eps={eps} nu*={ns} m={m}: Kc={Kc} Nf={Nf}: {v:.10e}  ({time.time()-t0:.2f}s)", flush=True)
