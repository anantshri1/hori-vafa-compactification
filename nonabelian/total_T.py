"""LOG 53 (C)-(D): mirror total T at eps = 0.2, eta-flux -64 (window_mirror.py fugacities) via the electric side and the
now-verified map I_el(n) = (x^{3(1-r)} w (y1y2)^2)^{n/2} I_mir(-n)."""
import numpy as np, elec_fourier as EF, time
eps = 0.2; x = np.exp(-eps/2); r = 0.3; ns = int(round(2*np.log(1/eps)/eps)); n = 4*ns
y3 = np.exp(1j*np.array([0.4, 1.1, -0.6])); y = list(y3) + [1/np.prod(y3)]; e = np.exp(0.9j)
pre = (x**(3*(1-r))*e*(y[0]*y[1])**2)**(n/2)
for M in (20, 40, 60):
    t0 = time.time(); E = EF.index(y, e, n, x, r, M, Nf=1024)
    print(f"eps={eps} n={n} M={M}: I_el = {E:.10e}  ->  T_mirror = {E/pre:.10e}  [{time.time()-t0:.0f}s]", flush=True)
