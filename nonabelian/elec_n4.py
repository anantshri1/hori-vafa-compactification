import numpy as np, fluxmap_n as F
y3 = np.exp(1j*np.array([0.4, 1.1, -0.6])); y = list(y3) + [1/np.prod(y3)]; e = np.exp(0.9j)
for n in (2, 4):
    E = F.el(y, e, n, 6); pre = (F.x**(3*(1-F.r))*e*(y[0]*y[1])**2)**(n/2)
    open(f'elec_n{n}.txt', 'w').write(f"{E.real!r} {E.imag!r} {pre.real!r} {pre.imag!r}\n"); print(n, E, pre, flush=True)
