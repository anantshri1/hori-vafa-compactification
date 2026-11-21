import numpy as np, elec_fourier as EF, gateA2 as G, itertools
x, r = 0.03, 0.45
for yang, wang in [([0.4, 1.1, -0.6], 0.9), ([-0.8, 0.5, 1.7], -0.7)]:
    y, e, fl = G.fl_of(yang, wang)
    for n in (2, 3, 4, 6):
        G.n = n
        Mt = G.mirror_total(fl, [range(-4, 5), range(-4, 5), range(-n - 4, 5), range(-4, 5)])
        E = EF.index(y, e, n, x, r, 10); pre = (x**(3*(1-r))*e*(y[0]*y[1])**2)**(n/2)
        print(f"point {yang}: n={n}: C = {E/Mt/pre:.8f}", flush=True)
