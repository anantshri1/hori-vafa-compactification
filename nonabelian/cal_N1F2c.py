"""Calibration step 3: thesis-frame mirror (Q: +1 with eta, Qt: -1) with free R-charges, eta and FI x-powers, and the mixed
gauge-eta CS term on/off. Electric side fixed: R=r, omega=w (verified symmetric-frame electric)."""
import numpy as np, itertools
from idx import theory_index
x, r = 0.3, 0.3; h = (1-r)/2
pts = [(np.exp(0.4j), np.exp(0.9j)), (np.exp(1.3j), np.exp(-0.7j)), (np.exp(2.2j), np.exp(1.6j))]
E = [theory_index(1, [((1,), {'y': 1}, r), ((1,), {'y': -1}, r)], [[-1]], {}, {0: ['om']}, {'y': y, 'om': w}, x, M=10, Ngrid=48) for y, w in pts]
res = []
for kf, t, RQ, RQt, u in itertools.product([0, -0.5], [0, -h, -2*h], [0, h, 2*h], [0, h, 2*h], [-2*h, -h, 0, h, 2*h]):
    Kf = {(0, 'eta'): kf} if kf else {}
    Mv = [theory_index(1, [((1,), {'eta': -1}, RQ), ((-1,), {}, RQt)], [[1]], Kf, {0: ['fi']},
                       {'eta': w**-1*x**t, 'fi': y**2*x**u}, x, M=10, Ngrid=48) for y, w in pts]
    rat = [a/b for a, b in zip(E, Mv)]
    res.append((max(abs(z/rat[0]-1) for z in rat), kf, t, RQ, RQt, u, rat[0]))
res.sort(key=lambda z: z[0])
for z in res[:6]: print(f"spread {z[0]:.1e}: Kf={z[1]}, eta=w^-1 x^{z[2]:.2f}, R=({z[3]:.2f},{z[4]:.2f}), FI=y^2 x^{z[5]:.2f}, ratio {z[6]:.6f}")
