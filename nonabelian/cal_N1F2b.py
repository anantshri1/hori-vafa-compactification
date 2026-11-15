"""Calibration step 2: add an x-power on the mirror FI fugacity (the image of a gauge-R mixing / imaginary shift of z),
scan u in multiples of (1-r)/2 together with the eta-fugacity power t."""
import numpy as np
from idx import theory_index
x, r = 0.3, 0.3
pts = [(np.exp(0.4j), np.exp(0.9j)), (np.exp(1.3j), np.exp(-0.7j)), (np.exp(2.2j), np.exp(1.6j))]
h = (1-r)/2
def el(y, w, p):
    return theory_index(1, [((1,), {'y': 1}, r), ((1,), {'y': -1}, r)], [[-1]], {}, {0: ['om']},
                        {'y': y, 'om': w*x**p}, x, M=14, Ngrid=64)
def mir(y, w, t, u, s):
    return theory_index(1, [((1,), {'eta': -1}, 0), ((-1,), {}, 0)], [[1]], {(0, 'eta'): -0.5}, {0: ['fi']},
                        {'eta': w**-1*x**t, 'fi': y**(2*s)*x**u}, x, M=14, Ngrid=64)
res = []
for p in [0, 2*h]:
    E = [el(y, w, p) for y, w in pts]
    for t in [k*h for k in range(-4, 5)]:
        for u in [k*h for k in range(-4, 5)]:
            Mv = [mir(y, w, t, u, 1) for y, w in pts]
            rat = [a/b for a, b in zip(E, Mv)]
            res.append((max(abs(z/rat[0]-1) for z in rat), p, t, u, rat[0]))
res.sort(key=lambda z: z[0])
for z in res[:5]: print(f"spread {z[0]:.1e}: el omega=w x^{z[1]:.2f}; mir eta-fug=w^-1 x^{z[2]:.2f}, FI=y^2 x^{z[3]:.2f}; ratio {z[4]:.6f}")
