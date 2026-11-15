"""Background-flux check of the engine on the verified abelian CS pair: I_A(n) = (w x^{3r-1})^{n/4} I_B(n_a = n/2)
(pair_fluxmap.py, symmetric frame). Scans the sign of the B-side axial flux."""
import numpy as np
from idx import theory_index
x, r = 0.3, 0.3
y, w = np.exp(0.4j), np.exp(0.9j)
for n in (2, -2, 4):
    A = theory_index(1, [((1,), {'y': 1}, r), ((1,), {'y': -1}, r)], [[-1]], {}, {0: ['w']}, {'y': y, 'w': w}, x, M=16, Ngrid=64, bg={'w': n})
    for sgn in (1, -1):
        B = theory_index(1, [((1,), {'a': 1}, (1-r)/2), ((-1,), {'a': 1}, (1-r)/2)], [[1]], {}, {0: ['wb']},
                         {'a': np.sqrt(w), 'wb': y**2}, x, M=16, Ngrid=64, bg={'a': sgn*n//2})
        pred = (w*x**(3*r-1))**(n/4)
        print(n, sgn, abs(A/(pred*B)-1), A/B, pred)
