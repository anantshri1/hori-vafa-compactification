"""Engine check: reproduce the already-verified abelian CS pair (symmetric frame, pair_map.py: wB = y^2, a = wA^{1/2},
rB = (1-rA)/2) with idx.theory_index. Must agree to ~1e-12 before the engine is used for U(2)."""
import numpy as np
from idx import theory_index
x, r = 0.3, 0.3
for y, w in [(np.exp(0.4j), np.exp(0.9j)), (np.exp(1.3j), np.exp(-0.7j))]:
    A = theory_index(1, [((1,), {'y': 1}, r), ((1,), {'y': -1}, r)], [[-1]], {}, {0: ['w']}, {'y': y, 'w': w}, x, M=14, Ngrid=64)
    B = theory_index(1, [((1,), {'a': 1}, (1-r)/2), ((-1,), {'a': 1}, (1-r)/2)], [[1]], {}, {0: ['wb']},
                     {'a': np.sqrt(w), 'wb': y**2}, x, M=14, Ngrid=64)
    print(A, B, abs(A/B-1))
