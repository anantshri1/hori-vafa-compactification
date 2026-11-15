"""Which side moved? electric (M,Ng) = (8,40),(10,56),(10,40); mirror Ng = 20, 26 at point 1, flux 2 / eta-flux -2."""
import numpy as np, bg_u2_grid as B
from idx import theory_index
x, r = B.x, B.r
y3 = np.exp(1j*np.array([0.4, 1.1, -0.6])); y = list(y3) + [1/np.prod(y3)]; e = np.exp(0.9j)
ch = [((1, 0), {f'y{j}': 1}, 1-r) for j in range(4)] + [((0, 1), {f'y{j}': 1}, 1-r) for j in range(4)]
for M, Ng in [(8, 40), (10, 40), (8, 56), (10, 56)]:
    E = theory_index(2, ch, [[-1, -1], [-1, -1]], {}, {0: ['w'], 1: ['w']}, {**{f'y{j}': y[j] for j in range(4)}, 'w': e*x**(2*r)},
                     x, M=M, Ngrid=Ng, vector=2, weyl=2, bg={'w': 2})
    print('E', M, Ng, E, flush=True)
for Ng in (20, 26):
    print('Mir', Ng, B.mir(y, e, -2, Ng), flush=True)
