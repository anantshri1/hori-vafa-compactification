"""Option (c): integrate out the B chiral made massive by eta-flux -n. Prediction (index zero-point algebra):
drop B chiral and the B-eta mixed CS; K_BB: 1 -> 1/2; extra FI fugacity on B: x^{-(1-r)/2}; overall C = (eta x^{1-r})^{n/2}."""
import numpy as np
from idx import theory_index
r = 0.3; x = 0.05
L, T, B, R = 0, 1, 2, 3
def f(h, t): Q = [0]*4; Q[h] += 1; Q[t] -= 1; return tuple(Q)
def one(a, s): Q = [0]*4; Q[a] = s; return tuple(Q)
def Kmat(kBB):
    K = np.eye(4).tolist(); K[B][B] = kBB
    for a, b, v in [(L, B, -.5), (L, T, -.5), (R, B, -.5), (R, T, -.5), (T, B, .5)]: K[a][b] += v; K[b][a] += v
    return K
rest = [(f(L, B), {}, r), (f(T, L), {}, r), (one(T, -1), {}, r), (f(R, B), {}, r), (f(T, R), {}, r), (f(B, T), {}, 2 - 2*r)]
y = [np.exp(0.4j), np.exp(1.1j), np.exp(-0.6j)]; y.append(1/np.prod(y)); e = np.exp(0.9j)
rat = [y[i+1]/y[i] for i in range(3)]
fl = {'eta': 1/e, 'f1': rat[0], 'fT': rat[1]*x**(r/2), 'fB': rat[1]*x**(-r/2), 'f3': rat[2], 'g': x**(-(1-r)/2)}
FIfull = {L: ['f1'], T: ['fT'], B: ['fB'], R: ['f3']}
