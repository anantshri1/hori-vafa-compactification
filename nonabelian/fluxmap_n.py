"""Gate A (LOG 51): general-n flux map. Conjecture: I_el(top flux n) = (x^{3(1-r)} w (y1 y2)^2)^{n/2} I_mir(eta-flux -n).
Stable regime of bg_u2_stable.py (r=0.45, x=0.03). Prints C(n) = I_el/I_mir/(x^{3(1-r)} w (y1y2)^2)^{n/2}."""
import numpy as np, refined as RF, idx, sys
from idx import theory_index
orig = idx.theory_index
RF.x, RF.r = 0.03, 0.45; x, r = RF.x, RF.r
def el(y, e, n, M):
    ch = [((1, 0), {f'y{j}': 1}, 1-r) for j in range(4)] + [((0, 1), {f'y{j}': 1}, 1-r) for j in range(4)]
    return theory_index(2, ch, [[-1, -1], [-1, -1]], {}, {0: ['w'], 1: ['w']}, {**{f'y{j}': y[j] for j in range(4)}, 'w': e*x**(2*r)},
                        x, M=M, Ngrid=96, vector=2, weyl=2, bg={'w': n})
def mir(y, e, n, M):
    def patched(*a, **k): k['bg'] = {'eta': -n}; return orig(*a, **k)
    RF.theory_index = patched
    v = RF.mirror(y, e, 1, M=M, Ng=24); RF.theory_index = orig; return v
for n in [int(a) for a in sys.argv[1:]] or [2, 4]:
    for yang, wang in [([0.4, 1.1, -0.6], 0.9), ([-0.8, 0.5, 1.7], -0.7)]:
        y3 = np.exp(1j*np.array(yang)); y = list(y3) + [1/np.prod(y3)]; e = np.exp(1j*wang)
        E = el(y, e, n, 6); pre = (x**(3*(1-r))*e*(y[0]*y[1])**2)**(n/2)
        print(f"n={n} point {yang}: C = {E/mir(y, e, n, 3)/pre:.6f}  (mirror M=2: {E/mir(y, e, n, 2)/pre:.6f})", flush=True)
