"""U(2) with [4,0] (thesis U(2)_{(0,-2)} per eq. u_n_chiral at N=2, F=4)  vs  its planar abelian mirror (eq. mirpartfun),
in the index, using the dictionary fixed by cal_N1F2c.py:  mirror chirals R=(1-r)/2 (nominal 0), adjoint R=1+r (nominal 2),
mixed gauge-eta term kept, eta-fugacity = w^{-1}, FI fugacities = ratios of y's (sign scanned).
The thesis formula's right-triangle bifundamental indices are ambiguous (z^{(3)}_2 does not exist): variants A and B scanned.
Results are appended to u2_index.out (checkpoint: finished variants are skipped on restart)."""
import numpy as np, json, os, sys
from idx import theory_index
x, r = float(sys.argv[1]) if len(sys.argv) > 1 else 0.15, 0.3
h = (1-r)/2; Radj = 1 + r
M_el, M_mir, Ng_el, Ng_mir = 8, 3, 32, 12
pts = [((np.exp(0.4j), np.exp(1.1j), np.exp(-0.6j)), np.exp(0.9j)), ((np.exp(-0.8j), np.exp(0.5j), np.exp(1.7j)), np.exp(-0.7j))]
def ys(t): y1, y2, y3 = t; return [y1, y2, y3, 1/(y1*y2*y3)]
def electric(yv, w):
    y = ys(yv)
    ch = [((1, 0), {f'y{j}': 1}, r) for j in range(4)] + [((0, 1), {f'y{j}': 1}, r) for j in range(4)]
    return theory_index(2, ch, [[-1, -1], [-1, -1]], {}, {0: ['w'], 1: ['w']}, {**{f'y{j}': y[j] for j in range(4)}, 'w': w},
                        x, M=M_el, Ngrid=Ng_el, vector=2, weyl=2)
# mirror variables: 0=z1, 1=z21, 2=z22, 3=z3
def mirror(yv, w, variant, s, s1, s3):
    y = ys(yv)
    K = np.eye(4).tolist()
    def mix(a, b, v): K[a][b] += v; K[b][a] += v
    mix(1, 2, 0.5)        # adjoint pairing e^{-i pi z21 z22}
    mix(0, 2, -0.5)       # left bifundamental e^{+i pi z1 z22}
    ch = [((0, -1, 1, 0), {}, Radj),
          ((1, -1, 0, 0), {}, h), ((-1, 0, 1, 0), {}, h)]
    if variant == 'A':
        ch += [((0, -1, 0, 1), {}, h), ((0, 0, 1, -1), {}, h)]; mix(3, 2, -0.5)
    else:
        ch += [((0, 0, -1, 1), {}, h), ((0, 1, 0, -1), {}, h)]; mix(3, 1, -0.5)
    # links missing a mixed CS term in the printed formula (parity requires K_ab = 1/2 sum Q_a Q_b mod 1): signs scanned
    mix(0, 1, 0.5*s1)
    if variant == 'A': mix(1, 3, 0.5*s3)
    else: mix(2, 3, 0.5*s3)
    ch += [((0, 1, 0, 0), {'eta': -1}, h), ((0, 0, -1, 0), {}, h)]
    fl = {'eta': 1/w, 'f1': (y[1]/y[0])**s, 'f2': (y[2]/y[1])**s, 'f3': (y[3]/y[2])**s}
    return theory_index(4, ch, K, {(1, 'eta'): -0.5}, {0: ['f1'], 1: ['f2'], 2: ['f2'], 3: ['f3']}, fl, x, M=M_mir, Ngrid=Ng_mir)
out = f'u2_index_x{x}.out'
done = set(open(out).read().split('\n')) if os.path.exists(out) else set()
E = [electric(yv, w) for yv, w in pts]
print('electric', E, flush=True)
import itertools
for variant, s, s1, s3 in itertools.product('AB', (1, -1), (1, -1), (1, -1)):
        tag = f'{variant}{s}{s1}{s3}'
        if any(l.startswith(tag + ' ') for l in done): continue
        Mv = [mirror(yv, w, variant, s, s1, s3) for yv, w in pts]
        rat = [a/b for a, b in zip(E, Mv)]
        line = f'{tag} ratios {rat} spread {abs(rat[1]/rat[0]-1):.2e}'
        print(line, flush=True); open(out, 'a').write(line + '\n')
