"""LOG 18: profile of the mirror index (eta-flux -n) along the B-flux line: slice m_B fixed, others |m|<=1 (M=1 box)."""
import numpy as np, sys
from idx import theory_index
import optc_defs as O
n = int(sys.argv[1]); x = float(sys.argv[2]); r = O.r
fl = dict(O.fl); fl['eta'] = np.exp(-0.9j); C = (fl['eta']*x**(1 - r))**(n/2)
ch = [(O.one(O.B, 1), {'eta': -1}, r)] + O.rest
for mB in range(2, -n - 3, -1):
    I = theory_index(4, ch, O.Kmat(1), {(O.B, 'eta'): -.5}, O.FIfull, fl, x, M=1, Ngrid=16, bg={'eta': -n},
                     signfn=lambda m, mB=mB: 1 if m[O.B] == mB else 0, mshift=[0, 0, mB, 0]) / C
    print(f"n={n} x={x} m_B={mB:+d} (m_B/n={mB/n:+.2f}): slice/C = {I:.4e}  |.|={abs(I):.3e}", flush=True)
