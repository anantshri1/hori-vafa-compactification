"""LOG 18: mirror index with eta-flux -n, flux box centred at m = (0,0,0,0) vs centred at m_B = -n (B chiral light).
Ratio of the two boxes and their eta dependence (relative to the zero-point prefactor (eta x^{1-r})^{n/2})."""
import numpy as np, sys
from idx import theory_index
import optc_defs as O
r, x = O.r, float(sys.argv[2]) if len(sys.argv) > 2 else 0.05
n = int(sys.argv[1]); M = 2
full_ch = [(O.one(O.B, 1), {'eta': -1}, r)] + O.rest
for th in (0.9, -0.5):
    fl = dict(O.fl); fl['eta'] = np.exp(-1j*th)
    args = (4, full_ch, O.Kmat(1), {(O.B, 'eta'): -.5}, O.FIfull, fl, x)
    C = (fl['eta']*x**(1 - r))**(n/2)
    I0 = theory_index(*args, M=M, Ngrid=20, bg={'eta': -n})/C
    for sh in ([0, 0, -n, 0], [-n//2, -n//2, -n, -n//2], [0, -n//2, -n, 0]):
        Is = theory_index(*args, M=M, Ngrid=20, bg={'eta': -n}, mshift=sh)/C
        print(f"x={x} n={n} theta={th}: box0/C={I0:.4e}  shift {sh}: /C={Is:.4e}  |shift/box0|={abs(Is/I0):.2e}", flush=True)
