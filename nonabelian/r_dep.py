"""Fix how the column-2 FI R-mixing c depends on r: electric (FI = eta x^{2r}) vs mirror at r = 0.3, unrefined flavours."""
import numpy as np, sys
import mirror_table as MT
from idx import theory_index
r = 0.3; MT.r = r; x = 0.05; e = np.exp(0.7j)
ch = [((1, 0), {}, 1-r)]*4 + [((0, 1), {}, 1-r)]*4
E = theory_index(2, ch, [[-1, -1], [-1, -1]], {}, {0: ['w'], 1: ['w']}, {'w': e*x**(2*r)}, x, M=8, Ngrid=40, vector=2, weyl=2)
print('E', E, flush=True)
for c in [float(a) for a in sys.argv[1:]]:
    m = MT.mirror(x, 1/e, 1, c, 0, M=2, Ng=20)
    print(f"c={c}: |E-M| = {abs(E-m):.2e}", flush=True)
