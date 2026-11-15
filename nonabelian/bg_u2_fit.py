"""Fit the U(2) flux-sector ratio to a monomial C w^a y1^b1 y2^b2 y3^b3 (y4 = 1/(y1 y2 y3)): vary one fugacity at a time
by 0.3 rad and read off the power from the phase change; then check the prediction at an independent point."""
import numpy as np, sys
import bg_u2 as G   # reuses el, mir (module's own scan is guarded below)
n = 2; nb = int(sys.argv[1])
base_y = np.array([0.4, 1.1, -0.6]); base_w = 0.9; d = 0.3
def ratio(yang, wang):
    y3 = np.exp(1j*np.array(yang)); y = list(y3) + [1/np.prod(y3)]; e = np.exp(1j*wang)
    return G.el(y, e, n)/G.mir(y, e, nb)
R0 = ratio(base_y, base_w); print('R0', R0, flush=True)
pw = {}
Rw = ratio(base_y, base_w + d); pw['w'] = np.angle(Rw/R0)/d; print('w', pw['w'], abs(Rw/R0), flush=True)
for j in range(3):
    yy = base_y.copy(); yy[j] += d
    Rj = ratio(yy, base_w); pw[f'y{j+1}'] = np.angle(Rj/R0)/d; print(f'y{j+1}', pw[f'y{j+1}'], abs(Rj/R0), flush=True)
