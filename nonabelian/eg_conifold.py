"""LOG 41: local conifold structure. At Gr vacuum {12}: Gr weights {a,b,c,d} with a+d = b+c =: f (conifold relation
xy = zw); the two resolved fixed points have weights {a,c,e,f} and {b,d,-e,f}, e = b-a. phi(w) = theta1(w - z)/theta1(w)
(y = e^{2 pi i z}). D_loc(a,b,c) = phi(f)[phi(a)phi(c)phi(e) + phi(b)phi(d)phi(-e)] - phi(a)phi(b)phi(c)phi(d).
Check: D = D_loc({12}) + D_loc({34}) with weights read off from the fixed-point data."""
import numpy as np, eg as G
def phi(w, z): return G.th1(np.exp(2j*np.pi*(w - z)))/G.th1(np.exp(2j*np.pi*w))
def Dloc(a, b, c, z):
    d = b + c - a; f = b + c; e = b - a
    return phi(f, z)*(phi(a, z)*phi(c, z)*phi(e, z) + phi(b, z)*phi(d, z)*phi(-e, z)) - phi(a, z)*phi(b, z)*phi(c, z)*phi(d, z)
M = G.M
for yv in (0.13, 0.31):
    z = yv; y = np.exp(2j*np.pi*z)
    # {12} = (0,1): Gr weights M_k - M_j, j in {0,1}, k in {2,3}: a = M2-M1? identify from LOG: a+d = b+c
    W01 = [M[2]-M[0], M[3]-M[0], M[2]-M[1], M[3]-M[1]]   # (2-0),(3-0),(2-1),(3-1): (2-0)+(3-1) = (3-0)+(2-1)
    W23 = [M[0]-M[2], M[1]-M[2], M[0]-M[3], M[1]-M[3]]
    def D_for(W):  # choose a,b,c with a+d=b+c: a=W[0], d=W[3], b=W[1], c=W[2]
        return Dloc(W[0], W[1], W[2], z)
    D = G.eg_tor(G.mu, y) - G.eg_gr(M, y)
    print(f"y=e^(2pi i {yv}): D = {D:.8f};  D_loc(12)+D_loc(34) = {D_for(W01)+D_for(W23):.8f}", flush=True)
