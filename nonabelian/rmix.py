"""LOG 33: R-mixing through the limit. Exact non-Pochhammer modulus of a mirror sector:
   log|prefactor| = ex(m) log x + sum_f ef_f(m) log|fl_f|  (zero points x^{(1-R)|nu|/2}, FI fugacities incl. x^{+-r/2}).
In the ladder sign region (B, L-B, R-B, T-L, T-R, -T with nu > 0; vertical B-T with nu < 0), its r-dependence must be
m-independent iff the R-mixing cancels (mu_B = mu_T exactly)."""
import numpy as np
from fractions import Fraction as Fr
import optc_defs as O
def logpref(m, r, x, n):
    chir = [((0,0,1,0), {'eta': -1}, r), ((1,0,-1,0), {}, r), ((-1,1,0,0), {}, r), ((0,-1,0,0), {}, r),
            ((0,0,-1,1), {}, r), ((0,1,0,-1), {}, r), ((0,-1,1,0), {}, 2-2*r)]
    fl = {'eta': 1.0, 'f1': 1.0, 'fT': x**(r/2), 'fB': x**(-r/2), 'f3': 1.0}     # moduli only
    FI = {0: ['f1'], 1: ['fT'], 2: ['fB'], 3: ['f3']}
    tot = 0.0; signs = []
    for Q, flc, R in chir:
        nu = sum(Q[a]*m[a] for a in range(4)) + sum(c*(-n) for f, c in flc.items())
        tot += (1 - R)*abs(nu)/2*np.log(x); signs.append(np.sign(nu))
    for a, fs in FI.items():
        for f in fs: tot += -m[a]*np.log(fl[f])
    return tot, signs
x, n = 0.9, 40
# a point in the ladder region: k-vector (L,T,B,R) = (-2,-1,-3,-2), background b = 4 (LOG 30), scaled
s = 6; m0 = np.array([-2, -1, -3, -2])*s; nb = 4*s
for r in (0.2, 0.3, 0.45):
    base, sg = logpref(m0, r, x, nb)  # eta flux sign: nu_B = m_B + 4s > 0
    grads = []
    for a in range(4):
        e = np.zeros(4, int); e[a] = 1
        grads.append(logpref(m0 + e, r, x, nb)[0] - base)
    print(f"r={r}: signs {sg}; d log|pref| / d m_a = " + " ".join(f"{g:+.6f}" for g in grads), flush=True)
