"""
Step 1a: exact 3d identity  U(1)_{1/2} + one chiral (charge 1, R-charge R)  <->  free chiral,
in Anant's convention, with k = 1/2 taken literally in Z_cl = (-1)^{k m} w^{-m} z^{-k m - n}.

Gauge charge of the bare monopole in sector m: z^{-|m|/2} (chiral zero-point) * z^{-m/2} (CS)
= z^{-(m+|m|)/2}  -> m < 0 sectors are neutral, so T_- (m = -1) is the candidate dual chiral,
with R-charge (1-R)/2 (from x^{(1-R)|m|/2}) and topological fugacity w (from w^{-m}).

The sign (-1)^{k m} = (-1)^{m/2} is ambiguous; we scan candidates and keep the one(s)
for which the identity holds to machine precision.
"""
import numpy as np
from index_lib import chir_ratio, free_chiral

R = 0.2
x = 0.3
N = 128                      # trapezoid points on |z| = 1
M = 80                       # flux cutoff
z = np.exp(2j * np.pi * (np.arange(N) + 0.5) / N)

# flux-sector integrals, independent of w and of the sign choice
sector = {}
for m in range(-M, M + 1):
    integrand = z ** (-(m + abs(m)) // 2) * chir_ratio(z, m, R, x)
    sector[m] = x ** ((1 - R) * abs(m) / 2) * np.mean(integrand)

signs = {
    "+1":        lambda m: 1.0,
    "(-1)^m":    lambda m: (-1.0) ** m,
    "e^{+i pi m/2}": lambda m: np.exp(1j * np.pi * m / 2),
    "e^{-i pi m/2}": lambda m: np.exp(-1j * np.pi * m / 2),
}

def I_E(w, sgn):
    return sum(sgn(m) * w ** (-m) * sector[m] for m in sector)

ws = np.exp(1j * np.array([0.3, 1.1, 2.5, -0.7]))
RT = (1 - R) / 2
duals = {
    "Chir(w)":   lambda w: free_chiral(w, 0, RT, x),
    "Chir(-w)":  lambda w: free_chiral(-w, 0, RT, x),
    "Chir(1/w)": lambda w: free_chiral(1 / w, 0, RT, x),
    "Chir(-1/w)": lambda w: free_chiral(-1 / w, 0, RT, x),
}

print(f"R = {R}, x = {x}, R_T = (1-R)/2 = {RT}\n")
print(f"{'sign (-1)^(m/2) choice':<18}{'dual':<12} max |I_E/I_dual - 1| over 4 values of w")
best = []
for sname, sgn in signs.items():
    for dname, dual in duals.items():
        err = max(abs(I_E(w, sgn) / dual(w) - 1) for w in ws)
        best.append((err, sname, dname))
        print(f"{sname:<18}{dname:<12} {err:.3e}")

best.sort()
print("\nBest:", best[0])
print("Tail check, |sector(+-M)|:", abs(sector[M]), abs(sector[-M]))
