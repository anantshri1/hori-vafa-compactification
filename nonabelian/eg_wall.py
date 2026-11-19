"""LOG 37 (option 3). Chambers of the ladder's toric GLSM: rays (1,0),(0,1),(1,1); the two geometric phases are
cone((1,0),(1,1)) and cone((1,1),(0,1)), separated by the ray (1,1) <-> t1 = t2 <-> rho = 1 (the reduction's value).
(i) EG on both sides of the wall; (ii) chi_y (nome -> 0) of the toric phase vs Gr(2,4); compare the difference with
chi_y(P^1) (chi_y(P^1) - 1), the expected contribution of a small resolution of a curve C = P^1 of conifold points."""
import numpy as np, eg as G
rng = np.random.default_rng(5)
for tqv in (0.15, 1e-7):
    G.tq = tqv
    for yv in (0.13, 0.31):
        y = np.exp(2j*np.pi*yv)
        a = G.eg_gr(G.M, y)
        b1 = G.eg_tor(G.mu, y, eta=np.array([1.0, 0.37])); b2 = G.eg_tor(G.mu, y, eta=np.array([0.37, 1.0]))
        Y = -y**-1  # chi_y variable in this normalisation (fit below)
        print(f"nome={tqv:g} y=e^(2pi i {yv}): Gr {a:.6f} | toric chamber A {b1:.6f} | chamber B {b2:.6f} | toric-Gr {b1-a:.6f}", flush=True)
# chi_y polynomials: fit EG(nome->0) as Laurent polynomial in y^{1/2}
G.tq = 1e-9
ys = np.exp(2j*np.pi*np.linspace(0.05, 0.95, 9))
def fit(f):
    A = np.array([[y**(k/2) for k in range(-4, 5)] for y in ys]); return np.linalg.lstsq(A, np.array([f(y) for y in ys]), rcond=None)[0]
cg = fit(lambda y: G.eg_gr(G.M, y)); ct = fit(lambda y: G.eg_tor(G.mu, y))
print("chi_y Gr(2,4)   coeffs y^{k/2}, k=-4..4:", np.round(cg.real, 6))
print("chi_y toric     coeffs y^{k/2}, k=-4..4:", np.round(ct.real, 6))
print("difference                              :", np.round((ct - cg).real, 6))
