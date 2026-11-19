"""LOG 36: elliptic genera. (A) Gr(2,4) GLSM, U(2) + 4 fundamentals (R=0), equivariant masses M_j, JK residues.
(B) the ladder's abelian GLSM ('inverse Hori-Vafa'): one chiral per ladder monomial, gauge U(1)^2 = relations of the
exponents, Q1 = (1,1,0,1,0,1), Q2 = (1,0,1,0,1,1) on (X, Y/X, Z/X, W/Y, W/Z, q/W); chiral masses m_v with E^T m = mu,
mu = (mu_B, mu_L, mu_R, mu_T) = (M3-M2, M2-M1, M4-M3, M3-M2). (2,2) chiral (R=0): theta1(y^{-1}X)/theta1(X).
JK: each pole chiral contributes 1 (cancels the vector zero modes); |det| for abelian."""
import numpy as np, itertools
tq = 0.15  # nome
def th1(z, N=40):  # theta1 with multiplicative argument z = e^{2 pi i u}
    n = np.arange(N); u = np.log(z)/(2j*np.pi)
    return 2*np.sum((-1)**n*tq**((n+.5)**2/2)*np.sin((2*n+1)*np.pi*u))
def ch(X, y): return th1(X/y)/th1(X)
def eg_gr(M, y):
    Mu = np.exp(2j*np.pi*M); tot = 0
    for j1, j2 in itertools.permutations(range(4), 2):
        z = [Mu[j1], Mu[j2]]; t = 1.0
        for a in range(2):
            for k in range(4):
                if k == [j1, j2][a]: continue
                t *= ch(Mu[k]/z[a], y)        # chiral X = z_a / M_k ?  use M_k/z_a convention for both sides
        for a, b in [(0, 1), (1, 0)]:
            t *= th1(z[a]/z[b])/th1(z[a]/(z[b]*y))   # W-bosons
        tot += t/2
    return tot
E = np.array([(1,0,0,0),(-1,1,0,0),(-1,0,1,0),(0,-1,0,1),(0,0,-1,1),(0,0,0,-1)], float)  # (X,Y,Z,W) = (B,L,R,T)
Qg = np.array([(1,1,0,1,0,1), (1,0,1,0,1,1)], float)
def eg_tor(mu, y, eta=np.array([1.0, 0.37])):
    m = np.linalg.lstsq(E.T, mu, rcond=None)[0]           # particular solution; gauge shifts irrelevant
    tot = 0
    for P in itertools.combinations(range(6), 2):
        QP = Qg[:, list(P)]; d = np.linalg.det(QP)
        if abs(d) < 1e-9: continue
        c = np.linalg.solve(QP, eta)
        if np.any(c < 0): continue                          # JK: eta in the cone of the pole charges
        # pole: Q_i . s + m_i = integer for i in P  -> s = -QP^{-T} m_P (+ lattice shifts: |det| solutions)
        sols = []
        for k in itertools.product(range(int(round(abs(d)))), repeat=2):
            s = np.linalg.solve(QP.T, -m[list(P)] + np.array(k))
            sols.append(s)
        for s in sols:
            t = 1.0
            for i in range(6):
                if i in P: continue
                X = np.exp(2j*np.pi*(Qg[:, i] @ s + m[i]))
                t *= ch(X, y)
            tot += t/abs(d)/len(sols)*len(sols)/abs(d) if False else t/abs(d)
    return tot
rng = np.random.default_rng(0)
M = rng.normal(0, .15, 4) + 1j*rng.normal(0, .05, 4); M -= M.mean()
mu = np.array([M[2]-M[1], M[1]-M[0], M[3]-M[2], M[2]-M[1]])
for y in (np.exp(2j*np.pi*1e-7), np.exp(2j*np.pi*0.13), np.exp(2j*np.pi*(0.21+0.04j))):
    a, b = eg_gr(M, y), eg_tor(mu, y)
    print(f"y={y:.4f}: EG Gr(2,4) = {a:.8f}   EG ladder-GLSM = {b:.8f}   ratio {b/a if abs(a)>1e-12 else np.nan:.6f}", flush=True)
