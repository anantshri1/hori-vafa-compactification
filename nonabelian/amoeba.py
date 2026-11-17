"""LOG 26, step 1-2: leading large-flux asymptotics. At u = eps|nu| fixed each chiral gives, on the real torus,
the pure phase exp(2i Im Li2(s a zeta)/eps), a = x^{|nu|+1}; CS/FI/zero-point monomials z^{ez} are phases e^{i ez.phi}
with ez ~ 1/eps. A sector is O(1) iff the phase is stationary somewhere on T^4:
   sum_i Q_ia [-2 log|1 - s_i a_i e^{i theta_i}|] + eps ez_a = 0,   theta_i = Q_i.phi + arg(flavour_i),
otherwise exponentially small. => support in flux space = solvability region (amoeba-type). Test vs hfjob_sharp.log."""
import sys, numpy as np
from fractions import Fraction as Fr
from scipy.optimize import least_squares
import optc_defs as O
def sector_data(m, x, n, fl):
    chirals = [(O.one(O.B, 1), {'eta': -1}, O.r)] + O.rest
    K = O.Kmat(1); Kf = {(O.B, 'eta'): -.5}; FI = O.FIfull; bg = {'eta': -n}; nn = 4
    ez = [Fr(0)]*nn
    for a in range(nn): ez[a] += sum(Fr(K[a][b])*int(m[b]) for b in range(nn))
    for (a, f), lev in Kf.items(): ez[a] += Fr(lev)*bg.get(f, 0)
    for a, fls in FI.items():
        for f in fls: ez[a] -= bg.get(f, 0)
    Q, nus, ph = [], [], []
    for Qc, flc, R in chirals:
        nu = int(sum(Qc[a]*m[a] for a in range(nn)) + sum(c*bg.get(f, 0) for f, c in flc.items()))
        for a in range(nn): ez[a] += Fr(-Qc[a]*abs(nu), 2)
        Q.append(Qc); nus.append(nu); ph.append(sum(c*np.angle(fl[f]) for f, c in flc.items()))
    # FI fugacities are phases^{-m}: no phi dependence. Flavour moduli -> 1 as x -> 1 (dropped at leading order).
    return np.array(Q, float), np.array(nus), np.array(ph), np.array([float(e) for e in ez])
def residual(m, x, n, fl, starts=12, rng=np.random.default_rng(0)):
    eps = -2*np.log(x)
    Q, nus, ph, ez = sector_data(m, x, n, fl)
    a = x**(np.abs(nus) + 1.0); s = np.where(nus % 2 == 0, 1.0, -1.0)
    def f(phi):
        th = Q @ phi + ph
        g = -2*np.log(np.abs(1 - s*a*np.exp(1j*th)))
        return Q.T @ g + eps*ez
    best = np.inf
    for _ in range(starts):
        sol = least_squares(f, rng.uniform(-np.pi, np.pi, 4))
        best = min(best, np.linalg.norm(sol.fun))
        if best < 1e-6: break
    return best
if __name__ == "__main__":
    n, x, mB, lo, hi = int(sys.argv[1]), float(sys.argv[2]), int(sys.argv[3]), int(sys.argv[4]), int(sys.argv[5])
    fl = dict(O.fl); fl['eta'] = np.exp(-0.9j)
    out = []
    for mL in range(lo, hi + 1):
        ok = 0; tot = 0
        for mT in (0, 2, 4):
            for mR in range(lo, hi + 1):
                tot += 1; ok += residual((mL, mT, mB, mR), x, n, fl) < 1e-3
        out.append((mL, ok/tot))
    print(f"n={n} x={x} m_B={mB}: solvable fraction vs m_L: " + " ".join(f"{l}:{p:.2f}" for l, p in out), flush=True)
