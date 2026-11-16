"""LOG item 12: 3d flux sectors of U(2)_{(0,-2)}[4,0] by residues (K-theoretic vortex sums), vs Gr(2,4) sectors by residues.
Poles: sigma_a = -mu_{j_a} - i(R/2 + |m_a|/2 + k_a) (closing below, as in gr24_res). Residues by 4-direction averaging.
3d sector = 0.5 (eps/2pi)^2 (-2 pi i)^2 sum Res_sigma; 2d sector = 0.5 (-2 pi i)^2 sum Res_sigma. ALL (j1,j2) pairs kept."""
import numpy as np
from scipy.special import loggamma
from idx import qpoch
import higgs_flux as H
R, mu, theta = H.R, H.mu, H.theta
def f3(S1, S2, eps, m, n):
    x = np.exp(-eps/2); q = x*x
    val = np.ones_like(S1, dtype=complex); logmon = np.zeros_like(S1, dtype=complex)
    for a, S in enumerate((S1, S2)):
        nu = m[a]
        for mj in mu:
            ze = q**(-1j*(S + mj))
            val *= qpoch(x**(2 - R + abs(nu))/ze, q)/qpoch(ze*x**(R + abs(nu)), q)
            logmon += -abs(nu)/2*(np.log(ze) + (R - 1)*np.log(x)) - abs(nu)/2*1j*np.pi*nu
    d = abs(m[0]-m[1]); s = (-1)**(m[0]-m[1])
    z12 = s*q**(-1j*(S1-S2))
    val *= (1 - s*z12*x**d)*(1 - s/z12*x**d)*x**(-d)
    Km = -(m[0] + m[1])
    for a, S in enumerate((S1, S2)):
        val *= ((-1)**m[a]*q**(-1j*S))**(Km - n)
    return val*np.exp(logmon)*np.exp(-1j*theta*(m[0] + m[1]))
def f2(S1, S2, m, kappa):
    lv = 0
    for a, S in enumerate((S1, S2)):
        for mj in mu:
            a_ = S + mj
            lv = lv + loggamma(R/2 - 1j*a_ + abs(m[a])/2) - loggamma(1 - R/2 + 1j*a_ + abs(m[a])/2)
    return np.exp(lv)*((S1-S2)**2 + (m[0]-m[1])**2/4)*np.exp(1j*kappa*(S1 + S2) - 1j*theta*(m[0] + m[1]))
D = np.array([1, -1, 1j, -1j])
def res2(f, s1, s2, h):
    d1, d2 = np.meshgrid(D*h, D*h, indexing='ij')
    return np.mean(f(s1 + d1, s2 + d2)*d1*d2)
def sectors(eps, m, n, kap, K=15, h=1e-4, samej=True):
    J1, J2, K1, K2 = [g.ravel() for g in np.meshgrid(range(4), range(4), range(K), range(K), indexing='ij')]
    s1 = -mu[J1] - 1j*(R/2 + abs(m[0])/2 + K1); s2 = -mu[J2] - 1j*(R/2 + abs(m[1])/2 + K2)
    d1, d2 = [g.ravel() for g in np.meshgrid(D*h, D*h, indexing='ij')]
    A = s1[:, None] + d1[None]; B = s2[:, None] + d2[None]; w = (d1*d2)[None]
    r3 = np.mean(f3(A, B, eps, m, n)*w, axis=1)*eps**2
    r2 = np.mean(f2(A, B, m, kap)*w, axis=1)*(2*np.pi)**2
    if not samej: keep = J1 != J2; r3, r2 = r3*keep, r2*keep
    per = {}
    for a in range(4):
        for b in range(4):
            sel = (J1 == a) & (J2 == b); per[(a, b)] = np.array([r3[sel].sum(), r2[sel].sum()])
    return 0.5*r3.sum(), 0.5*r2.sum(), per
if __name__ == "__main__":
    import sys
    for eps in [float(e) for e in sys.argv[1:]] or (0.1, 0.05):
        n = int(round((1.0 - 8*np.log(1 - np.exp(-eps)))/eps/2))*2; kap = -eps*n - 8*np.log(1 - np.exp(-eps))
        rows = []
        for m in [(0, 0), (1, 0), (0, 1), (1, 1), (2, 0), (-1, 0), (1, -1)]:
            t3, t2, per = sectors(eps, m, n, kap)
            rows.append((m, t3/t2, per))
        base = rows[0][1]
        print(f"eps={eps} n={n} kappa={kap:.3f}  base ratio {base:.4e}")
        for m, rt, per in rows:
            pv = per[(0, 1)]; pd = per[(0, 0)]
            print(f"  m={m}: ratio/base {rt/base:.4f}   vac(0,1) {pv[0]/pv[1]/base:.4f}   same-j(0,0) 3d {abs(pd[0]):.2e} 2d {abs(pd[1]):.2e}", flush=True)
