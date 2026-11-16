"""Gr(2,4) GLSM S^2 flux sectors by residues (kappa < 0, contour closed below): poles of
G(R/2 - i(sigma_a+mu_j) + |m_a|/2) at sigma_a = -mu_{j_a} - i(R/2 + |m_a|/2 + k_a), j_1 != j_2 (Vandermonde kills j_1=j_2).
Each pole: 2 pi (-1)^k / k!  times the rest; Weyl 1/2 cancels the two orderings of (j_1, j_2) so we sum ordered pairs /2."""
import numpy as np
from scipy.special import loggamma, gammaln
import higgs_flux as H
R, mu, theta = H.R, H.mu, H.theta
def sector2d_res(m, kappa, K=40):
    tot = 0j
    for j1 in range(4):
        for j2 in range(4):
            if j1 == j2: continue
            for k1 in range(K):
                for k2 in range(K):
                    js, ks = (j1, j2), (k1, k2)
                    sig = [-mu[js[a]] - 1j*(R/2 + abs(m[a])/2 + ks[a]) for a in range(2)]
                    lv = 0j
                    for a in range(2):
                        for j in range(4):
                            arg = sig[a] + mu[j]
                            if j == js[a]:
                                lv += np.log(2*np.pi) + ks[a]*np.log(-1+0j) - gammaln(ks[a] + 1)
                            else:
                                lv += loggamma(R/2 - 1j*arg + abs(m[a])/2)
                            lv -= loggamma(1 - R/2 + 1j*arg + abs(m[a])/2)
                    vd = (sig[0]-sig[1])**2 + (m[0]-m[1])**2/4
                    tot += np.exp(lv)*vd*np.exp(1j*kappa*(sig[0] + sig[1]) - 1j*theta*(m[0] + m[1]))
    return 0.5*tot/(2*np.pi)**0   # measure d^2 sigma (no 1/2pi) as in sector2d
if __name__ == "__main__":
    import sys
    for eps in (0.1, 0.05):
        n = int(round((1.0 - 8*np.log(1 - np.exp(-eps)))/eps/2))*2; kap = -eps*n - 8*np.log(1 - np.exp(-eps))
        ms = [(0, 0), (1, 0), (1, 1), (2, 0), (2, 1), (-1, 0)]
        rr = []
        for m in ms:
            rr.append(H.sector3d(eps, m, n)/sector2d_res(m, kap))
        print(f"eps={eps} kappa={kap:.3f}:", ' '.join(f"{z/rr[0]:.4f}" for z in rr), f" base {rr[0]:.3e}", flush=True)
