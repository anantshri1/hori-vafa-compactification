"""LOG 13: do the j1=j2 residues cancel in the sum over flux sectors? Sum over m in [-M,M]^2 of same-j and distinct-j parts."""
import numpy as np, res3d as T
for eps in (0.1, 0.05):
    n = int(round((1.0 - 8*np.log(1 - np.exp(-eps)))/eps/2))*2; kap = -eps*n - 8*np.log(1 - np.exp(-eps))
    for M in (2, 4, 6):
        S = np.zeros((2, 2), complex)  # [3d,2d] x [samej, distinct]
        for m1 in range(-M, M+1):
            for m2 in range(-M, M+1):
                t3, t2, per = T.sectors(eps, (m1, m2), n, kap, K=12)
                for i in range(2):
                    sj = 0.5*sum(per[(j, j)][i] for j in range(4)); tot = (t3, t2)[i]
                    S[i] += [sj, tot - sj]
        print(f"eps={eps} M={M}: 3d samej/distinct {abs(S[0,0]/S[0,1]):.2e}  2d {abs(S[1,0]/S[1,1]):.2e}  full 3d/2d {S[0].sum()/S[1].sum():.5e}", flush=True)
