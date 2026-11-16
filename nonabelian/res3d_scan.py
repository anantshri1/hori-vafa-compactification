import numpy as np, res3d as T
for eps in (0.2, 0.1, 0.05, 0.025):
    n = int(round((1.0 - 8*np.log(1 - np.exp(-eps)))/eps/2))*2; kap = -eps*n - 8*np.log(1 - np.exp(-eps))
    K = 15 if eps > 0.06 else 12
    out = {}
    for m in [(0, 0), (1, 0), (-1, 0), (1, 1), (1, -1), (2, 0)]:
        t3, t2, per = T.sectors(eps, m, n, kap, K=K)
        out[m] = (t3, t2, per)
    b = out[(0, 0)]
    print(f"eps={eps} n={n} kappa={kap:.3f}", flush=True)
    for m, (t3, t2, per) in out.items():
        tot = (t3/t2)/(b[0]/b[1])
        vac = [(per[v][0]/per[v][1])/(b[2][v][0]/b[2][v][1]) for v in [(0, 1), (1, 2), (2, 3), (0, 3)]]
        sj3 = sum(per[(j, j)][0] for j in range(4))/(2*t3); sj2 = sum(per[(j, j)][1] for j in range(4))/(2*t2)
        print(f"  m={m}: tot {tot:.3f}  per-vac " + " ".join(f"{v:.3f}" for v in vac) + f"  samej frac 3d {abs(sj3):.2e} 2d {abs(sj2):.2e}", flush=True)
