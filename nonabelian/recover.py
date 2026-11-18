import numpy as np, itertools, pairsum as PS
rng = np.random.default_rng(9); q = PS.q; c = PS.c; lad = PS.lad
for trial in range(4):
    a, bL, bR = rng.normal(size=3)*0.6 + 1j*rng.normal(size=3)*0.6
    mu = np.array([a, bL, bR, a]); print("mu_B=mu_T", np.round(a,3), "mu_L", np.round(bL,3), "mu_R", np.round(bR,3))
    P = PS.crit(c, mu, tries=8000)
    if len(P) != 6: print("found", len(P)); continue
    t = np.array([q*np.exp(lad[5] @ np.log(p)) for p in P])
    err, pairs, S = PS.best_pairing(t)
    # recover 4 numbers: pair sums; pick pairs (0,x): a0 from three sums
    T = {}
    for i, j in itertools.combinations(range(6), 2): pass
    # brute force: choose 4 numbers x with all 6 pair sums = t (as multiset)
    best = (np.inf, None)
    for perm in itertools.permutations(range(6)):
        s12, s13, s14, s23 = t[perm[0]], t[perm[1]], t[perm[2]], t[perm[3]]
        x1 = (s12 + s13 - s23)/2; x = np.array([x1, s12 - x1, s13 - x1, s14 - x1])
        pred = np.sort_complex(np.array([x[i]+x[j] for i, j in itertools.combinations(range(4), 2)]))
        e = np.linalg.norm(pred - np.sort_complex(t))
        if e < best[0]: best = (e, x)
    x = best[1]; ex = np.poly(x)
    # candidate masses from column FIs: differences (mu_L, mu_B, mu_R) in some order, sign s
    out = []
    for order in itertools.permutations([bL, a, bR]):
        for s in (1, -1):
            m = s*np.concatenate([[0], np.cumsum(order)]); m = m - m.mean() + x.mean()
            em = np.poly(m)
            out.append((np.linalg.norm(em[1:4] - ex[1:4]), [np.round(z, 3) for z in order], s, ex[4] - em[4]))
    o = min(out, key=lambda z: z[0])
    print(f"pair-sum err {err:.1e}; recover err {best[0]:.1e}; best mass map err e1..e3 {o[0]:.2e} order {o[1]} sign {o[2]}; e4 diff {o[3]:.4f} vs q {q:.4f}", flush=True)
