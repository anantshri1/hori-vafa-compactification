"""LOG 39 (option 1b). Fixed-point (JK residue) terms: 8 for the toric small resolution, 6 for Gr(2,4).
Step 1: which terms coincide exactly (same tangent weights)? Step 2: is D = EG_res - EG_Gr carried entirely by the
unmatched ones, and do those form the local geometry of the exceptional locus?"""
import numpy as np, itertools, eg as G
def gr_terms(M, y):
    Mu = np.exp(2j*np.pi*M); out = {}
    for j1, j2 in itertools.combinations(range(4), 2):
        tot = 0
        for (a1, a2) in [(j1, j2), (j2, j1)]:
            z = [Mu[a1], Mu[a2]]; t = 1.0
            for a in range(2):
                for k in range(4):
                    if k != [a1, a2][a]: t *= G.ch(Mu[k]/z[a], y)
            for a, b in [(0, 1), (1, 0)]: t *= G.th1(z[a]/z[b])/G.th1(z[a]/(z[b]*y))
            tot += t/2
        out[(j1, j2)] = tot
    return out
def tor_terms(mu, y, eta=np.array([1.0, 0.37])):
    m = np.linalg.lstsq(G.E.T, mu, rcond=None)[0]; out = {}
    nm = ['X', 'Y/X', 'Z/X', 'W/Y', 'W/Z', 'q/W']
    for P in itertools.combinations(range(6), 2):
        QP = G.Qg[:, list(P)]; d = np.linalg.det(QP)
        if abs(d) < 1e-9 or np.any(np.linalg.solve(QP, eta) < 0): continue
        s = np.linalg.solve(QP.T, -m[list(P)]); t = 1.0
        for i in range(6):
            if i not in P: t *= G.ch(np.exp(2j*np.pi*(G.Qg[:, i] @ s + m[i])), y)
        out[(nm[P[0]], nm[P[1]])] = t/abs(d)
    return out
G.tq = 0.15
for yv in (0.13, 0.31):
    y = np.exp(2j*np.pi*yv); g = gr_terms(G.M, y); t = tor_terms(G.mu, y)
    print(f"y=e^(2pi i {yv})")
    used = set()
    for kt, vt in t.items():
        match = [kg for kg, vg in g.items() if abs(vt - vg) < 1e-9*(1 + abs(vg))]
        print(f"  toric {kt}: {vt:.6f}   equals Gr vacuum {match}")
        used |= set(match)
    un_g = {k: v for k, v in g.items() if k not in used}
    un_t = sum(v for k, v in t.items() if not any(abs(v - gv) < 1e-9*(1+abs(gv)) for gv in g.values()))
    print(f"  unmatched Gr vacua {list(un_g)}: sum {sum(un_g.values()):.6f};  unmatched toric sum {un_t:.6f};  D = {sum(t.values())-sum(g.values()):.6f}")
