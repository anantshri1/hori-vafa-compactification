"""LOG 38: JK pole classes of the ladder toric GLSM relative to the wall ray (1,1):
S_int: pairs with charges {(1,0),(0,1)} (eta on the wall is strictly inside); S_A: {(1,0),(1,1)}; S_B: {(0,1),(1,1)};
S_same: identical charges (det 0, excluded). Chambers: A = S_int + S_A, B = S_int + S_B. Test combinations vs Gr(2,4)."""
import numpy as np, itertools, eg as G
def classes(mu, y):
    m = np.linalg.lstsq(G.E.T, mu, rcond=None)[0]; S = {'int': 0, 'A': 0, 'B': 0}; cnt = {'int': 0, 'A': 0, 'B': 0}
    for P in itertools.combinations(range(6), 2):
        QP = G.Qg[:, list(P)]; d = np.linalg.det(QP)
        if abs(d) < 1e-9: continue
        cs = {tuple(G.Qg[:, i].astype(int)) for i in P}
        key = 'int' if cs == {(1, 0), (0, 1)} else ('A' if cs == {(1, 0), (1, 1)} else 'B')
        s = np.linalg.solve(QP.T, -m[list(P)]); t = 1.0
        for i in range(6):
            if i not in P: t *= G.ch(np.exp(2j*np.pi*(G.Qg[:, i] @ s + m[i])), y)
        S[key] += t/abs(d); cnt[key] += 1
    return S, cnt
for tqv in (0.15, 1e-7):
    G.tq = tqv
    for yv in (0.13, 0.31, 0.21+0.04j):
        y = np.exp(2j*np.pi*yv); gr = G.eg_gr(G.M, y); S, cnt = classes(G.mu, y)
        combos = {'int': S['int'], 'int+A': S['int']+S['A'], 'int+(A+B)/2': S['int']+(S['A']+S['B'])/2, 'A': S['A'], 'B': S['B'],
                  'A+B-int': S['A']+S['B']-S['int']}
        print(f"nome={tqv:g} y={yv}: Gr {gr:.6f}; counts {cnt}; " + "; ".join(f"{k}: {v:.6f}" for k, v in combos.items()), flush=True)
