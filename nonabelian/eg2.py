import numpy as np, itertools, eg as G
def terms(mu, y, eta=np.array([1.0, 0.37])):
    m = np.linalg.lstsq(G.E.T, mu, rcond=None)[0]; out = []
    for P in itertools.combinations(range(6), 2):
        QP = G.Qg[:, list(P)]; d = np.linalg.det(QP)
        if abs(d) < 1e-9 or np.any(np.linalg.solve(QP, eta) < 0): continue
        s = np.linalg.solve(QP.T, -m[list(P)]); t = 1.0
        for i in range(6):
            if i not in P: t *= G.ch(np.exp(2j*np.pi*(G.Qg[:, i] @ s + m[i])), y)
        out.append((P, t/abs(d)))
    return out
y = np.exp(2j*np.pi*0.13); a = G.eg_gr(G.M, y); T = terms(G.mu, y)
names = ['X', 'Y/X', 'Z/X', 'W/Y', 'W/Z', 'q/W']
print("EG Gr(2,4) =", np.round(a, 6))
for P, t in T: print(" pole", [names[i] for i in P], np.round(t, 6))
best = min(((abs(sum(T[i][1] for i in S) - a), S) for S in itertools.combinations(range(8), 6)))
print("best 6-subset |sum - EG_Gr| =", best[0], "subset", [[names[i] for i in T[j][0]] for j in best[1]])
# also individual Gr(2,4) vacuum terms
