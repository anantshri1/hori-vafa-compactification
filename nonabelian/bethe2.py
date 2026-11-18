import numpy as np, bethe as Bm
from bethe import solve, dedupe, rng, y, w, Qs, K, f, t
for sig in (+1, -1):
    def F_el(l):
        z = np.exp(l); s = z[0]*z[1]
        return np.array([-(z[a]**(2*sig)) - w*s*np.prod(1 - y*z[a]) for a in range(2)])/(1 + np.abs(z)**4 + np.abs(z)**-4)
    el = dedupe([s for s in solve(F_el, 2, 4000) if abs(s[0]-s[1]) > 1e-5*(1+abs(s[0]))], lambda s: np.sort_complex(s))
    Kp = K + sig*0.5*Qs.T @ Qs
    def F_mir(l):
        X = f*np.exp(Qs @ l)
        lhs = np.exp(Kp @ l)*np.prod((1 - X)[:, None]**(-Qs), axis=0)
        return (lhs - t)/(1 + np.abs(lhs))
    ms = dedupe(solve(F_mir, 4, 8000), lambda s: s)
    print(f"sigma={sig:+d}: electric {len(el)}, mirror {len(ms)}", flush=True)
