"""LOG 19 diagnostic: where in (m_L, m_T, m_R) does the weight sit at fixed m_B? Full box, per-sector |contribution|.
Usage: python diag.py n x mB lo hi Ng"""
import sys, numpy as np, optc_defs as O, idx_jax
n, x, mB, lo, hi, Ng = int(sys.argv[1]), float(sys.argv[2]), int(sys.argv[3]), int(sys.argv[4]), int(sys.argv[5]), int(sys.argv[6])
r = O.r; fl = dict(O.fl); fl['eta'] = np.exp(-0.9j); C = (fl['eta']*x**(1 - r))**(n/2)
ch = [(O.one(O.B, 1), {'eta': -1}, r)] + O.rest
rng = range(lo, hi + 1)
rT = range(int(sys.argv[7]), int(sys.argv[8]) + 1) if len(sys.argv) > 8 else rng
secs = [(a, b, mB, d) for a in rng for b in rT for d in rng]
vals = idx_jax.theory_index(4, ch, O.Kmat(1), {(O.B, 'eta'): -.5}, O.FIfull, fl, x, M=0, Ngrid=Ng, bg={'eta': -n}, sectors=secs, per_sector=True, chunk=64)
res = {s: v/C for s, v in zip(secs, vals)}
tot = sum(res.values()); print(f"DIAG n={n} x={x} m_B={mB} box[{lo},{hi}]^3: slice/C={tot:.6e}")
for i, nm in [(0, 'L'), (1, 'T'), (3, 'R')]:
    R_ = rT if nm == 'T' else rng
    marg = {v: sum(abs(z) for s, z in res.items() if s[i] == v) for v in R_}
    print(f"  marginal |.| in m_{nm}: " + " ".join(f"{v}:{marg[v]:.2e}" for v in R_))
for s, z in sorted(res.items(), key=lambda t: -abs(t[1]))[:12]:
    print(f"  top sector {s}: {z:.3e}")
