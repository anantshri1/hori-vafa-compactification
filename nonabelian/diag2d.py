"""LOG 25: sharpening test. At fixed m_B and m_T in a list, full 2d (m_L, m_R) grid; prints the m_L marginal (summed over m_R
and the m_T list) and the diagonal profile, plus an aliasing check (3 extreme sectors at Ng and Ng+8).
Usage: python diag2d.py n x mB Llo Lhi Ng mT1,mT2,..."""
import sys, time, numpy as np, optc_defs as O, idx_jax
n, x, mB, lo, hi, Ng = int(sys.argv[1]), float(sys.argv[2]), int(sys.argv[3]), int(sys.argv[4]), int(sys.argv[5]), int(sys.argv[6])
mTs = [int(t) for t in sys.argv[7].split(',')]
r = O.r; fl = dict(O.fl); fl['eta'] = np.exp(-0.9j); C = (fl['eta']*x**(1 - r))**(n/2)
ch = [(O.one(O.B, 1), {'eta': -1}, r)] + O.rest
def run(secs, ng): return idx_jax.theory_index(4, ch, O.Kmat(1), {(O.B, 'eta'): -.5}, O.FIfull, fl, x, M=0, Ngrid=ng, bg={'eta': -n}, sectors=secs, per_sector=True, chunk=32)
corners = [(lo, mTs[0], mB, lo), (hi, mTs[-1], mB, hi), (lo, mTs[0], mB, hi)]
a, b = run(corners, Ng), run(corners, Ng + 8)
print(f"[{time.strftime('%H:%M:%S')}] ALIAS n={n} x={x} m_B={mB}: " + " ".join(f"{s}:{abs(u/C):.2e}/{abs(v/C):.2e}" for s, u, v in zip(corners, a, b)), flush=True)
rng = range(lo, hi + 1)
secs = [(l, t, mB, rr) for t in mTs for l in rng for rr in rng]
res = {s: v/C for s, v in zip(secs, run(secs, Ng))}
print(f"[{time.strftime('%H:%M:%S')}] D2 n={n} x={x} m_B={mB} m_T={mTs} L,R in [{lo},{hi}] Ng={Ng}: partial/C={sum(res.values()):.6e}")
marg = {l: sum(abs(z) for s, z in res.items() if s[0] == l) for l in rng}
print("  m_L marginal: " + " ".join(f"{l}:{marg[l]:.2e}" for l in rng))
print("  diagonal |.| at m_T=" + str(mTs[len(mTs)//2]) + ": " + " ".join(f"{l}:{abs(res[(l, mTs[len(mTs)//2], mB, l)]):.2e}" for l in rng), flush=True)
