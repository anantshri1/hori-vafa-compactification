"""GPU scan (LOG 19): slice profile of the mirror index (eta-flux -n) along m_B, with L,T,R windows co-moving:
centre (round(cL*mB), round(cT*mB), mB, round(cR*mB)), half-width M. Usage: python scan.py n x M Ng cL cT cR"""
import sys, time, numpy as np, optc_defs as O, idx_jax
n, x, M, Ng = int(sys.argv[1]), float(sys.argv[2]), int(sys.argv[3]), int(sys.argv[4])
cL, cT, cR = map(float, sys.argv[5:8])
r = O.r; fl = dict(O.fl); fl['eta'] = np.exp(-0.9j); C = (fl['eta']*x**(1 - r))**(n/2)
ch = [(O.one(O.B, 1), {'eta': -1}, r)] + O.rest
t0 = time.time(); tot = 0j
for mB in range(2, -n - 3, -1):
    c = (int(round(cL*mB)), int(round(cT*mB)), int(round(cR*mB)))
    secs = [(c[0]+a, c[1]+b, mB, c[2]+d) for a in range(-M, M+1) for b in range(-M, M+1) for d in range(-M, M+1)]
    I = idx_jax.theory_index(4, ch, O.Kmat(1), {(O.B, 'eta'): -.5}, O.FIfull, fl, x, M=0, Ngrid=Ng, bg={'eta': -n}, sectors=secs)/C
    tot += I
    print(f"n={n} x={x} M={M} Ng={Ng} c=({cL},{cT},{cR}) m_B={mB:+d} m_B/n={mB/n:+.3f} slice/C={I:.6e} |.|={abs(I):.4e} t={time.time()-t0:.0f}s", flush=True)
print(f"TOTAL n={n} x={x} M={M} c=({cL},{cT},{cR}): sum/C={tot:.8e}", flush=True)
