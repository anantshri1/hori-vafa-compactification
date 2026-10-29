"""
Step 1f: CFS limit 2 analogue for U(1)_{1/2} + Phi <-> free chiral:
topological fugacity w = e^{i theta} FIXED (not -> 1), background flux n = 0.

Exact sector formulas (step 1c):
  neutral m = -n:  S = w^n x^{(1-R)n/2} F_n,
  charged m = +n:  S = w^{-n} x^{(1-R)n/2} (-1)^n x^{(R+n)n} G_n.
Dual: Chir(w, 0, R_T). Leading asymptotics (step 0 formula A with |a| = 1):
  Chir ~ exp(2 i Cl2(theta)/eps) * (1 - w)^{-(1-R_T)/2} (1 - 1/w)^{-(1-R_T)/2},  Cl2(theta) = Im Li2(e^{i theta}).
Question: where does the electric weight sit, and does the dual asymptotic formula hold?
"""
import mpmath as mp
mp.mp.dps = 50

R, theta = mp.mpf('0.2'), mp.mpf('0.7')
RT = (1 - R) / 2
w = mp.e**(1j * theta)

def P_list(q, K):
    out = [mp.mpf(1)]
    for k in range(1, K + 1):
        out.append(out[-1] * (1 - q**k))
    return out

for eps in [mp.mpf(e) for e in ['0.4', '0.2', '0.1', '0.05']]:
    q = mp.e**(-eps); x = mp.e**(-eps/2)
    K = int(12/float(eps)) + 60
    P = P_list(q, K + 400)
    Nmax = int(60/float(eps*(1-R)/4))
    S = {}
    for n in range(0, Nmax):
        u2 = x**(2*n)
        S[-n] = w**n * x**((1-R)*n/2) * mp.fsum((-1)**k * q**(k*(k+1)/2) * u2**k / P[k]**2 for k in range(K))
    for n in range(1, int(mp.sqrt(80/eps)) + 5):
        u2 = x**(2*n)
        S[n] = w**(-n) * x**((1-R)*n/2) * (-1)**n * x**((R+n)*n) * \
               mp.fsum((-1)**j * q**(j*(j+1)/2) * u2**j / (P[j]*P[j+n]) for j in range(K))
    tot = mp.fsum(S.values())
    dual = mp.qp(x**(2-RT)/w, q) / mp.qp(w*x**RT, q)
    lead = mp.e**(2j*mp.im(mp.polylog(2, w))/eps) * ((1-w)*(1-1/w))**(-(1-RT)/2)
    fixed = abs(mp.fsum(S[m] for m in S if abs(m) <= 3)) / abs(tot)
    absS = {m: abs(S[m]) for m in S}
    meanabs = mp.fsum(abs(m)*a for m, a in absS.items()) / mp.fsum(absS.values())
    mpeak = max(absS, key=absS.get)
    print(f"eps={float(eps):<5} |I_E/dual-1|={float(abs(tot/dual-1)):.1e}  "
          f"|dual/leading-1|={float(abs(dual/lead-1)):.2e}  share |m|<=3: {float(fixed):.3f}  "
          f"eps<|m|>={float(eps*meanabs):.2f}  peak |S_m| at m={mpeak} (eps m={float(eps*mpeak):.2f})")
