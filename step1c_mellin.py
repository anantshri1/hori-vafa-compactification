"""
Step 1c: analytic small-circle limit of U(1)_{1/2} + chiral, and numerical confirmation.

Exact (verified): neutral sector m = -n contributes
    S_{-n} = w^n x^{(1-R)n/2} F_n,   F_n = sum_k (-1)^k q^{k(k+1)/2} u_n^{2k} / (q;q)_k^2,  u_n = x^n
and charged sectors m = +n (n >= 1), by the same q-binomial manipulation (coefficient of z^n),
    S_{+n} = w^{-n} x^{(1-R)n/2} (-1)^n x^{(R+n)n} G_n,
    G_n = sum_j (-1)^j q^{j(j+1)/2} u_n^{2j} / ((q;q)_j (q;q)_{j+n}).

Small eps (q = e^{-eps}, x = e^{-eps/2}, w = q^{-i mu}):
  (q;q)_k ~ k! eps^k  =>  F_n ~ sum_k (-1)^k (u/eps)^{2k}/k!^2 = J_0(2 u_n / eps).
  The weight sits at u_n ~ eps, i.e. |m| ~ (2/eps) log(1/eps): NOT at fixed flux.
  Set t = u_n/eps, dn = (2/eps) dt/t, x^{(1-R)n/2} = (eps t)^{R_T}, w^n = (eps t)^{-2i mu}, R_T = (1-R)/2:
     I_E ~ (2/eps) eps^{s} int_0^inf dt t^{s-1} J_0(2t),   s = R_T - 2 i mu
         = eps^{s-1} Gamma(s/2)/Gamma(1-s/2)
         = eps^{R_T-1-2i mu} Gamma(R_T/2 - i mu)/Gamma(1 - R_T/2 + i mu)
  which is exactly the 2d free chiral of the dual (step 0, test C).
"""
import mpmath as mp
mp.mp.dps = 60

R, mu = mp.mpf('0.2'), mp.mpf('0.35')
RT = (1 - R) / 2
s = RT - 2j * mu

def qpoch_list(q, K):
    out = [mp.mpf(1)]
    for k in range(1, K + 1):
        out.append(out[-1] * (1 - q**k))
    return out

def F(u2, q, P, K):
    return mp.fsum((-1)**k * q**(k*(k+1)/2) * u2**k / P[k]**2 for k in range(K))

def G(u2, q, n, P, K):
    return mp.fsum((-1)**j * q**(j*(j+1)/2) * u2**j / (P[j] * P[j+n]) for j in range(K))

print("Check 1: F_n vs J_0(2 x^{n+1}/eps) at eps = 0.05 (Bessel regime of each neutral sector)")
eps = mp.mpf('0.05'); q = mp.e**(-eps); x = mp.e**(-eps/2); P = qpoch_list(q, 400)
for n in [40, 80, 120, 160, 240]:
    u = x**n
    Fn = F(u*u, q, P, 300)
    print(f"  n={n:4d}  t=x^n/eps={float(u/eps):6.2f}   F_n={mp.nstr(Fn,8):>14}   "
          f"J0(2x^(n+1)/eps)={mp.nstr(mp.besselj(0, 2*x**(n+1)/eps),8):>14}")

print("\nCheck 2: full electric index vs the 2d prediction  eps^{s-1} G(s/2)/G(1-s/2)")
print("  (and the size of the charged m>0 sectors relative to the total)")
for eps in [mp.mpf(e) for e in ['0.4', '0.2', '0.1', '0.05', '0.025']]:
    q = mp.e**(-eps); x = mp.e**(-eps/2); w = mp.e**(1j*eps*mu)
    K = int(12/float(eps)) + 60
    P = qpoch_list(q, K + 400)
    Nmax = int(60/float(eps*(1-R)/4))
    neutral = mp.fsum(w**n * x**((1-R)*n/2) * F(x**(2*n), q, P, K) for n in range(0, Nmax))
    charged = mp.fsum(w**(-n) * x**((1-R)*n/2) * (-1)**n * x**((R+n)*n) * G(x**(2*n), q, n, P, K)
                      for n in range(1, int(mp.sqrt(80/eps)) + 5))
    tot = neutral + charged
    pred = eps**(s - 1) * mp.gamma(s/2) / mp.gamma(1 - s/2)
    # exact dual (3d) for reference: Chir(w, 0, R_T)
    dual = mp.qp(x**(2-RT)/w, q) / mp.qp(w*x**RT, q)
    print(f"  eps={float(eps):<6} |I_E/dual-1|={float(abs(tot/dual-1)):.1e}   "
          f"|I_E/pred-1|={float(abs(tot/pred-1)):.3e}   |charged|/|total|={float(abs(charged/tot)):.3e}")
