"""B in A's good Higgs limit: y = q^{-i mu}, wA = e^{i theta}, nA = (kap + 4 log(1-q))/eps (even).
B sees: wB = y^2 = q^{-2i mu} (FI scaled), a = e^{i theta/2}, axial flux nax = nA/2 (large).
Diagnostics: exact identity with contact term; distribution of B's gauge flux m; u-values of Q, Qt."""
import numpy as np
from pair_limit import sectors
from tsu2_lib import qinf

rA, mu, theta, kap = 0.3, 0.25, 0.8, -2.0
rB = (1 - rA)/2

def B_sectors(a, nax, x, M, N):
    """B with axial flux: chiral fluxes m+nax, -m+nax. Direct per-m evaluation with the parity trick per chiral."""
    q = x*x
    z = np.exp(2j*np.pi*(np.arange(N)+0.5)/N)
    def chir_tab(fug, R, numax):
        out = {}
        for p in (0, 1):
            s = (-1.0)**p
            cn = s/fug*x**(2-R+p); cd = s*fug*x**(R+p)
            Pn, Pd = qinf(cn, q), qinf(cd, q)
            fn = np.ones(N, complex); fd = np.ones(N, complex)
            j, nu = 0, p
            while nu <= numax:
                out[nu] = (Pn/fn)/(Pd/fd)
                fn = fn*(1 - cn*q**j); fd = fd*(1 - cd*q**j)
                j += 1; nu += 2
        return out
    numax = M + abs(nax)
    TQ = chir_tab(z*a, rB, numax); TT = chir_tab(a/z, rB, numax)
    S = {}
    for m in range(-M, M+1):
        f1, f2 = m + nax, -m + nax
        zp = (-abs(f1) + abs(f2))//2; ap = -(abs(f1) + abs(f2))//2
        S[m] = a**ap * x**((1-rB)*(abs(f1)+abs(f2))/2) * np.mean(z**(zp + m) * TQ[abs(f1)] * TT[abs(f2)])
    return S

if __name__ == "__main__":
    for eps, N in [(0.2, 1024), (0.1, 2048), (0.05, 4096)]:
        x = np.exp(-eps/2); q = x*x; Lq = np.log(1-q)
        nA = 2*int(round((kap + 4*Lq)/(2*eps))); nax = nA//2
        y = q**(-1j*mu); wA = np.exp(1j*theta); a = np.exp(1j*theta/2)
        M = int(70/eps)
        ch = [(lambda z: z*y, 1, rA), (lambda z: z/y, 1, rA)]
        SA = sectors(ch, lambda m: -(m+abs(m)) - nA, x, M, N, lambda m: x**((1-rA)*abs(m)))
        IA = sum(wA**(-m)*SA[m] for m in SA)
        SB = B_sectors(a, nax, x, M, N)
        IB = sum((y*y)**(-m)*SB[m] for m in SB)
        contact = (wA * x**(3*rA - 1))**(nA/4)
        ab = {m: abs(SB[m]) for m in SB}; tot = sum(ab.values())
        mpk = max(ab, key=ab.get)
        mean = sum(abs(m)*v for m, v in ab.items())/tot
        fixed = abs(sum((y*y)**(-m)*SB[m] for m in SB if abs(m) <= 3))/abs(IB)
        print(f"eps={eps:<5} nA={nA:<5} nax={nax:<5} |I_A/(contact I_B)-1|={abs(IA/(contact*IB)-1):.1e}  "
              f"B: fixed-flux share {fixed:.3f}, eps<|m|> {eps*mean:.3f}, peak m={mpk} (eps m={eps*mpk:+.3f}); "
              f"u = x^|nax| / eps = {x**abs(nax)/eps:.3f}")
