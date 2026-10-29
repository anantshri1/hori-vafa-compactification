"""Phase 2: direct evaluation of I_B in the limit x=e^{-eps/2}, y=q^{-i mu}, omega=e^{i th}, n even large,
with kappa := -(eps*n + 4 ln eps)  (eps solved from n).
Uses the Phase-1(b) map. Two evaluators:
  brute(): original |nu| chirals from idx.py
  chain(): same integrand, signed-flux form, Pochhammer ratios updated m -> m+-2 by single factors."""
import numpy as np, sys, time
from scipy.optimize import brentq
from idx import chir, poch

def eps_from_n(n, kap):
    return brentq(lambda e: e*n + 4*np.log(e) + kap, 1e-8, 10)

def setup(r, mu, kap, th, n):
    eps = eps_from_n(n, kap)
    return eps, -eps/2, (1-r)/2

def brute(r, mu, kap, th, n, N, mlo, mhi):
    eps, lx, R = setup(r, mu, kap, th, n)
    a1, a2 = eps*mu, -eps*mu; la = 1j*th; lwB = 1j*(a1-a2) - la/2
    phi = 2*np.pi*np.arange(N)/N; lz = 1j*phi
    pref = np.exp(n*(0.5j*th - 0.25j*a1 + 0.75j*a2) + n*(3*r-1)/4*lx)
    c = {}
    for m in range(mlo, mhi+1):
        f = np.exp(-m*lwB)*np.exp((m+n/2)*lz)*chir(lz+la, m+n, R, lx)*chir(-lz, -m, R, lx)
        c[m] = f.mean()*pref
    return eps, c

def chain(r, mu, kap, th, n, N, tol=1e-17, verbose=False):
    eps, lx, R = setup(r, mu, kap, th, n)
    x = np.exp(lx); q = x*x
    phi = 2*np.pi*np.arange(N)/N; z = np.exp(1j*phi)
    w = np.exp(1j*th)
    C = (-1)**n*x**n/w; Cp = (-1)**n*x**n*w
    def Zm(m): return z*(-1)**m*x**(-float(m))
    def exact(m):
        Z = Zm(m)
        P2 = poch(Z*x**(2-R), q)/poch(np.conj(Z)*x**R, q)
        P1 = poch(C*x**(2-R)/Z, q)/poch(Cp*x**R/np.conj(Z), q)
        return P1, P2
    glob = x**(n*r)
    def cm(m, P1, P2):
        return glob*np.exp(-1j*eps*mu*(n+2*m))*np.mean(P1*P2)
    c = {}
    m0 = -n//2
    for start in (m0, m0+1):
        P1, P2 = exact(start); c[start] = cm(start, P1, P2)
        # upward
        p1, p2, m = P1.copy(), P2.copy(), start
        small = 0
        while True:
            Zo = Zm(m); Zn = Zm(m+2)
            p2 = p2*(1 - Zn*x**(2-R))/(1 - np.conj(Zn)*x**R)
            p1 = p1*(1 - Cp*x**R/np.conj(Zo))/(1 - C*x**(2-R)/Zo)
            m += 2; c[m] = cm(m, p1, p2)
            small = small+1 if abs(c[m]) < tol else 0
            if small > 20: break
        p1, p2, m = P1.copy(), P2.copy(), start
        small = 0
        while True:
            Zo = Zm(m); Zn = Zm(m-2)
            p2 = p2*(1 - np.conj(Zo)*x**R)/(1 - Zo*x**(2-R))
            p1 = p1*(1 - C*x**(2-R)/Zn)/(1 - Cp*x**R/np.conj(Zn))
            m -= 2; c[m] = cm(m, p1, p2)
            small = small+1 if abs(c[m]) < tol else 0
            if small > 20: break
    return eps, c

if __name__ == "__main__":
    p = (0.3, 0.25, -2.0, 0.8)
    n = 60; N = 4096
    t = time.time(); eps, cc = chain(*p, n, N); t1 = time.time()
    mlo, mhi = min(cc), max(cc)
    print('eps', eps, 'm-range', mlo, mhi, 'chain sum', sum(cc.values()), f'{t1-t:.1f}s')
    eps, cb = brute(*p, n, N, mlo, mhi)
    print('brute sum', sum(cb.values()), f'{time.time()-t1:.1f}s')
    print("max |c_chain - c_brute|", max(abs(cc.get(m,0)-cb[m]) for m in cb))
