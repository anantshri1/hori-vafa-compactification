"""LOG 47: exact mirror flux sectors at LARGE flux via Fourier convolution (no 4d angle grid).
Sector integrand = C * z^{ez} * prod_i F_i(zeta_i), zeta_i = (-1)^{nu_i} f_i z^{Q_i}, F_i = qpoch ratio (zero points and
CS/FI monomials are in ez, C exactly as in idx_jax). F_i(theta) = sum_k c_i(k) e^{ik theta}; the T^4 zero mode imposes
sum_i k_i Q_i = -ez: k = k0 + t1 (T1 cycle) + t2 (T2 cycle) + t3 (B,L-B,T-L,-T cycle), t in Z^3.
Cost per sector ~ (2K+1)^3 with K the Fourier cutoff (small in the Bessel regime), independent of |m|."""
import numpy as np
from fractions import Fraction as Fr
from idx import qpoch
import optc_defs as O
CH = [(O.one(O.B, 1), {'eta': -1}, O.r)] + O.rest          # B, L-B, T-L, -T, R-B, T-R, B-T
Q = np.array([c[0] for c in CH], int)
CYC = np.array([[0,1,1,0,0,0,1], [0,0,0,0,1,1,1], [1,1,1,1,0,0,0]])   # integer kernel of Q^T
assert np.all(CYC @ Q == 0)
BASIS = [0, 1, 3, 4]                                           # B, L-B, -T, R-B: unimodular block
Qb = Q[BASIS].T.astype(float); assert abs(abs(np.linalg.det(Qb)) - 1) < 1e-9
def bookkeeping(m, x, fl, n):
    K = O.Kmat(1); Kf = {(O.B, 'eta'): -.5}; FI = O.FIfull; bg = {'eta': -n}
    lf = {f: np.log(v) for f, v in fl.items()}
    ez = [Fr(0)]*4; ef = {f: Fr(0) for f in fl}; ex = Fr(0)
    for a in range(4): ez[a] += sum(Fr(K[a][b])*int(m[b]) for b in range(4))
    for (a, f), lev in Kf.items(): ef[f] += Fr(lev)*int(m[a]); ez[a] += Fr(lev)*bg.get(f, 0)
    for a, fs in FI.items():
        for f in fs: ef[f] -= int(m[a]); ez[a] -= bg.get(f, 0)
    nus = []
    for Qc, flc, R in CH:
        nu = int(sum(Qc[a]*m[a] for a in range(4)) + sum(c*bg.get(f, 0) for f, c in flc.items()))
        for a in range(4): ez[a] += Fr(-Qc[a]*abs(nu), 2)
        for f, c in flc.items(): ef[f] += Fr(-c*abs(nu), 2)
        ex += (Fr(R) - 1)*Fr(-abs(nu), 2); nus.append(nu)
    assert all(e.denominator == 1 for e in ez)
    logc = sum(float(ef[f])*lf[f] for f in fl) + float(ex)*np.log(x)
    return np.array([int(e) for e in ez]), np.array(nus), np.exp(logc), lf
def coeffs(nu, R, phase, x, Nf):
    q = x*x; th = 2*np.pi*np.arange(Nf)/Nf
    zeta = (1 if nu % 2 == 0 else -1)*phase*np.exp(1j*th); a = abs(nu)
    F = qpoch(x**(2 - R + a)/zeta, q)/qpoch(zeta*x**(R + a), q)
    return np.fft.fft(F)/Nf                                    # c(k) at index k mod Nf
def sector(m, x, fl, n, Kc=24, Nf=256):
    ez, nus, C, lf = bookkeeping(m, x, fl, n)
    ph = [np.exp(sum(c*lf[f] for f, c in flc.items())) for _, flc, _ in CH]   # flavour fugacity (complex) of each chiral
    c = [coeffs(nus[i], CH[i][2], ph[i], x, Nf) for i in range(7)]
    k0 = np.zeros(7, int); k0[BASIS] = np.round(np.linalg.solve(Qb, -ez)).astype(int)
    assert np.all(Q.T @ k0 == -ez)
    t = np.arange(-Kc, Kc + 1); T1, T2, T3 = np.meshgrid(t, t, t, indexing='ij')
    tot = np.ones_like(T1, dtype=complex)
    for i in range(7):
        k = k0[i] + CYC[0, i]*T1 + CYC[1, i]*T2 + CYC[2, i]*T3
        tot = tot*c[i][np.mod(k, Nf)]
    return C*tot.sum()

LADDER = [0, 1, 2, 3, 4, 5]          # all but the vertical B-T
SYM = True   # symmetric (pure-phase) LG approximant x^{1+nu}(zeta - 1/zeta): the 2d S^2 form; avoids spurious growth
def coeffs_lg(nu, R, phase, x, Nf):
    """Ladder-LG approximant of one chiral on the circle: exp[(zeta x^{R+nu} - x^{2-R+nu}/zeta)/(1-q)] with SIGNED nu
    (i.e. always the monomial w^Q, whatever the sign of nu)."""
    q = x*x; th = 2*np.pi*np.arange(Nf)/Nf
    zeta = (1 if nu % 2 == 0 else -1)*phase*np.exp(1j*th)
    with np.errstate(over='ignore'):
        F = np.exp((zeta*x**(R + nu) - x**(2 - R + nu)/zeta)/(1 - q)) if not SYM else np.exp(x**(1 + nu)*(zeta - 1/zeta)/(1 - q))
    return np.fft.fft(F)/Nf
def sector_lg(m, x, fl, n, Kc=16, Nf=256):
    """Same sector with the 6 ladder chirals replaced by their LG approximant (signed nu) and the vertical chiral by 1;
    zero points / CS / FI monomials and constants kept EXACT (identical to sector())."""
    ez, nus, C, lf = bookkeeping(m, x, fl, n)
    ph = [np.exp(sum(c*lf[f] for f, c in flc.items())) for _, flc, _ in CH]
    c = []
    for i in range(7):
        if i in LADDER:
            # signed nu: |w^Q| = x^{nu}; the exact factor used |nu| -- here we keep the ladder monomial
            c.append(coeffs_lg(nus[i], CH[i][2], ph[i], x, Nf))
        else:
            e = np.zeros(Nf, complex); e[0] = 1.0; c.append(e)
    k0 = np.zeros(7, int); k0[BASIS] = np.round(np.linalg.solve(Qb, -ez)).astype(int)
    t = np.arange(-Kc, Kc + 1); T1, T2, T3 = np.meshgrid(t, t, t, indexing='ij')
    tot = np.ones_like(T1, dtype=complex)
    for i in range(7):
        k = k0[i] + CYC[0, i]*T1 + CYC[1, i]*T2 + CYC[2, i]*T3
        tot = tot*c[i][np.mod(k, Nf)]
    return C*tot.sum()
