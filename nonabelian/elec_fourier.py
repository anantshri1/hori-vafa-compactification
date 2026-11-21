"""LOG 53: electric U(2)_{(0,-2)}[4,0] index with topological background flux n, by 1D Fourier convolution per sector
(no 2d grid, no cancellation problem). Same bookkeeping as idx.theory_index (vector=2, weyl=2)."""
import numpy as np, itertools
from fractions import Fraction as Fr
from idx import qpoch
def index(y, e, n, x, r, M, Nf=512):
    q = x*x; lf = {f'y{j}': np.log(y[j]) for j in range(4)}; lw = np.log(e*x**(2*r))
    th = 2*np.pi*np.arange(Nf)/Nf; tot = 0j
    for m in itertools.product(range(-M, M + 1), repeat=2):
        ez = [Fr(-(m[0] + m[1]))]*2                                   # K = -J
        ez = [ez[0] - n, ez[1] - n]                                    # FI background: z_a^{-n}
        ew = -(m[0] + m[1]); ef = {f: Fr(0) for f in lf}; ex = Fr(0)
        c = []
        for a in range(2):
            F = np.ones(Nf, complex); nu = int(m[a])
            for j in range(4):
                ez[a] += Fr(-abs(nu), 2); ef[f'y{j}'] += Fr(-abs(nu), 2); ex += (Fr(1 - r).limit_denominator(1000) - 1)*Fr(-abs(nu), 2)
                zeta = (1 if nu % 2 == 0 else -1)*np.exp(1j*th)*y[j]
                R = 1 - r
                F *= qpoch(x**(2 - R + abs(nu))/zeta, q)/qpoch(zeta*x**(R + abs(nu)), q)
            c.append(np.fft.fft(F)/Nf)
        d = abs(m[0] - m[1]); s = (-1)**(m[0] - m[1]); ex += -d                 # two W factors, x^{-d/2} each
        wl = {0: 1 + x**(2*d), 1: -s*x**d, -1: -s*x**d}                         # coefficient of z1^l z2^{-l}
        assert all(v.denominator == 1 for v in ez)
        k1, k2 = -int(ez[0]), -int(ez[1])
        val = sum(w*c[0][(k1 - l) % Nf]*c[1][(k2 + l) % Nf] for l, w in wl.items())
        const = np.exp(sum(float(ef[f])*lf[f] for f in lf) + ew*lw)*x**float(ex)
        tot += const*val
    return tot/2
if __name__ == "__main__":
    x, r = 0.03, 0.45
    y3 = np.exp(1j*np.array([0.4, 1.1, -0.6])); y = list(y3) + [1/np.prod(y3)]; e = np.exp(0.9j)
    for n in (2, 4):
        for M in (6, 10):
            print(f"n={n} M={M}: electric (Fourier) = {index(y, e, n, x, r, M):.10e}", flush=True)
    print("grid reference: n=2 3.079700354e-03+5.62e-08j ; n=4 1.38772856e-06+5.52e-08j")
