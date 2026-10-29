"""2d functions: CP^1 S^2 partition function (residues) and the B-side limiting integral."""
import mpmath as mp
mp.mp.dps = 30

def Z_CP1(r, mu, kap, th, M=None, K=None):
    """Z = sum_m e^{-i th m} int dsig/2pi e^{i kap sig} prod_{s=+-} G(r/2 - i(sig+s mu) + |m|/2)/G(1-r/2 + i(sig+s mu) + |m|/2)
    closed in the lower half plane (kap<0)."""
    r, mu, kap, th = map(mp.mpf, (r, mu, kap, th))
    M = M or int(60/abs(kap)) + 20; K = K or int(60/abs(kap)) + 20
    tot = mp.mpc(0)
    for m in range(-M, M+1):
        h = mp.mpf(abs(m))/2
        def zch(u):
            return mp.gamma(r/2 - 1j*u + h)/mp.gamma(1 - r/2 + 1j*u + h)
        for s in (+1, -1):
            for k in range(K):
                sig = -s*mu - 1j*(r/2 + h + k)
                u_other = sig - s*mu
                res = (-1)**k/mp.factorial(k)/mp.gamma(1 - r/2 + 1j*(sig + s*mu) + h)
                tot += mp.exp(-1j*th*m)*mp.exp(1j*kap*sig)*res*zch(u_other)
    return tot

def ZB_bessel(r, mu, kap, th):
    """4 e^{kap r/2} int_0^inf ds cos(4 mu s) J0(2 sqrt2 e^{kap/4} sqrt(cosh 2s - cos th))"""
    r, mu, kap, th = map(mp.mpf, (r, mu, kap, th))
    c = 2*mp.sqrt(2)*mp.exp(kap/4)
    f = lambda s: mp.cos(4*mu*s)*mp.besselj(0, c*mp.sqrt(mp.cosh(2*s) - mp.cos(th)))
    # J0 argument grows like e^s: oscillation; integrate piecewise between zeros-ish
    # break points: where c*e^s/sqrt2 ~ pi k  ->  use intervals in s of geometric growth
    pts = [mp.mpf(0)]
    s = mp.mpf(0)
    while s < 60:
        X = c*mp.sqrt(mp.cosh(2*s) - mp.cos(th))
        ds = min(mp.mpf(0.5), mp.pi/(2*max(X, 1)) , mp.pi/(8*max(abs(4*mu), mp.mpf(1e-3))))
        s += ds; pts.append(s)
    val = mp.quad(f, pts)
    return 4*mp.exp(kap*r/2)*val, pts[-1]

def ZB_rho(r, mu, kap, th):
    """Same integral written in rho = e^{kap/4+s}:
    Z_B = 4 e^{kap r/2} int_{sqrt b}^inf (drho/rho) cos(4 mu log(rho/sqrt b)) J0(2|rho - b e^{i th}/rho|),  b=e^{kap/2}"""
    r, mu, kap, th = map(mp.mpf, (r, mu, kap, th))
    b = mp.exp(kap/2); sb = mp.sqrt(b)
    f = lambda p: mp.cos(4*mu*mp.log(p/sb))*mp.besselj(0, 2*mp.sqrt(p**2 + b**2/p**2 - 2*b*mp.cos(th)))/p
    # finite part up to P0 with plain quadrature, tail with quadosc
    P0 = sb + 20
    head = mp.quad(f, mp.linspace(sb, P0, 81))
    tail = mp.quadosc(f, [P0, mp.inf], omega=2)
    return 4*mp.exp(kap*r/2)*(head + tail)
