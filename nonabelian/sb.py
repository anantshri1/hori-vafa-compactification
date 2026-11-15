"""Double sine s_b at b=1 (Q=2), thesis convention: chiral of R-charge r, charge +1 under U(1)_x -> s_b(iQ/2(1-r) - x).
Integral rep: log s_b(x) = -i int_0^inf dw/w [ sin(2xw)/(2 sinh(bw) sinh(w/b)) - x/w ],  |Im x| < Q/2.
Checked against the functional equations s_b(x - ib/2)/s_b(x + ib/2) = 2 cosh(pi b x) (convention-dependent sign tested in sb_test)."""
import mpmath as mp
mp.mp.dps = 20
b = mp.mpf(1); Q = b + 1/b
def log_sb(x):
    f = lambda w: (mp.sin(2*x*w)/(2*mp.sinh(b*w)*mp.sinh(w/b)) - x/w)/w
    return -1j*mp.quad(f, [0, 1, 5, 20, mp.inf])
def sb(x): return mp.e**log_sb(x)
if __name__ == "__main__":
    for x in [0.1+0.2j, -0.3+0.1j]:
        r1 = sb(x - 1j*b/2)/sb(x + 1j*b/2)
        print(x, mp.nstr(r1, 10), " 2cosh(pi b x) =", mp.nstr(2*mp.cosh(mp.pi*b*x), 10), " 1/that =", mp.nstr(1/(2*mp.cosh(mp.pi*b*x)), 10))
