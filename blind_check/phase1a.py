"""Phase 1(a): exact identity I_A = I_B at n=0."""
import numpy as np
from idx import u1_index
one = lambda m: 1
alt = lambda m: 1 if m % 2 == 0 else -1

def IA(r, a1, a2, th, lx, sign=one, M=150, N=1024):
    # Q1,Q2 charge +1, fugacities y_i = e^{i a_i}, R=r; omega=e^{i th}; k=+1
    ch = [(1, 1j*a1, 0, r), (1, 1j*a2, 0, r)]
    return u1_index(ch, 1, 1j*th, 0, lx, sign, M, N)

def IB(RQ, RQt, la, lb, lwB, lx, sign=one, M=150, N=1024):
    ch = [(1, la, 0, RQ), (-1, lb, 0, RQt)]
    return u1_index(ch, -1, lwB, 0, lx, sign, M, N)

def IB_map(r, a1, a2, th, lx, sign=one):
    # symmetric frame: R_Q=R_Qt=(1-r)/2, a=b=omega^{1/2}(y1y2)^{-1/4}, omega_B=y1/y2
    lab = 1j*th - 0.5j*(a1+a2)
    return IB((1-r)/2, (1-r)/2, lab/2, lab/2, 1j*(a1-a2), lx, sign)

def IB_asym(r, a1, a2, th, lx, t, sign=one):
    # equivalent gauge frame: z -> z x^t:  R_Q=(1-r)/2+t, R_Qt=(1-r)/2-t, omega_B -> omega_B x^{-t}
    lab = 1j*th - 0.5j*(a1+a2)
    return IB((1-r)/2+t, (1-r)/2-t, lab, 0, 1j*(a1-a2) - lab/2 - t*lx, lx, sign)

def IB_naive(r, a1, a2, th, lx, sign=one):
    # naive S^3 translation: R=1-r both, Q carries omega*(y1y2)^{-1/2}, Qt nothing, omega_B=y1/y2
    lab = 1j*th - 0.5j*(a1+a2)
    return IB(1-r, 1-r, lab, 0, 1j*(a1-a2), lx, sign)

if __name__ == "__main__":
    pts = [(0.3, 0.7, -1.9, 0.8, np.log(0.35)),
           (0.41, 2.3, 0.4, -2.2, np.log(0.5)),
           (0.17, -0.6, 1.1, 2.9, np.log(0.25)),
           (0.62, 1.3, 2.9, 0.37, np.log(0.42))]
    for (r, a1, a2, th, lx) in pts:
        A = IA(r, a1, a2, th, lx)
        B = IB_map(r, a1, a2, th, lx)
        Bw = IB((1-r)/2, (1-r)/2, (1j*th-0.5j*(a1+a2))/2, (1j*th-0.5j*(a1+a2))/2, 1j*(a2-a1), lx)  # Weyl image
        Bs = IB_asym(r, a1, a2, th, lx, 0.4*(1-r))
        Bn = IB_naive(r, a1, a2, th, lx)
        Bss = IB_map(r, a1, a2, th, lx, sign=alt)
        Ass = IA(r, a1, a2, th, lx, sign=alt)
        print(f"pt r={r} x={np.exp(lx):.3f}: I_A={A:.15f}")
        print(f"   |I_A-I_B(map)|={abs(A-B):.2e}  Weyl wB=y2/y1: {abs(A-Bw):.2e}  shifted frame t=0.4(1-r): {abs(A-Bs):.2e}")
        print(f"   naive map: I_B={Bn:.10f}  |diff|={abs(A-Bn):.2e}")
        print(f"   both signs (-1)^m: {abs(Ass-Bss):.2e};  A sign 1 vs B (-1)^m: {abs(A-Bss):.2e}")
