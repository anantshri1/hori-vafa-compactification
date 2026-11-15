"""Diagnostic 4: resolve electric - mirror by topological charge (Fourier in w) at two x values; the charge and x-scaling of
the leading mismatch identify the culprit operator."""
import numpy as np, u2_defs as U
import sys
CFG = (sys.argv[1], int(sys.argv[2]), int(sys.argv[3]), int(sys.argv[4])) if len(sys.argv) > 4 else ('A', 1, 1, 1)
yv = U.pts[0][0]; nth = 8
for x in (0.04, 0.07):
    U.x = x; U.M_mir = 2; U.Ng_mir = 12
    th = 2*np.pi*np.arange(nth)/nth
    E = np.array([U.electric(yv, np.exp(1j*t)) for t in th]); Mv = np.array([U.mirror(yv, np.exp(1j*t), *CFG) for t in th])
    cE = np.fft.fft(E)/nth; cM = np.fft.fft(Mv)/nth   # coefficient of w^{-k} at index k
    for k in [0, 1, 2, nth-1, nth-2]:
        kk = k if k < nth//2 else k - nth
        print(f"x={x} w^{-kk:+d}: E {cE[k]:.6f}  M {cM[k]:.6f}  diff {abs(cE[k]-cM[k]):.2e}", flush=True)
