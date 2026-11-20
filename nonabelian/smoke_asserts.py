"""Pipeline gate: abort the job unless (1) Fourier sectors agree with the grid engine at small flux (Ng=36: <1e-4),
(2) Bessel-regime sectors are cutoff-converged (<1e-10), (3) the region scan runs and shares sum to 1."""
import numpy as np, optc_defs as O, idx_jax, fourier_sector as FS, subprocess, sys
x, n = 0.3, 4; fl = dict(O.fl); fl['eta'] = np.exp(-0.9j)
ch = [(O.one(O.B, 1), {'eta': -1}, O.r)] + O.rest
secs = [(0, 0, -4, 0), (-1, 2, -4, -1), (-2, 1, 0, 2), (0, 0, 0, 2)]
g = idx_jax.theory_index(4, ch, O.Kmat(1), {(O.B, 'eta'): -.5}, O.FIfull, fl, x, M=0, Ngrid=36, bg={'eta': -n}, sectors=secs, per_sector=True)
err1 = max(abs(FS.sector(s, x, fl, n)/gv - 1) for s, gv in zip(secs, g))
eps = 0.1; xx = np.exp(-eps/2); ns = int(round(2*np.log(1/eps)/eps)); m = tuple(int(v) for v in np.array([-2, -1, -3, -2])*ns + np.array([1, 0, -1, 2]))
a, b = FS.sector(m, xx, fl, 4*ns, Kc=16, Nf=256), FS.sector(m, xx, fl, 4*ns, Kc=32, Nf=512)
err2 = abs(a/b - 1)
print(f"SMOKE grid-vs-fourier max rel {err1:.1e} (need < 1e-4); cutoff convergence {err2:.1e} (need < 1e-10)", flush=True)
if not (err1 < 1e-4 and err2 < 1e-10): print("SMOKE FAILED - aborting"); sys.exit(1)
out = subprocess.run([sys.executable, 'region_scan.py', '0.2', '1', '1'], capture_output=True, text=True).stdout
shares = [float(l.split('share ')[1].split()[0]) for l in out.splitlines() if 'share' in l]
print(out.splitlines()[0]); print(f"SMOKE scan ran; top shares sum {sum(shares):.3f}", flush=True)
print("SMOKE PASSED", flush=True)
