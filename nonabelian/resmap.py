"""LOG 66: residual map, eps = 0.3: exact - model (Bessel switch and branch) binned by (m_T, m_B) and by each diagonal's
regime (nu <= 0 / crossover 0 < nu < nu* / Bessel nu >= nu*). Progress + partial maps printed every 200 tasks."""
import sys, numpy as np, itertools, time, os, json
from multiprocessing import Pool
import fourier_sector as FS
eps, r = 0.3, 0.3; x = np.exp(-eps/2); ns = int(round(2*np.log(1/eps)/eps)); n = 4*ns; h, Tlo, Thi = 8, -48, 12
nst = 2*np.log(1/eps)/eps
FS.set_r(r)
y3 = np.exp(1j*np.array([0.4, 1.1, -0.6])); y = list(y3) + [1/np.prod(y3)]; e = np.exp(0.9j); rat = [y[i+1]/y[i] for i in range(3)]
fl = {'eta': 1/e, 'f1': rat[0], 'fT': rat[1]*x**(r/2), 'fB': rat[1]*x**(-r/2), 'f3': rat[2]}
def reg(nu): return 'inv' if nu <= 0 else ('cross' if nu < nst else 'bessel')
def task(t):
    mT, mL = t; bins = {}
    for mB, mR in itertools.product(range(-n - h, h + 1), repeat=2):
        m = (mL, mT, mB, mR)
        ex = FS.sector(m, x, fl, n, 12, 128)
        if abs(ex) < 1e-300: continue
        bs = FS.sector_bessel(m, x, fl, n, 12, 128, 1.0); br = FS.sector_branch(m, x, fl, n, 12, 128)
        _, nus, _, _ = FS.bookkeeping(m, x, fl, n)
        keys = [f"TB {mT//8*8},{mB//8*8}", "diag " + ','.join(reg(nus[i]) for i in (1, 2, 4, 5))]
        for k in keys:
            a = bins.setdefault(k, [0j, 0j, 0j]); a[0] += ex; a[1] += ex - bs; a[2] += ex - br
    return bins
if __name__ == "__main__":
    ncpu = len(os.sched_getaffinity(0)); tasks = [(mT, mL) for mT in range(Tlo, Thi + 1) for mL in range(-n - h, h + 1)]
    print(f"{len(tasks)} tasks on {ncpu} workers; nu* = {nst:.2f}", flush=True); tot = {}; t0 = time.time()
    def report(final=False):
        Ex = sum(v[0] for k, v in tot.items() if k.startswith('TB'))
        print(("FINAL" if final else "PARTIAL") + f" exact total {Ex:.6e}", flush=True)
        for pre in ('diag', 'TB'):
            items = sorted([(k, v) for k, v in tot.items() if k.startswith(pre)], key=lambda t: -abs(t[1][1]))[:12]
            for k, v in items:
                print(f"   {k:32s} exact {v[0]:.4e}  ex-Bessel {v[1]:.4e} ({abs(v[1])/abs(Ex):.3f} of |T|)  ex-branch {v[2]:.4e} ({abs(v[2])/abs(Ex):.3f})", flush=True)
    with Pool(ncpu) as pool:
        for i, b in enumerate(pool.imap_unordered(task, tasks)):
            for k, v in b.items():
                a = tot.setdefault(k, [0j, 0j, 0j]); a[0] += v[0]; a[1] += v[1]; a[2] += v[2]
            if (i + 1) % 400 == 0: print(f"  {i+1}/{len(tasks)} [{time.time()-t0:.0f}s]", flush=True); report()
    report(True)
