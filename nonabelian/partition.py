"""LOG 54: smooth partition of unity over ALL flux space, by branch pattern of the 7 chirals.
piece_p = sum_m prod_i s_i(nu_i)^{[p_i=+]} (1 - s_i(nu_i))^{[p_i=-]} * S(m),  s(nu) = (1 + tanh(nu/w))/2.
Exact pieces and ladder-LG pieces; exact pieces must sum to the full exact total T (from the electric side).
Usage: python partition.py eps r x_override(or 0) n_override(or 0) hLRB hT w"""
import sys, numpy as np, itertools, time
from multiprocessing import Pool
import optc_defs as O, fourier_sector as FS, elec_fourier as EF
eps, r, xo, no, h, hT, w = float(sys.argv[1]), float(sys.argv[2]), float(sys.argv[3]), int(sys.argv[4]), int(sys.argv[5]), int(sys.argv[6]), float(sys.argv[7])
Tlo = int(sys.argv[8]) if len(sys.argv) > 8 else -hT
VARIANT = sys.argv[9] if len(sys.argv) > 9 else 'ladder'   # 'ends': diagonals ladder-LG, B and -T exact, vertical 1   # lower end of the T range (LOG 55: weight extends to m_T << -12)
x = xo if xo > 0 else np.exp(-eps/2)
ns = int(round(2*np.log(1/eps)/eps)); n = no if no > 0 else 4*ns
FS.set_r(r)
yang, wang = [0.4, 1.1, -0.6], 0.9
y3 = np.exp(1j*np.array(yang)); y = list(y3) + [1/np.prod(y3)]; e = np.exp(1j*wang); rat = [y[i+1]/y[i] for i in range(3)]
fl = {'eta': 1/e, 'f1': rat[0], 'fT': rat[1]*x**(r/2), 'fB': rat[1]*x**(-r/2), 'f3': rat[2]}
def slab(task):
    mT, mL = task; out = {}
    for mB, mR in itertools.product(range(-n - h, h + 1), repeat=2):
        m = (mL, mT, mB, mR)
        ex = FS.sector(m, x, fl, n, Kc=12, Nf=128)
        if abs(ex) < 1e-300: continue
        lg = FS.sector_bessel(m, x, fl, n, Kc=12, Nf=128, frac=1.0) if VARIANT == 'bessel' else FS.sector_branch(m, x, fl, n, Kc=12, Nf=128) if VARIANT == 'branch' else (FS.sector_lg_end(m, x, fl, n, Kc=12, Nf=128, ends_exact=True) if VARIANT == 'ends' else FS.sector_lg(m, x, fl, n, Kc=12, Nf=128))
        _, nus, _, _ = FS.bookkeeping(m, x, fl, n)
        s = 0.5*(1 + np.tanh(np.array(nus)/w))
        # patterns with non-negligible weight: each chiral contributes both branches where s in (1e-6, 1-1e-6)
        opts = [[('+', si)] if si > 1 - 1e-9 else ([('-', 1 - si)] if si < 1e-9 else [('+', si), ('-', 1 - si)]) for si in s]
        for combo in itertools.product(*opts):
            wt = np.prod([c[1] for c in combo]); p = ''.join(c[0] for c in combo)
            a = out.setdefault(p, [0j, 0j]); a[0] += wt*ex; a[1] += wt*lg
    return out
if __name__ == "__main__":
    t0 = time.time(); tot = {}
    import os, json
    ncpu = len(os.sched_getaffinity(0)); tasks = [(mT, mL) for mT in range(Tlo, hT + 1) for mL in range(-n - h, h + 1)]
    print(f"{len(tasks)} tasks on {ncpu} workers", flush=True)
    with Pool(ncpu) as pool:
        for i, d in enumerate(pool.imap_unordered(slab, tasks)):
            for p, v in d.items():
                a = tot.setdefault(p, [0j, 0j]); a[0] += v[0]; a[1] += v[1]
            if (i + 1) % 200 == 0 or i + 1 == len(tasks):
                Tx = sum(v[0] for v in tot.values()); Tm = sum(v[1] for v in tot.values())
                print(f"  {i+1}/{len(tasks)} tasks [{time.time()-t0:.0f}s]  running exact total {Tx:.6e}  running MODEL total {Tm:.6e}", flush=True)
                json.dump({p: [v[0].real, v[0].imag, v[1].real, v[1].imag] for p, v in tot.items()}, open('partition_ckpt.json', 'w'))
    Tex = sum(v[0] for v in tot.values()); Tlg = sum(v[1] for v in tot.values())
    pre = (x**(3*(1-r))*e*(y[0]*y[1])**2)**(n/2); Tel = EF.index(y, e, n, x, r, 40 if n > 10 else 10, Nf=1024)/pre
    print(f"eps={eps} x={x:.4f} r={r} n={n} box L,R,B in [{-n-h},{h}], T in [{Tlo},{hT}], w={w}")
    print(f"COVERAGE: sum of exact pieces {Tex:.8e}  vs  electric-side total {Tel:.8e}  rel {abs(Tex/Tel-1):.2e}")
    print(f"LG total (same lattice) {Tlg:.8e}   |LG/T - 1| = {abs(Tlg/Tel-1):.2e}")
    for p, v in sorted(tot.items(), key=lambda t: -abs(t[1][0]))[:14]:
        print(f"  piece {p}: exact {v[0]:.6e} (|.|/|T| {abs(v[0])/abs(Tel):.3f}) | LG {v[1]:.6e} | |ex-LG|/|T| {abs(v[0]-v[1])/abs(Tel):.3f}")
