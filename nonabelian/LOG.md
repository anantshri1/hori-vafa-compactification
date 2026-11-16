# Non-abelian checks: U(2)_{(0,-2)} with [4,0] <-> planar abelian mirror (started 15 Nov 2026)

Reference: thesis (Unitary_Mirrors_1.tex): eqs u_n_chiral, mirpartfun, Figure planar_mirror_general and caption,
Table tab:SQCD_index (N=2, F=4, r=1/5).

## Results so far
1. Index engine (idx.py) checked: abelian CS pair (engine_check.py, 4e-6, truncation-limited) and Aharony duality
   U(2)_0 with 2 flavours <-> mesons + V_+- (engine_check_u2.py, 5e-5 and converging).
2. Electric U(2)_{(0,-2)} [4,0] vs thesis reference index (table_scan2.py): 2.4e-5 (= table truncation O(x^{21/5}))
   with CS matrix K = -J (from u_n_chiral), chirals R = 1-r, FI fugacity = eta x^{2r}. The x^{2r} (R-topological mixing)
   fixed at r=0.3 by requiring the dressed monopole at R = 4r = R(meson of mirror).
3. Mirror from the thesis figure (refined.py, mirror_table.py): nodes L, T, B, R; diagonals R=r, vertical R=2-2r;
   self CS +1, mixed -1/2 (diagonal), +1/2 (vertical), -1/2 for eta-B; orientation fixed by closed triangles and the
   meson path eta -> B -> L -> T -> flavour; eta-fugacity eta^{-1}; column-2 FI R-mixing x^{+-r/2} (T: +, B: -).
   - vs thesis table (r=1/5, x=0.05): 2.5e-5 (table truncation), c = r/2 unique (c=0, 0.2 give 1e-2).
   - vs electric directly, r=0.3 unrefined: 1.8e-5 at c = r/2 (c=0.1, 0.2, 0.35 give 5e-3 .. 2e-2).
   - vs electric, refined flavours, r=0.3, x=0.05: 1.7e-5, 7.7e-6 at two points (numerical floor: mirror M=2, Ng=20).
   - running: x=0.12, M=3 (refined_x012.log), to resolve the FI orientation s, which x=0.05 cannot see.

## Errors of mine caught (for the methods log)
- W-boson factor missing the (-1)^{m_a-m_b} that the (-1)^nu chiral convention implies (Aharony check off by 3%; fixed).
- Read the thesis S^3 formula (mirpartfun) literally: its right-triangle indices contain a non-existent z^{(3)}_2, and it
  omits mixed CS terms on two links (parity requires K_ab = 1/2 sum Q_a Q_b mod 1). The figure + caption are correct.
- Grid aliasing (Ng=12) produced a 3e-3 floor that I first read as a physics mismatch; Ng=20 gives 2e-5.
- Calibrating on N=1 alone did not fix the R-mixing of the multi-node columns (it has no column with 2 nodes).
4. Higgs limit, zero-flux sector (higgs_m00.py): 3d sector of U(2)_{(0,-2)} [4,0] -> Gr(2,4) GLSM sector (Vandermonde
   eps^2 sigma_12^2, Gamma ratios, (1-q)^{R-1-2i(sigma+mu)}) with |ratio-1| = 0.43, 0.11, 0.057 at eps = 0.2, 0.1, 0.05: O(eps).
   Next: nonzero flux sectors need a background topological flux scaled ~1/eps (as for CP^1) to give a finite 2d FI.
- x=0.12 refined run was killed when the cloud environment restarted; not needed (x=0.05 result stands).
5. Background fluxes added to idx.py (FI: z_a^{-n_f}; mixed gauge-flavour: z_a^{K n_f}; chirals: nu += charge*n_f).
   Engine check (bg_check.py): abelian CS pair I_A(n) = (w x^{3r-1})^{n/4} I_B(n_a = n/2) reproduced to 1e-10 for n = +-2, 4.
6. U(2) sector map with background topological flux (bg_u2.py, n = 2): eta-flux in {+-2, +-1} alone does NOT give a
   monomial ratio (best: eta-flux -2, residual 1e-2 after allowing a w-power). Missing: flavour/column background fluxes
   and/or an x^{...} R-mixing prefactor; the ratio test must also allow y-powers. Open.
7. Flux-2 sector map (bg_u2_fit.py, bg_u2_check.py): with mirror eta-flux -n, the ratio I_el(n=2)/I_mir is close to a
   monomial C w (y1 y2)^2 (powers read off: w 1.02, y1 1.99, y2 1.97, y3 -0.01; the y1 y2 factor matches the thesis
   background term e^{2 pi i eta (X_1 + X_2)}), C ~ x^{2.1}. But C differs by 0.9% between two points, unchanged by
   M = 2 -> 3. Either an extra ingredient (another background flux / R-mixing) or grid aliasing (Ng=20) with flux.
   Next: Ng = 24/28 to separate the two; then flux sectors of the Higgs limit.
8. The 1% residual of item 7 is NUMERICAL, on the electric side (bg_diag.py): with background flux n=2 the electric
   flux sum is unstable in the cutoff: (M,Ng) = (8,40) 1.5674e-3, (8,160) 1.5664e-3, (9,160) 1.5652e-3+2.4e-5i,
   (10,200) 1.05e-3-2.1e-4i. Large-flux sectors have exponentially large integrands whose gauge projection cancels
   catastrophically in double precision; grid refinement does not cure it. The mirror is stable (Ng 20 vs 26: 1e-8).
   Fix needed: exact projection (Laurent/series in x, or high precision) for the electric side before the map can be
   pinned beyond ~1e-3. The monomial C w (y1 y2)^2 with eta-flux -n stands as the working hypothesis.
9. Flux-2 sector map CONFIRMED in a stable regime (bg_u2_stable.py; r=0.45, x=0.03, electric M=6 Ng=96 where M=5,6 agree
   to 3e-5, mirror M=3): I_el(top flux 2) = x^{3(1-r)} w (y1 y2)^2 I_mir(eta-flux -2), C = 0.0030711, 0.0030719 at two
   points vs x^{3(1-r)} = 0.0030697 (agreement ~3e-4; the earlier r=0.3, x=0.05 value C ~ x^{2.1} = x^{3(1-r)} too).
   Item 8's diagnosis stands for M >= 7 (cancellation), but the 1% in item 7 was mostly flux truncation at r=0.3.
10. Higgs-limit flux sectors (higgs_flux.py): first attempt INCONCLUSIVE, method flawed. Compared the exact 3d sector
    with a 2d Gr(2,4) sector integrated on a finite sigma-window; ratio/ratio(0,0) is not sector-independent and does not
    converge (scan over kappa sign, theta sign, (-1)^{sum m}). Cause identified: with 4 chirals of R=0.7 per Cartan the
    2d integrand behaves like |sigma|^{4(R-1)+2} = |sigma|^{0.8} at large sigma (Vandermonde included), so the window
    integral is meaningless; the 2d side must be evaluated by residues (6 Higgs vacua of Gr(2,4), vortex sums), as was
    needed for CP^1. Next: residue-sum Gr(2,4) sectors; running coefficient is 8 (4 chirals x (1-q)^{-2 i sigma}), not 4.
11. gr24_res.py: Gr(2,4) sectors by residues implemented (6 ordered vacua pairs, vortex sums). Comparison with the 3d sectors
    still fails, and the reason is on the 3d side: higgs_flux.sector3d integrates sigma over a finite window [-7,7], but the
    3d contour is the full circle, sigma in [-pi/eps, pi/eps]. With background flux n ~ 500 the full-circle grid needs
    ~1e6 points with q-Pochhammers at q -> 1 (thousands of factors): too slow by brute force.
    Plan: evaluate the 3d sectors by residues as well (poles of the Pochhammers inside |z|=1 = K-theoretic vortex sums at
    the same 6 Higgs vacua) and compare residue by residue with the 2d vortex sums. This is cheap and is the clean form
    of the Higgs-limit statement.
12. Higgs-limit flux sectors by residues (res3d.py, res3d_scan.py). 3d sectors of U(2)_{(0,-2)}[4,0] with background top.
    flux n ~ 1/eps evaluated as K-theoretic vortex sums (poles sigma_a = -mu_{j_a} - i(R/2+|m_a|/2+k_a), closed below),
    Gr(2,4) sectors at the same poles; residues by 4-direction averaging, vectorised (8 s per eps).
    RESULT: residue by residue, for each distinct-flavour vacuum (j1,j2), ratio_3d/2d(m) / ratio(0,0) -> (-1)^{m1+m2}
    with O(eps) corrections. E.g. m=(1,0), vacuum (2,3): -0.534, -0.734, -0.857, -0.926 at eps = 0.2, 0.1, 0.05, 0.025;
    m=(2,0): 0.227, 0.479, 0.694, 0.834. Sign = theta -> theta + pi (the usual shift). Sector totals converge more slowly.
    Caveat: poles with j1 = j2 and k1 != k2 are NOT zero sector by sector (Vandermonde does not vanish there) and are large
    (|same-j sum| >> |total|), so they cancel against something; 3d and 2d same-j pieces do not yet agree pole by pole.
    Open: show the same-j pieces cancel in the sum over sectors (expected: vortex-sum form has distinct j only).
    Coulomb-limit ladder test: NOT started this check-in.
13. Same-flavour (j1=j2) residues cancel in the sum over flux sectors (samej.py), in 3d and 2d separately:
    |samej|/|distinct| summed over m in [-M,M]^2 = 1e-2, 1e-7, 1e-11 for M = 2, 4, 6 (eps = 0.1; eps = 0.05 same trend).
    So the LOG 12 caveat is resolved: per-sector same-j pieces are artefacts of the |m| pole labelling and drop out of the
    full partition function, leaving the distinct-flavour vortex sums, which match vacuum by vacuum. Higgs limit
    U(2)_{(0,-2)}[4,0] -> Gr(2,4) GLSM: established at the level of all flux sectors.
14. Coulomb limit of the planar mirror: CONTINUUM FOUND (coulomb_dir.py). Zero gauge flux, eta-flux -n (image of the
    electric Higgs-limit background). The B chiral carries |nu| = n, so its modulus is sigma-independent on the real
    contour (it is the one that becomes massive). Large-sigma falloff |f| ~ s^p along rays (eps=0.02, s=4..32):
      (1,1,1,1): p = -0.77..-1.0   (2d prediction -0.7: only the -T chiral charged)
      (1,0,0,0): -1.39 (pred -1.4) ; (0,1,0,0): -1.72 (pred -1.7) ; (0,0,1,0) and (1,0,1,1): -1.0..-1.28 (pred -1.0)
      generic (1,2,3,4): -3.4 (pred 5(r-1)+(1-2r) = -3.1)
    2d counting rule |Gamma(R/2-iQ.s)/Gamma(1-R/2+iQ.s)| ~ |Q.s|^{R-1} reproduced. Non-integrable cone: the 2-plane
    sigma_L = sigma_B = sigma_R (coords (a,b,a,a)); there only T-L, -T, T-R, B-T are charged, p = 3(r-1)+(1-2r) = -1.7 > -2.
    So the 2d mirror has a 2-dimensional continuum (non-compact Coulomb directions), lifted at most by oscillating FI/CS
    phases. Per agreement with Anant: STOPPED here to ask how to treat it.
15. Anant (16 Nov): try all three treatments of the continuum, (a) LG with continuum variables, (b) regulate and take the
    limit, (c) integrate out the massive B chiral first; evaluate each independently, also via the elliptic genus.
    Option (c), 3d step DONE (optc.py, x=0.05, r=0.3, M=3): mirror with eta-flux -n equals, up to O(x^n),
      (eta x^{1-r})^{n/2} I_red,  I_red = quiver without the B chiral and without the B-eta mixed CS, K_BB: 1 -> 1/2,
      extra FI fugacity x^{-(1-r)/2} on B (R-mixing). Check: |full/(red C)| = 1 - 3e-7, 1 - 8e-8 (n=6, 8); the phase was
      exactly e^{-i n theta}, i.e. a sign typo in the eta power of C in the script (fixed, docstring was right).
    Parity: K_BB = 1/2 with B charges (L-B, R-B, B-T) -> 1/2*3 = 3/2 = 1/2 mod 1, consistent.
    Consequence: the eta background does NOT leave an FI on B (chiral zero point and mixed CS cancel); the would-be
    continuum sigma_L = sigma_B = sigma_R must be controlled by the CS levels of the reduced quiver. Next: its 2d limit.
16. Power counting of the reduced quiver's 2d Coulomb integrand (cones.py -> cones.out), all subspaces V cut out by
    vanishing chiral weights; convergence needs exponent sum_{charged}(R_i - 1) < -dim V.
    CORRECTION to LOG 14: the full 4d integral is itself divergent for r > 0 (exponent 3r-4, margin 3r); I wrongly called
    the generic direction integrable (the measured -3.4 at (1,2,3,4) matches 3r-4 = -3.1, and -3.1 > -4).
    Pattern: every divergent subspace has margin r, 2r or 3r (cured by r < 0), EXCEPT four lines with exponent exactly -1
    for any r: (1,0,1,1), (1,0,1,0), (0,0,1,1), (0,0,1,0). On each, the charged set is two diagonals + the vertical B-T,
    i.e. a superpotential triangle: R-charges sum to 2, so (r-1)+(r-1)+(1-2r) = -1 identically. Log-divergent,
    r-independent: lifted only by real FI / twisted masses, not by R-charge continuation.
    Option (b) verdict so far: R-continuation to r < 0 regulates everything except the 4 triangle lines; those need a real
    FI (|y| != 1). Option (a) remark: Hori-Vafa dualising each 2d chiral (6 Y's) with 4 sigma constraints leaves 2 LG
    variables (3 with B), fewer than dim Gr(2,4) = 4 -- suggests the extra directions sit at sigma ~ 1/eps (bulk divergence).
17. (check-in #9) Where the Gr(2,4) physics must sit on the mirror. Argument from LOG 15: at fixed gauge flux and fixed
    sigma, the reduced quiver (option c) carries NO eta dependence beyond the overall (eta x^{1-r})^{n/2}, yet the Gr(2,4)
    S^2 partition function depends non-trivially on its FI (kappa, theta) = image of eta. So the fixed-sigma, small-flux
    Coulomb limit cannot reproduce Gr(2,4); the identity full = C * red (LOG 15) used mirror M=3 < n and so excluded
    the sectors where it fails. The relevant region is m_B ~ -n, where nu_B = m_B + n is small: there the B chiral is
    LIGHT and the three chirals touching B in triangles (L-B, R-B, B-T) carry |nu| ~ n, i.e. become HEAVY. Integrating
    those out instead (option c'), the light theory is: node B with the eta-charged chiral; nodes L,T,R with T-L, -T, T-R;
    plus induced CS/FI from three heavy chirals (k shifts of +-1/2 and n-dependent FIs on L, T, R). The superpotential
    triangles are broken in this region, which is also where the r-independent log lines of LOG 16 came from.
    This also explains the bulk divergence of LOG 16: the 2d integrand at fixed flux has its weight at |sigma| ~ 1/eps.
    Next (concrete): evaluate the 3d mirror index restricted to sectors m_B = -n + m'_B, |m'| <= M, by shifting
    backgrounds (B gets FI-type background, L,T,R get mixed-CS-induced ones), compare (i) its eta dependence and (ii) its
    eps -> 0 limit with the Gr(2,4) sectors. 4 nodes x M~6-8 x Ng^4 is the first job where a GPU (Anant's HF credits) helps.
    Elliptic genus: not started; only meaningful once the light theory is identified.
