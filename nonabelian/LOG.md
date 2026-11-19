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
18. (check-in #10) Shifted flux sectors (idx.py gained mshift=; shifted.py, slices.py). Mirror with eta-flux -n:
    - Boxes centred at m_B = -n contribute as much as, or more than, the box at 0: |shift/box0| = 0.83 (x=0.05),
      1.49 (x=0.3) for n=4. So LOG 15's "full = C * red" used a truncation (M=3 < n) that misses O(1) of the index.
    - Beyond the monomial C = (eta x^{1-r})^{n/2}, eta dependence of each box is x-suppressed (4e-5 at x=0.05,
      ~4e-3 at x=0.3); with bg=0 it is 1e-2 at x=0.05. Must be re-examined as x -> 1.
    - Profile along the B-flux line (n=8, x=0.3, others |m|<=1, slices_n8_x03.log): |slice/C| =
      0.04 (m_B=+2), 0.12, 0.58 (0), 3.6, 5.1, 6.8, 8.9, 11.6, 15.2, 19.7, 23.1 (m_B=-8 = -n), 3.3, 0.04 (m_B=-10).
      Support is the interval m_B in [-n, 0] with sharp edges and a rotating phase across it.
    Interpretation: the 3d precursor of the 2d continuum is a LINE of flux sectors of length n ~ 1/eps; with
    u = eps m_B it becomes a COMPACT variable u in [-eps n, 0] (finite interval), carrying an oscillating weight --
    i.e. a 2d twisted-chiral direction on a finite interval rather than a non-compact Coulomb branch. Caveat: L, T, R
    fluxes truncated at |m| <= 1; their couplings to m_B (mixed CS -1/2, +1/2) probably require shifts ~ m_B/2 too.
    Next: (1) widen L,T,R windows co-moving with m_B (shift m_L, m_R ~ m_B/2?); (2) n = 8, 12, 16 to check the profile
    collapses as a function of m_B/n; (3) fit the slice phase per step against eta/theta and y. GPU worthwhile for (1)-(2).
19. GPU scan (HF Jobs T4, log hfjob_scan.log; JAX engine idx_jax.py agrees with idx.py to 5e-16 and the GPU reproduces
    the CPU slices digit for digit). Mirror, eta-flux -n, x = 0.3, slices along m_B:
    - ROBUST: support is m_B in [-n, 0] for n = 8, 12, 16, sharp edges (drop x10 at -n-1, x100 at -n-2).
    - ROBUST: slices near m_B ~ 0 (|m_B| <= 3) converged in the L,T,R window (identical for M = 3, 4, c = 0, 0.5).
    - NOT CONVERGED: the interior. Totals n=8: c=0 M=3 63-17i, M=4 77-16i; c=0.5 M=3 -1-78i, M=4 -16-12i.
      The L,T,R fluxes spread with |m_B|. A CPU diagnostic (n=4, m_B=-2) puts m_L, m_R near m_B and m_T near 0,
      i.e. co-moving coefficients ~ (1, 0, 1), not (0.5, 0.5, 0.5) as I guessed.
    - Profile at fixed x is a function of m_B (n only sets the edge), not of m_B/n.
    DESIGN FLAW of this batch (mine): the m_B/n collapse is only meaningful in the joint limit x = e^{-eps/2} -> 1 with
    n ~ 1/eps; at fixed x = 0.3 (eps ~ 2.4) it tests nothing about 2d. Next batch must scale x with n.
    Next GPU job (run_diag.sh): full (m_L,m_T,m_R) box at fixed m_B to locate the weight, at x = 0.3 and x = 0.78.
20. Flux-window diagnostic (HF job, 16 Nov; second run hit the 60m timeout: per-sector jit calls were unbatched, fixed
    with idx_jax per_sector=True). n=8, x=0.3, m_B=-4, full box [-9,5]^3, slice/C = 29.1-16.9i:
    marginal in m_L (= m_R exactly, L<->R symmetry): spread over [-5, 0], peak at -3; falls 1e-2 at -6 and at +2.
    marginal in m_T: peaked at 0, spread over [-3, +1].
    Pattern: m_B <= m_L, m_R <= 0 (approximately), m_T ~ 0: an INTERLACING pattern of flux intervals, nested in the
    interval m_B in [-n, 0]. This is Gelfand-Tsetlin-like (column 1: B; column 2: L, R; top: T / flavour), i.e. the
    shape of a ladder diagram -- suggestive of the Givental/BCFKvS ladder for Gr(2,4). Speculative; check at x -> 1.
21. Diagnostic at x = 0.78 (hfjob_diag2.log; n=8, Ng=16): the boxes did NOT contain the weight. m_T marginal peaks at the
    upper box edge (m_B=-4: 23.8 at +2 of box [-8,3]; m_B=-7: 86.6 at +2 of [-10,2]); L, R marginals barely decay at
    either edge (4.4 at +3, 8.3 at -10). L, R bulk still inside [m_B, 0] (peaks -1..-2 and -3) but with long tails.
    Interlacing (LOG 20) is therefore NOT confirmed at x=0.78: the spread grows as x -> 1 in every direction.
    Reading: all four fluxes spread ~ 1/eps, so the 2d continuum is 4-dimensional (u_a = eps m_a, a = L,T,B,R),
    consistent with LOG 16 (whole 4d Coulomb integral divergent) -- and 4 = dim_C Gr(2,4), which would resolve the
    LOG 16 count (2 LG variables from Hori-Vafa with sigma constraints, vs 4 needed). Hypothesis, not shown.
    Next job: asymmetric box L,R in [-14,6], T in [-4,14] at m_B=-4 (run_diag.sh).
22. Wide asymmetric box at x=0.78 (hfjob_diag3.log; n=8, m_B=-4, L,R in [-14,6], T in [-4,14], Ng=16):
    - L, R: sharp LOWER edge near -n (1.1 at -8, 0.12 at -9, 7e-3 at -10, 1.6e-8 at -14); bulk on [-7, 0];
      slow UPPER tail (4.7 at +6). Upper tail is real: single sectors (4,2,-4,-1), (6,2,-4,-1) identical at Ng=16,20,24.
    - T: bump on [-4, 6], peak at +2..+3, 1e-5 by +9; the RISE at m_T >= 11 (16.6 at +14) is GRID ALIASING:
      sector (-1,14,-4,-1) = 5.8e-2 (Ng=16), 2.3e-9 (Ng=20), 5.1e-15 (Ng=24). So the Ng=16 slice totals at x=0.78 are
      contaminated; any production run needs Ng >= 20-24 (rule: Ng > 2 max|gauge power| in the box).
    Picture so far: B in [-n,0], L,R bounded below by -n with a decaying upper tail, T in a fixed window around +2.
    Weight occupies a bounded region in all four flux directions whose size grows with n: a polytope-like support.
    Whether its shape is the Gelfand-Tsetlin polytope of Gr(2,4) is the question; not established.
23. Aliasing-safe polytope scan (hfjob_scan2.log; n=8, x=0.78, Ng=20, L,R in [-9,6], T in [-4,7]; 2.2 sectors/s on T4,
    ~23.5 min per slice, ~$0.80 total). m_B=-4 vs the Ng=16 run: slice -6.56+0.55i vs -6.52+1.28i (aliasing ~10%).
    Lower edge of the m_L (= m_R) marginal, taken where it drops to 1% of its peak:
      m_B =  0: ~ -5.5 | m_B = -2: ~ -7 | m_B = -4: ~ -9 | m_B = -6, -8: below the box (truncated).
    -> the L,R lower edge TRACKS m_B with slope ~ 0.9-1, with a constant offset ~ -5 (smearing at x=0.78, eps ~ 0.5).
    Upper side of L,R: slow tail beyond +6 in every slice, roughly m_B-independent (truncated by the box).
    m_T: same window [-4,7], peak at +2, in every slice: independent of m_B.
    Pattern: m_B (in [-n,0]) <= m_L, m_R <= (fixed upper scale), m_T fixed. This is qualitatively the interlacing of a
    Gelfand-Tsetlin pattern (bottom entry bounds the next row from below; the top row is fixed). Not yet quantitative:
    edges are smeared over ~5 units at x = 0.78, and the upper L,R edge is not contained.
    Slice totals for m_B = -6, -8 are unreliable (box truncation).
24. Both L,R edges contained (hfjob_scan3.log; n=8, x=0.78, Ng=20, L,R in [-14,12], T in [-5,8]; 78 min/slice, ~$1.6).
    Edges of the m_L (= m_R) marginal at 1% of peak:
                lower      upper     | T peak, T upper(1%)
      m_B =  0: -5.3       +8.5      | +2, ~+6
      m_B = -4: -9.2       +8.3      | +2, ~+6
      m_B = -8: -12.7      +8.0      | +2, ~+6
    Lower edge slope d(edge)/d(m_B) = 0.93; upper edge independent of m_B (drift 0.5 over 8 units); T window fixed.
    => at x = 0.78: m_B <= m_L, m_R <= const (~ T scale), smeared by ~5 units on each side: the interlacing
    inequalities of a Gelfand-Tsetlin pattern, with B the bottom entry and the T/flavour data as the fixed top row.
    Slice totals with contained boxes: m_B=0 -0.939+0.169i, -4 -6.18+0.58i (was -6.56+0.55i in [-9,6]),
    -8 33.5-3.1i (was 37.0+3.8i truncated).
    Open: (i) do the smeared edges sharpen as x -> 1 with n ~ 1/eps (must, if this is a 2d polytope)? (ii) is the
    upper bound of L,R really set by T (vary the T-related fugacity / R-mixing), (iii) match to the Gr(2,4) GT polytope
    inequalities and to the BCFKvS ladder potential.
25. Edge-sharpening test (hfjob_sharp.log, sharp_edges.txt; 2d (m_L, m_R) grids at m_T in {0,2,4}; ~4.3 h T4, ~$1.8).
    Aliasing check: corner sectors at Ng and Ng+8 all < 1e-9 of the bulk -> grids fine.
    m_L marginal, edges at 1% of peak, edge width = distance between the 10% and 0.01% levels:
      (n,x,eps)        m_B   lo edge (offset)   hi edge   width lo/hi   width/n lo/hi
      (8, .78, .497)    0    -5.17 (-5.17)      +8.23     2.5 / 3.2     0.32 / 0.40
      (8, .78, .497)   -8   -12.68 (-4.68)      +7.82     2.7 / 3.6     0.34 / 0.45
      (12,.85, .325)    0    -7.06 (-7.06)     +10.91     3.1 / 3.1     0.26 / 0.26
      (12,.85, .325)  -12   -18.30 (-6.30)     +10.14     6.4 / 6.3     0.53 / 0.53
    - At m_B = 0 the edges sharpen relative to n (0.32 -> 0.26 lo, 0.40 -> 0.26 hi).
    - At m_B = -n they do NOT: width/n grows 0.34 -> 0.53. That slice is strongly oscillatory at n=12 (dips, e.g. 55 at
      m_L=-14 between 293 and 283), which inflates the 10% crossing; but the growth is not obviously an artefact.
    - Edge OFFSETS in u = eps m: lo 2.57 -> 2.29 (m_B=0), 2.33 -> 2.05 (m_B=-n); hi 4.09 -> 3.55. A smearing that
      vanishes in 2d would shrink these by eps ratio 0.65; observed ~0.87-0.89. So most of the offset is a genuine finite
      feature in u, i.e. in 2d the L,R support would be ~ [u_B - 2.3, ~3.5], not the sharp GT inequality u_B <= u_L.
    VERDICT: INCONCLUSIVE, leaning against a literal Gelfand-Tsetlin polytope with edges at u_L = u_B. The qualitative
    structure (lower edge tracks m_B, upper edge fixed, T window fixed) survives at n=12. Two points are not a trend;
    further GPU spending is not the best next step. Next: derive the large-flux asymptotics analytically (dilogarithm /
    Bessel saddle in u for each chiral, as for the abelian pair) and predict the edge positions, then compare.
26. Analytic step 1-2 (amoeba.py, amoeba_par.py; am_*.log). Leading large-flux asymptotics: at u = eps|nu| fixed each chiral
    is, on the real torus, the pure phase exp(2i Im Li2(s a e^{i theta})/eps), a = x^{|nu|+1}; CS/FI/zero-point monomials are
    phases e^{i ez.phi}, ez ~ 1/eps. A sector is O(1) iff the phase is stationary somewhere on T^4:
        sum_i Q_ia [-2 log|1 - s_i a_i e^{i theta_i}|] + eps ez_a = 0   (4 real equations, 4 angles),
    otherwise exponentially small. Prediction: flux support = solvability region (amoeba-type: |1 - a e^{i theta}| ranges
    over [1-a, 1+a]). Compared with the measured m_L marginals (10%-of-peak edges, hfjob_sharp.log):
                         predicted (sharp)   measured (10%)
      n=8,  x=.78, m_B=0:    [-3, 4]         [-4.1, 6.9]
      n=8,  x=.78, m_B=-8:   [-10, 4]        [-11.7, 6.1]
      n=12, x=.85, m_B=0:    [-4, 8]         [-5.9, 9.6]
      n=12, x=.85, m_B=-12:  [-15, 7]        [-17.0, 8.7]
    - Upper edge predicted m_B-independent (4, 4; 8, 7): matches. Lower edge slope predicted 0.88 (n=8), 0.92 (n=12);
      measured 0.95, 0.93. Measured edges always 1-3 units OUTSIDE the predicted ones, as expected for stationary-phase
      (Airy-type) tails at a 10% threshold. Not a fit: no free parameters.
    - Consequence: the finite edge offsets in u (LOG 25) are the boundary of this stationarity region, a genuine 2d object
      (amoeba of the effective mirror's 1 - a_i zeta_i factors, i.e. |e^{-Y_i}| = e^{-u_i/2} in Hori-Vafa variables), not
      a smeared GT inequality u_L >= u_B. Interlacing-like shape (lower edge tracks B, upper fixed) comes out of it.
    Next (cheap, CPU only): (i) the region at fixed u as eps -> 0 (does it converge?); (ii) write the stationarity equations
    as Re of the critical equations of W = sum Li2 + CS + FI, take the 2d limit, compare with the BCFKvS ladder potential.
27. (i) eps -> 0 at fixed n*eps = 4 (amoeba_lim.py, am_lim.log; u-grid 0.5, T/R sampled on u-grids):
      predicted u_L support, m_B = 0:   n=8 [-1.5, 2.0], n=32 [-1.5, 2.5], n=64 [-1.5, 2.5]
                             m_B = -n:  n=8 [<=-4 (grid), 2.0], n=32 [<=-4 (grid)], n=64 [-5.0, 2.5]
    -> the region CONVERGES in u (to grid resolution): a fixed 2d region with finite edge offsets (lower edge ~ u_B - 1..1.5).
    (ii) Structure, derived: within a region of fixed signs of nu_i, X_i = x^{|nu_i|} zeta_i^{+-1} is a monomial in
    w_a = x^{m_a} z_a (times flavour), and the stationarity equations are the phase (Im) halves of dW = 0 for
         W = sum_i (+-) Li2(X_i) + (1/2) sum K_ab log w_a log w_b + (FI + background) . log w,
    i.e. the 3d twisted superpotential of the mirror on R^2 x S^1 (Bethe / K-theoretic form). So at fixed n*eps the scans sit
    in a "3d on a circle" regime: Li2 stays Li2, NOT the 2d Hori-Vafa form.
    CORRECTION to my framing in LOG 18-25: the genuine 2d limit (as for the abelian pair, LOG/project log step 1c) is the
    Bessel regime x^{|nu|} ~ eps, i.e. n ~ ln(1/eps)/eps with kappa = -eps n - c ln eps fixed. There Li2(X)/eps -> X/eps
    = e^{-Y} with Y = u/2 + ln eps: a Hori-Vafa LG with one exponential per light chiral (6, or 7 with B), constrained
    by the sigma (gauge) equations. Only THAT LG can be compared with the BCFKvS/Givental ladder of Gr(2,4).
    Next: (a) derive this LG explicitly (variables, constraints, CS-induced terms) and count independent variables;
    (b) fetch the exact BCFKvS / Marsh-Rietsch Gr(2,4) ladder potential and compare term by term; (c) if useful, one GPU
    slice in the Bessel regime (n ~ ln(1/eps)/eps) to confirm the weight moves there.
28. Analytic step 4 (ladder_cmp.py). Bessel-regime LG of the planar mirror: variables = the 4 gauge variables w_a (as for
    the abelian pair, where one node gave W = Y + Lambda/Y), one monomial w^{Q_i} per chiral:
      B: w_B | L-B: w_L/w_B | R-B: w_R/w_B | T-L: w_T/w_L | T-R: w_T/w_R | -T: 1/w_T | vertical B-T: w_B/w_T
    Gr(2,4) ladder (EHX hep-th/9605225 App. B; BCFKvS): W = X + Y/X + Z/X + W/Y + W/Z + q/W, 4 variables, 4 critical
    points (EHX: Y = Z = X^2/2, W = X^3/4, X^4 = 4e^t, Hessian 8e^t -- checked here, c=0 gives 4).
    MATCH: X = w_B, Y = w_L, Z = w_R, W = w_T, q <-> FI on T (eta/flavour). The six non-vertical chirals give EXACTLY the
    six ladder monomials, and the identification is the ladder's own shape (B bottom, L/R middle, T top, q at the end).
    EXTRA: the vertical chiral B-T gives a 7th monomial c X/W. With it, W_c has 8 critical points for generic c
    (c = 0.7e^{0.3i}, +-1), vs 4 for the bare ladder and chi(Gr(2,4)) = 6 (Marsh-Rietsch's compactified mirror).
    Possible readings (NOT decided): (a) 8 = 6 + 2, the 2 extra being topological vacua of the U(1) CS level -2, which
    ARW say decouple as massive; (b) the vertical term is absent/suppressed in the Bessel regime (it has R = 2-2r and sits
    in the superpotential triangles); (c) its coefficient is fixed by the reduction to a special value. The coefficient c
    and which chirals are simultaneously in the Bessel regime are not yet derived.
29. Anant (18 Nov): bets on reading (a) (CS level -2 leaves 2 extra vacua), cannot rule out (b),(c); do (a) first.
    While setting up (a), a derivation that bears on (b) (algebra, not yet checked numerically):
    each chiral's ratio of q-Pochhammers contains BOTH X_+ = x^{nu} zeta = w^{Q} and X_- = x^{nu}/zeta = x^{2nu} w^{-Q};
    in the Bessel regime x^{|nu|} ~ eps only the one with the small x-power survives at O(1) after the 1/eps:
    w^{Q} for nu > 0, w^{-Q} for nu < 0, with coefficient ~ x^{|nu|}/eps.
    In a superpotential triangle (e.g. L-B, T-L, B-T) the three fluxes sum to zero, so one of them has the opposite sign
    and |nu| = sum of the other two: x^{|nu|} ~ eps^2, its monomial is O(eps) and DROPS. Which one drops depends on the
    flux region. Region "vertical drops in both triangles" gives exactly the 6-term ladder W_tor (LOG 28). Other regions
    give other Laurent charts (the vertical B-T term present, one diagonal absent). Conjecture: these are cluster-type charts
    of one mirror, and the union is the compactified (Marsh-Rietsch) mirror with chi = 6 vacua; the CS-level-2 topological
    vacua then sit outside the Bessel regime (decoupled), as in ARW. This merges (a) and (b).
    Next: enumerate the sign regions/charts, check they glue by cluster mutations, count vacua per chart and of the union,
    and compute the electric 3d Bethe-vacuum count to see where the +2 lives.
30. Charts and vacuum counts (charts.py, count_crit.py).
    - Crude chart model (fluxes in units of nu*, chiral survives iff |k.Q + b| = 1) is too coarse: it misses the ladder
      chart unless the eta background is b = 4 (then all six ladder chirals survive and B-T drops, as LOG 29 says),
      and its other charts have <= 6 terms with 0 not interior. Kouchnirenko volumes are only upper bounds here.
    - KEY numerical fact: the vacuum count of the 6-term ladder depends on its coefficient invariant
      rho = c_LB c_TL / (c_RB c_TR) (the two B -> T paths, via L and via R):
         ladder, rho = 1 (EHX/BCFKvS normalisation):   4 critical points (2 at infinity, EHX)
         ladder, generic coefficients:                  8 critical points
         ladder + vertical X/W, generic:                12
      In the reduction rho is NOT 1 in general: the two paths carry different column FIs = flavour fugacity ratios of the
      electric theory, i.e. rho ~ twisted-mass data. So with generic electric masses the Bessel-regime ladder chart has
      8 = 6 + 2 vacua: chi(Gr(2,4)) = 6 equivariant vacua plus 2, matching reading (a) (CS-level-2 topological vacua),
      with the vertical term dropped by the triangle rule (reading (b)-mechanism). Massless limit rho -> 1: 2 vacua
      run to infinity (the familiar EHX deficit 4 = 6 - 2 of the bare ladder) -- to be checked: which 2?
    Next: (1) electric 3d Witten index / Bethe-vacuum count of U(2)_{(0,-2)}[4,0] (expect 8 if (a) holds);
    (2) derive rho and q from the reduction (FI/fugacity map) and check the 8 critical values against Gr(2,4) equivariant
    quantum cohomology (6) + 2 topological; (3) as rho -> 1 identify which vacua go to infinity.
31. Step 6 (bethe_exact.py, count_mass.py, values.py).
    (1) Electric Bethe vacua by exact elimination (resultants, random rational fugacities): symmetric quantisation with
        chiral factor X^{+1/2}/(1-X) gives 6 vacua (resultant degree 17, 12 ordered solutions z1 != z2); the other sign
        X^{-1/2}/(1-X) gives 20. The two are different real-mass chambers; the Higgs/GLSM chamber is the one with 6 =
        chi(Gr(2,4)): NO extra 2 vacua there. (Random-Newton counts for the 4-node mirror were not reliable: 12, 28 and an
        electric 207/33 from a mis-normalised rerun; discarded.)
    (2) RETRACTION of LOG 30's reading: twisted masses do not enter the Bessel-regime LG through the coefficient ratio rho
        but as LINEAR terms mu.log w (column FIs ~ mass differences, as in Hori-Vafa). Counting with that:
           ladder, unit coefficients, mu = 0:            4  (EHX)
           ladder, unit coefficients, mu generic:        6  = chi(Gr(2,4))   (two choices of mu, both 6)
           ladder + vertical X/W, mu generic:            10
        So: ladder + masses = 6 vacua = electric Higgs-chamber Bethe count. The vertical term must be absent (10 != 6),
        which is what the triangle rule (LOG 29) gives. Reading (a) (+2 CS vacua) is NOT supported in this chamber; the
        2 "extra" of LOG 30 were an artefact of putting masses in the coefficients. Reading (b) is.
    (3) Critical values: first attempt to match q dW/dq spectra of ladder+mu.log w with Gr(2,4) GLSM (-(sigma_i+sigma_j))
        over permutations/signs of the mass map FAILED (best error ~5): normalisation of mu vs m and of q not fixed yet.
        Open.
32. Step 7: VACUUM-LEVEL MATCH (pairsum.py, recover.py). Normalisation-free test first: Gr(2,4) equivariant QC has
    q dW/dq = -(r_i + r_j) at its 6 vacua, r = roots of prod(sigma - m_j) = q, so the 6 values must split into 3
    complementary pairs with equal sums.
    - Ladder + mu.log w with 4 INDEPENDENT mu: fails (rel. err 0.6-0.7). With the column constraint mu_B = mu_T (B and T
      sit in the same column of the planar quiver and carry the same column FI y3/y2): holds to 2e-13 (mu_B = -mu_T or
      only one of them: fails, 0.2-1.25).
    - Recovering the 4 numbers from the pair sums and fitting masses: with m = cumulative sums of the column FIs in the
      order (mu_L, mu_{B=T}, mu_R) (overall sign +-1, shift fixed by e1), the 4 numbers satisfy
          e1, e2, e3 = those of (m_1..m_4)  to 2e-13 - 6e-13,   and   e4(x) - e4(m) = q exactly (1.1974+0.5062i).
      i.e. they ARE the roots of prod_j (sigma - m_j) = q, and the 6 LG vacua are exactly the 6 pairs: the equivariant
      quantum cohomology vacua of Gr(2,4), vacuum by vacuum, 3 random mass points.
    (Some solver runs found 5 of 6 critical points; only complete runs used.)
    CONCLUSION (chain): 3d U(2)_{(0,-2)}[4,0] <-> planar abelian mirror (exact index identity) -> Higgs limit: Gr(2,4)
    GLSM (all flux sectors) ; Coulomb (Bessel) limit of the mirror: one monomial per light chiral, triangle rule drops the
    vertical chiral, column FIs -> twisted masses as log terms with mu_B = mu_T: the BCFKvS/Givental ladder LG with masses,
    whose 6 vacua are exactly Gr(2,4)'s. Reading (b); no extra CS vacua in this chamber.
    NOT yet done: derive the Bessel-regime LG directly from exact 3d sectors numerically (the monomial rule and triangle rule
    are derived analytically, checked only through this vacuum structure); coefficients/q from the map; mirror 3d Bethe
    count by exact elimination; the -/+ chamber with 20 vacua.
33. Coefficient step (1): R-MIXING CANCELS EXACTLY (rmix.py). In the Bessel regime the only R-dependent pieces that
    survive at finite u = eps m are linear in u: chiral zero points x^{(1-R_i)|nu_i|/2} (sign_i (R_i-1) Q_ia /4 per node)
    and the FI R-mixing x^{-+r/2} on T,B (+-r/4). Analytic sum in the ladder sign region (vertical B-T has nu < 0):
       L: 0,  R: 0,  T: -r/4 (chirals) + r/4 (FI) = 0,  B: +r/4 (chirals) - r/4 (FI) = 0.
    Numerically (exact prefactor, x = 0.9, deep in the ladder region): d log|prefactor| / d m_a = 0 to 1e-7 on all four
    nodes for r = 0.2, 0.3, 0.45. So mu_B = mu_T holds EXACTLY: no imaginary twisted mass survives on any node.
    It needs BOTH 3d quirks:
      - vertical R = 2-2r (fixed by the superpotential triangles): with R_vert = 1-r the residue is (0,-7/3,+7/3,0) eps r/4;
      - FI R-mixing x^{+-r/2} on T,B (fixed by the monopole superpotential, R[V] = 2): without it the residue is
        (0,-1,+1,0) eps r/4, i.e. exactly along the monopole's charge vector (T-, B+).
    Reading: the monopole superpotential's 2d imprint is precisely the absence of a twisted mass along w_B/w_T; the
    vertical chiral's monomial is O(eps) (triangle rule) but its ZERO-POINT is what the monopole R-mixing cancels against.
    (2) Monopole as a twisted-chiral exponential: in the index it enters only through the fugacity constraint; no
    O(1) monomial w_B/w_T appears in the ladder region (vertical monomial O(eps); a surviving one would give 10 vacua).
    Whether a 2d e^{-Y}-type term from the monopole sits at higher order is NOT derived.
    (3) q and monomial coefficients from the map: next.
34. Coefficient step (3) (bessel_chir.py). Single chiral in the Bessel regime x^{nu} ~ eps, checked against the exact
    q-Pochhammer ratio (6 angles):
       (x^{2-R+nu}/zeta; q)/(zeta x^{R+nu}; q) = exp[(zeta x^{R+nu} - x^{2-R+nu}/zeta)/(1-q)] (1 + O(eps))
    max rel. err 8.8e-2, 4.3e-2, 2.2e-2, 1.1e-2 at eps = 0.2, 0.1, 0.05, 0.025 (the (-1)^nu sign lives inside zeta; a
    first version that also multiplied by s double-counted it, fixed).
    With w = x^{nu} zeta: X+ = x^R w, X- = x^{2-R} wbar, so each chiral gives (x^R w^Q - x^{2-R} wbar^Q)/(1-q): the
    W - Wbar of the S^2 LG. CORRECTION to LOG 29: X- is NOT suppressed; it is the antiholomorphic partner. The
    triangle rule (a chiral whose |nu| is the sum of the other two is O(eps)) is unaffected: it is a statement about x^{|nu|}.
    Coefficients: c_i = f_i x^{R_i}/(1-q) -> f_i/eps. Rescaling invariants of the 6-term ladder:
      rho = c_LB c_TL/(c_RB c_TR) = x^{2r}/x^{2r} = 1 EXACTLY (the diagonal chirals carry no flavour fugacity; the masses
            live in the column FIs, i.e. in the mu.log w terms) -> the EHX/BCFKvS normalisation, as used in LOG 32;
      q   = c_B c_LB c_TL c_{-T} = eta^{-1} x^{n+4r}/(1-q)^4 ->  ln q = -eps n/2 - 4 ln eps - i theta_eta + O(eps).
    Electric side (LOG 9-12): kappa_el = -eps n - 8 ln(1-q). So ln q_LG = kappa_el/2: the one-loop running coefficients
    4 (mirror) and 8 (electric) agree up to the same factor 2 as eps n vs eps n/2, i.e. the sigma normalisation.
    Monopole exponential at subleading order: no trace found; not further investigated.
35. Potential from exact sectors (sector_pt.py). Exact 3d mirror sectors deep in the ladder region, fluxes
    m = (-2,-1,-3,-2) nu* + delta (delta random in [-3,3]^4), eta-flux 4 nu*, x^{nu*} ~ eps, random gauge angles (8 points):
    exact product of the 7 chiral q-Pochhammer ratios vs exp[sum over the 6 ladder chirals (X+ - X-)/(1-q)]:
       max |exact/LG - 1| = 0.71, 0.35, 0.20 at eps = 0.2, 0.1, 0.05   (O(eps))
       vertical B-T chiral alone: max |value - 1| = 0.48, 0.22, 0.10  (O(eps): triangle rule confirmed on exact numbers)
    So pointwise, the exact 3d integrand in this region is e^{W - Wbar} of the ladder at O(eps), the vertical chiral
    contributing only its (exact, R-cancelling) zero point.
36. Elliptic genus (eg.py, eg2.py; nome 0.15, R=0, equivariant masses M_j, mu = (M3-M2, M2-M1, M4-M3, M3-M2)).
    (A) Gr(2,4) GLSM by JK residues: Witten index (y -> 1) = 6 exactly; EG computed at two y's.
    (B) The ladder's "inverse Hori-Vafa" GLSM: one chiral per ladder monomial, gauge U(1)^2 = the two exponent relations,
        charges (1,1),(1,0),(0,1),(1,0),(0,1),(1,1) on (X, Y/X, Z/X, W/Y, W/Z, q/W), masses from E^T m = mu.
        Witten index = 8 (JK chamber eta = (1, 0.37); the mirror chamber gives 8 too), EG ratio to Gr(2,4):
        1.33, 1.135-0.030i, 1.041-0.035i at y = 1, e^{2 pi i 0.13}, e^{2 pi i (0.21+0.04i)}. No 6-subset of its 8 JK
        residues reproduces the Gr(2,4) EG (best mismatch 2.1 on |EG| ~ 25).
    Interpretation: this smooth toric GLSM has two FIs, which in its HV mirror are the two ladder invariants q and rho.
    Generic rho -> 8 vacua (matches LOG 30's count for generic coefficients); the 3d reduction fixes rho = 1 EXACTLY
    (LOG 34), where 2 vacua run to infinity: rho = 1 is a WALL of the toric GLSM (non-compact Coulomb direction; the
    toric degeneration of Gr(2,4) is singular there). So the EG of the reduced mirror is NOT the EG of a smooth toric GLSM;
    it is the EG of the ladder LG on the wall, i.e. a non-compact LG (Troost 2010; Ashok-Doroud-Troost 2013), which
    is generically non-holomorphic. The vacuum match (LOG 32) is consistent with this; the EG test is NOT passed by this
    route and needs the non-compact computation.
37. Option 3, wall-crossing (eg_wall.py). The ladder's toric GLSM has rays (1,0),(0,1),(1,1); its two geometric phases
    cone((1,0),(1,1)) and cone((1,1),(0,1)) are separated by the ray (1,1), i.e. t1 = t2 <-> rho = 1, the reduction's value.
    - EG identical in both chambers (to all printed digits, nome 0.15 and 1e-7): no jump across the wall (flop invariance).
    - chi_y (nome -> 0), rigid (mass-independent; checked on 3 mass points, one sample at y=e^{2 pi i 0.13} off by a
      branch issue):  chi_y[Gr(2,4)] = y^{-2}(1 + y + 2y^2 + y^3 + y^4)   (Hodge 1,1,2,1,1)  -- matches the known answer;
                      chi_y[toric]   - chi_y[Gr(2,4)] = y^{-1} + y = y^{-2}(y + y^3)     (1.369095, -0.736250 at two y,
                      exact),
      i.e. the toric phase has Hodge numbers 1,2,2,2,1: exactly one extra class in H^{1,1} and its dual in H^{3,3}.
    Reading: the toric phases are the two small resolutions (related by a flop) of the singular toric degeneration X_0 of
    Gr(2,4); the extra (1,1) class is the exceptional curve, whose volume ~ log rho vanishes on the wall rho = 1. The two
    "extra" vacua of LOG 30/36 are these two classes; at rho = 1 the curve shrinks and they run off to infinity (non-compact
    Coulomb direction). The 3d reduction lands exactly on the conifold wall, where the vacuum count is chi(Gr(2,4)) = 6.
    Consistent with the duality; not yet a proof that the wall LG has the Gr(2,4) EG (option 1 = the non-compact computation).
38. Option 1, first attempt (eg_onwall.py): naive on-wall JK. Splitting the 12 non-degenerate pole pairs of the toric
    GLSM by their charge cone relative to the wall ray (1,1): S_int {(1,0),(0,1)} (4 poles), S_A {(1,0),(1,1)} (4),
    S_B {(0,1),(1,1)} (4); S_A = S_B identically. No combination (S_int, S_int + S_A = either chamber,
    S_int + (S_A+S_B)/2, ...) reproduces the Gr(2,4) EG at nome 0.15 or 1e-7 (e.g. y=e^{2 pi i 0.13}: Gr 25.12-3.67i,
    S_int 51.6+49.8i, chamber 28.41-4.92i). As expected: the holomorphic part of a non-compact EG is not a JK subset;
    it needs the Troost / Ashok-Doroud-Troost treatment of the continuum (Appell-Lerch-type completion).
    STATUS: option 1 is a genuine research computation, not done. Stopped to ask Anant.
    Cheaper intermediate test available: D(tau, y, M) = EG[small resolution] - EG[Gr(2,4)] has chi_y limit y^{-1}+y
    exactly (LOG 37); check whether D at finite nome equals the local elliptic genus of the exceptional geometry
    (P^1-fibration over the conifold curve C = P^1, normal bundle O(-1)+O(-1)), equivariantly.
39. Option 1b (eg_local.py), fixed-point decomposition of the full equivariant EG (nome 0.15, two y, generic masses):
    the 8 JK/fixed-point terms of the toric small resolution vs the 6 of Gr(2,4):
      - FOUR toric fixed points coincide EXACTLY (to 1e-9) with four Gr(2,4) vacua:
          (Y/X, Z/X) = {M2,M4}, (Y/X, W/Z) = {M2,M3}, (Z/X, W/Y) = {M1,M4}, (W/Y, W/Z) = {M1,M3}
        (Gr vacua labelled by the pair of flavours, 0-based (1,3),(1,2),(0,3),(0,2)).
      - The other FOUR toric fixed points, (X, Y/X), (X, W/Y), (Y/X, q/W), (W/Y, q/W), replace exactly TWO Gr vacua,
        {M1,M2} and {M3,M4}.
      - Hence D = EG_res - EG_Gr = [4 exceptional-locus terms] - [2 Gr terms], at finite nome, equivariantly
        (3.291843-1.254030i and -2.462352+1.218492i at the two y; chi_y limit y^{-1}+y, LOG 37).
    Geometric reading: the toric degeneration X_0 of Gr(2,4) is singular along a curve C = P^1 with two torus-fixed points;
    those are exactly the Gr vacua {12} and {34}; the small resolution replaces each by the two ends of an exceptional P^1
    (2 -> 4 fixed points, +2 to the index, Hodge 1,2,2,2,1). Away from C the two spaces agree fixed point by fixed point.
    So at the level of the full equivariant EG, the reduced mirror (ladder on the wall rho = 1) differs from Gr(2,4) only
    at the two vacua {12}, {34}, where the exceptional curve shrinks. Option 1b PASSES as a consistency check.
    What remains for option 1a: show that the non-compact wall theory turns the 4 exceptional terms into the 2 Gr terms
    (+ a non-holomorphic remainder). This is now a LOCAL problem: two copies of a local conifold-type model.
40. Option 1a, setup (no numerics yet). Locally near each of the two Gr vacua {12}, {34} the toric degeneration X_0 is
    (conifold) x C, its small resolution (resolved conifold) x C, and Gr(2,4) the smoothing (deformed conifold) x C.
    So the local test is a conifold-transition statement for elliptic genera: EG(resolved) - EG(Gr) at that vacuum
    (2 exceptional terms - 1 Gr term) must equal what the non-compact wall theory removes. Natural model: the
    near-conifold SCFT is the SL(2,R)/U(1) cigar (N=2 Liouville) tensored with spectators (Ghoshal-Vafa), whose EG is
    holomorphic part = Jacobi form x Appell-Lerch sum A_{1,k}, plus a non-holomorphic completion; literature reading:
    holomorphic (mock) part = states localised on the blown-up 2-cycle, non-holomorphic part = bulk continuum.
    Troost 1004.3649 not fetchable here (allowlist); ar5iv / 1311.0918 (ADT GLSM) / 1706.02576 / 1905.05207 reachable.
    Plan: (i) local equivariant EGs of resolved and deformed conifold x C from the fixed-point data already in eg_local.py;
    (ii) cigar EG at k = 1 with the Appell-Lerch sum (ADT 1311.0918 App. A conventions); (iii) test
    [2 exceptional terms] - [Gr term] = - (holomorphic cigar piece x spectator), vacuum by vacuum, at finite nome.
41. Option 1a, step (i) DONE (eg_weights.py, eg_conifold.py). Tangent weights at the fixed points show the local conifold
    exactly: at Gr vacuum {12} the weights {a,b,c,d} = {M3-M1, M4-M1, M3-M2, M4-M2} obey a+d = b+c =: f (the weight
    relation of xy = zw, i.e. Gr is locally the smoothing xy - zw = s with s of weight f along the curve C); the two
    replacing toric fixed points are (Y/X, q/W) with weights {a, c, e, f} and (W/Y, q/W) with {b, d, -e, f}, e = b - a
    (exceptional P^1 weight +-e, curve direction f). Same at {34} with (X, Y/X), (X, W/Y).
    Universal local function, phi(w) = theta1(w - z)/theta1(w):
       D_loc(a,b,c) = phi(f)[phi(a)phi(c)phi(e) + phi(b)phi(d)phi(-e)] - phi(a)phi(b)phi(c)phi(d),  d = b+c-a, f = b+c, e = b-a
    and the full equivariant difference EG_res - EG_Gr = D_loc({12}) + D_loc({34}) EXACTLY (8 digits, nome 0.15, two y).
    So the whole EG question is: does the non-compact (wall, FI = 0) local conifold theory, U(1) with charges (1,1,-1,-1)
    x the spectator C_s, have holomorphic EG = phi(a)phi(b)phi(c)phi(d) (the smooth Gr point), i.e. does its continuum
    remove exactly D_loc? Steps (ii)-(iii) (cigar / Appell-Lerch at k = 1) not done yet.
42. Option 1a, steps (ii)-(iii): VERDICT = the target fails as formulated (analytic, plus LOG 38/39 numerics).
    - Local model: U(1) with charges (1,1,-1,-1) and flavour weights a, d (charge +1) and b, c (charge -1) with the
      conifold relation a + d = b + c. Its EG integrand is elliptic in the gauge variable (y^{sum Q} = 1), so the sum of all
      residues vanishes: JK(eta > 0) = JK(eta < 0) = the 2-term resolved answer (seen numerically: S_A = S_B, LOG 38).
      A plain chiral GLSM has no tau-bar dependence, so there is no Appell-Lerch completion to subtract: the
      holomorphic "wall" EG of this model is just the resolved one. The continuum cannot remove D_loc. (The ADT
      non-holomorphy needs a non-compact compensator / Stueckelberg field, which this model does not have.)
    - Vacua: at FI = 0 the 2d Coulomb equation (s+a)(s+d) = (s+b)(s+c) becomes s(a+d-b-c) + (ad - bc) = 0; with the
      conifold relation a+d = b+c it has NO finite solution (ad != bc generically): both local vacua run to infinity even
      with twisted masses. This is the mechanism behind 8 -> 6 (LOG 30/36/37): the reduction lands at FI = 0 AND the
      equivariant weights automatically satisfy the conifold relation.
    So: vacuum count of the 2d endpoint = 6 = chi(Gr(2,4)), but its (naive, holomorphic) EG = EG(small resolution)
    (Witten index 8), not EG(Gr(2,4)); the 2 extra states are continuum states at infinity. The EG of the 2d limit theory
    is not a continuous function of the parameters at this point. The 3d index identity is exact, so this is an
    order-of-limits / non-compactness issue of the 2d endpoint, not a failure of the 3d duality.
    Stopped to ask Anant how to read this.
43. Mirror 3d Bethe-vacuum count, EXACT (Singular std + vdim, saturated by z_a and 1 - X_i, random rational fugacities,
    two seeds agree; bme_singular.py, bme_scan.py, bme_scan*.log). Per-chiral quantisation sign sigma_i (chiral factor
    X^{sigma_i/2}/(1-X), K' = K + 1/2 sum sigma_i Q_i Q_i^T), all 2^7 = 128 assignments:
      counts that occur: 4,5,6,7,8,9,10,12,14,16,20,24,28. Uniform sigma = +1: 12; uniform -1: 28 (so the earlier random-
      Newton 12 and 28 were right; I had wrongly discarded them).
      Electric (exact, LOG 31): Higgs chamber 6, other chamber 20.
      - 20 occurs on the mirror (2 L<->R-symmetric assignments + 2 others).
      - 6 occurs ONLY for 8 assignments that break L<->R: vertical B-T = +1, (B, -T) = (+,-) or (-,+), and exactly one
        of the four diagonals flipped relative to the others' pattern (list in bme_scan2.log).
    So the electric Higgs-chamber index 6 IS reproduced by the mirror in a chamber where the real masses distinguish the
    L and R columns (natural: an ordering of the SU(4) real masses distinguishes the columns). The exact map
    electric real-mass chamber -> mirror sign assignment is NOT derived; only existence of a matching mirror chamber.
    Consequence for the EG puzzle (LOG 42): the 3d genus-1 index in the Higgs chamber is 6, on both sides. The 2 runaway
    states of the 2d endpoint are not counted by the 3d theory; the naive endpoint EG (index 8 = the resolution) over-counts.
    Supports reading (a): endpoint EG = EG(resolution) - D_loc{12} - D_loc{34} = EG(Gr(2,4)), with the subtraction
    meaning "drop the states at infinity of the conifold Coulomb branch". Index-level only; not a derivation of the
    y-, tau-dependent EG.
44. Item 1 (chamber map) RESOLVED, and LOG 31/43 CORRECTED (bme_constr.py, Singular, seeds 3 and 8).
    - The chiral factor X^{+1/2}/(1-X) with the index's bare K is the theory defined by the index (symmetric quantisation,
      zero point zeta^{-|nu|/2}); X^{-1/2}/(1-X) differs by a factor X^{-1} per chiral, i.e. CS levels shifted by -Q Q^T:
      a DIFFERENT theory, not another real-mass chamber. So electric "20" (LOG 31) and mirror "28/22" are other theories.
    - My mirror counts in LOG 43 used generic chiral fugacities, which VIOLATE the superpotential: the planar triangles
      force f_LB f_TL f_BT = 1 and f_RB f_TR f_BT = 1. With these imposed:
         mirror, symmetric quantisation: 6 vacua (both seeds; with t_T = +t_B, -t_B, or t_T free)  = electric 6.
         controls: constraints off -> 12 (the LOG 43 number); constraints on, sigma = -1 -> 22.
    So the 3d Witten index matches, 6 = 6, in the theory itself; no special/L-R-breaking chamber is needed. The LOG 43
    reading ("6 only in L<->R-breaking chambers") was an artefact of ignoring the superpotential. The monopole constraint
    (t_T vs t_B) does not change the count. The genus-1 index of both 3d theories is 6 = chi(Gr(2,4)); the 2 runaway
    states of the 2d endpoint are not counted in 3d (supports reading (a) of LOG 42, index level).
45. Item 2(a), chart analysis (charts2.py). Exact chiral factor in the Bessel regime, by sign of nu (derived from the
    q-Pochhammer, LOG 34 extended): nu > 0: exp[(x^R w^Q - x^{2-R} wbar^Q)/eps]; nu < 0: exp[(x^R / wbar^Q - x^{2-R}/w^Q)/eps].
    So the holomorphic term is w^Q on one side and -1/w^Q on the other: a piecewise (tropical / mutation-type) W, not one
    global Laurent polynomial. In the ladder region the vertical chiral has nu ~ -2 nu*, its term is -1/w^{Q_BT} ~ O(eps):
    the true W there is ladder + O(eps) (T/B) -- consistent with LOG 35's pointwise check.
    Charts (one chiral dropped per triangle): ladder (B-T dropped) has 6 critical points (masses on); the four charts
    with B-T active and one of {L-B,T-L}, one of {R-B,T-R} dropped have 2 each (5 monomials).
    With the reduction's coefficients the alternative charts' invariant is ~ x^n/eps ~ eps^3 (ladder q = O(1)); their
    2 critical points drift to |log w| ~ 8.4, 12.0, 14.8 at eps = 0.1, 0.03, 0.01 (one lost by the solver at 0.01),
    i.e. ~ 3 ln(1/eps) away from the ladder region. Whether they are self-consistent saddles of the exact 3d integrand
    (they would add vacua beyond 6, contradicting the exact 3d count of 6, LOG 44) is NOT settled: needs the raw-region
    check (signs and sizes of all nu_i at those points). Expectation: inconsistent (exit their chart). Open.
    Item 2(b) (numerical region sums) not started.
