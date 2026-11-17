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
