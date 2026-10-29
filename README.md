# Reproducibility package: 3d → 2d reductions (26–27 Oct 2026)

Requirements: Python 3, `numpy`, `scipy`, `mpmath`. Run each script from this folder: `python3 <script>`.
Runtimes range from seconds to ~20 min (noted). Every claim below is "computed" (a script prints it) unless marked otherwise.

## Claim → script → what you should see

| # | Claim | Script | Expected output |
|---|---|---|---|
| 1 | Single q-Pochhammer / chiral asymptotics; flux sector ν localises at ζ = (−1)^ν | `step0_chiral.py` | all `|ratio-1|` halve as eps halves |
| 2 | U(1)_{1/2}+Φ ↔ free chiral holds only with an integer CS sign, not (−1)^{m/2} | `step1a_identity.py` | `+1, Chir(w)` at ~1e-14; e^{±iπm/2} rows O(1) |
| 3 | FI-scaled limit: fixed flux carries <1% (CS side) | `step1b_smallbeta.py` (~10 min) | `share |m|<=3` → 0.008 at eps 0.05 |
| 4 | Bessel–Mellin mechanism, O(eps²) | `step1c_mellin.py` (~15 min) | `|I_E/pred-1|` 7e-4, 1.8e-4, 4.5e-5, 1.1e-5 |
| 5 | Background-flux contact terms for U(1)_{1/2} pair | `step1d_flux_identity.py` | first scan fails (as reported); see log for the fitted prefactor |
| 6 | T[SU(2)] exact self-mirror I(Y,W,a)=I(W,Y,1/a) | `tsu2_lib.py` + snippet in log | 1e-15 |
| 7 | T[SU(2)] Higgs copy → T*CP¹ GLSM; Coulomb copy → continuum | `tsu2_limit2.py` (~15 min) | sector errors O(eps); copy-B profile collapse to 4 digits |
| 8 | T[SU(2)] 2d identity GLSM = Coulomb gas | `tsu2_coulomb.py` (~20 min) | ~3e-4 (limited by slow GLSM flux sum) |
| 9 | CS pair exact 3d map (symmetric gauge) | `pair_map.py` (~20 s) | best candidate err ~1e-15: wB=y^2, a=wA^{1/2}, rB=(1-rA)/2 |
| 10 | Same duality in Anant's gauge (both R = 1-r, η on Q̃, A-side η = w x^{1-r}) | `pair_hisgauge.py` + snippet in log | ~4e-15 |
| 11 | Operator check: triplet currents at x², dressed-monopole triplet at x^{3+r} | `pair_ops.py` | coefficients −2, +1 as reported |
| 12 | Flux map (w x^{3r-1})^{nA/4}, nax = nA/2 | `pair_fluxmap.py` | `nA/2` rows have constant |R| |
| 13 | n>0 is the CP¹ Higgs phase; sectors → GLSM at O(eps) (after measure fix) | `pair_sign.py` (~10 min) | ratio/eps → 1 |
| 14 | CS pair 2d identity converges O(eps) | `pair_2d.py` (~15 min) | 6.1e-2, 3.0e-2, 1.5e-2 |
| 15 | **Pure 2d**: Z_S²[CP¹] = e^{rκ/2} ∫dρ e^{−2iμρ} J0(4e^{κ/4}√(sinh²(ρ/2)+sin²(θ/2))) | `cp1_exact.py` (~1 min) | rel diff ~1e-16 at 4 points |
| 16 | T*CP¹ elliptic genus (JK): η-independent, →2, = −HLM(3.43) | `eg_tcp1.py`, `eg_hlm.py` | 1e-21; 2.000015; ratio −1.0 |
| 17 | HLM Taub-NUT genus → ALE only as g²→∞ | `eg_hlm.py` | 0.16 … 0.0073 |
| 18 | CP¹ genus defined only at y = ±1; vanishes at y = −1 | `eg_cp1.py` | boundary integral O(1) except z=1/2; E ≈ 1e-16 |
| 19 | Hankel split of J0 is NOT the vacuum basis (Re parts agree, Im parts don't) | `cs_vacua.py` (~2 min) | e.g. Im Z_+ = −2.184 vs Im F1/2 = +0.402 |
| 20 | Vacuum by vacuum: Z_s = e^{rκ/2} π/sin(πa_s) I_{−a_s}(x)I_{−a_s}(x̄); Z = e^{rκ/2}[2(KĪ+IK̄) + (4 sin πa/π)KK̄]; I_{−a} = I_a + (2 sin πa/π)K_a | `cs_blocks.py` (~3 min) | 1e-17 (A), 1e-20 (B); \|Z_+\| ~ 1/μ as μ→0, sum finite |
| 21 | Direct S² contour integral = residue sum (fixes 2π and sign of the GLSM form) | `cp1_contour_full.py` (~10 min) | 1e-16 at all 4 points (old `cp1_contour.py` used the wrong oscillation period) |
| 22 | Blind subagent rederivation of claims 9–15 | `blind_check/` + `blind_subagent_prompt.md` | same map, same J0 formula, same values at 4 points |

## NOT computed (derivations/arguments only — treat as unchecked)
- Copy B ≅ two NS5 on a circle (KK-summed one-loop metric), T-dual to HLM's Taub-NUT family (Tong hep-th/0204186).
- Radius dictionary g² = R e₃d²/(2πk), k an unfixed O(1) normalisation. (My first version, g² = π e₃d² R with sinh(2πRr), was WRONG: caught by comparing with Tong (1.4) and a numerical Poisson check.)
- Z_{S²} independent of the Taub-NUT radius.
- Novelty of claim 15 (literature check pending).
- That the off-diagonal 2 in claim 20 is the Cecotti–Vafa CP¹ soliton number (normalisation-dependent; unchecked).

## Known errors of mine caught during the work (for the method write-up)
Wrong (−1)^{km} reading; dilogarithm regime misplaced (1b); "second CFS limit" mislabel; double-counted measure;
wrong axial phase sign; misread of [1,1]; "A's mesons" claim; n<0 phase misidentified as Higgs; dropped dφ = ε dσ;
J0 scale e^{−κ/4} (should be e^{+κ/4}); 4e-4 residual was quadrature, not physics;
first vacuum test used the Hankel split (wrong basis — an angular-contour split, not thimbles);
contour check treated a chirp (phase ~ −4σ log σ) as a fixed-frequency oscillation (1e-5 error).
