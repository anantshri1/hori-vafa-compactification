# Blind verification task: a 3d N=2 Chern–Simons mirror pair and its 2d limit
(draft by Claude, 27 Oct, for Anant to review before launch; the subagent gets ONLY this text)

You are checking a calculation independently. You have Python with numpy, scipy and mpmath. You have no access to
any earlier work on this problem; derive everything yourself. Report failures and dead ends honestly. A clear
"I could not make it work, here is where it breaks" is worth more than a forced match.
Write one script per phase and put the numerical evidence in your report.

## Conventions (use these exactly)
Superconformal index on S²×S¹: x = e^{−β}, q = x².
One chiral of R-charge R, total fugacity ζ, in flux sector ν ∈ ℤ:
  Chir[ζ, ν, R] = (ζ x^{R−1})^{−|ν|/2} · (ζ^{−1}(−1)^ν x^{2−R+|ν|}; x²)_∞ / (ζ(−1)^ν x^{R+|ν|}; x²)_∞
U(1) gauge theory, gauge fugacity z, gauge flux m, contour |z| = 1, measure dz/(2πi z), sum over m ∈ ℤ.
CS level k and topological symmetry (fugacity ω, background flux n) contribute
  Z_cl = (sign) · ω^{−m} z^{−km−n}
Warning: for the sign, the literal (−1)^{km} can fail. Fix the sign by demanding a consistent result. Combine all
half-integer powers (chiral zero-point factors and CS) before evaluating.

## The pair
Theory A: U(1) at CS level +1, two chirals Q1, Q2 of gauge charge +1.
Theory B: U(1) at CS level −1, one chiral Q of charge +1, one chiral Q̃ of charge −1.
The S³ (real-mass) map, as known:
- A: real masses X1, X2 for Q1, Q2; FI parameter η; R-charge r for both chirals.
- B: FI parameter X2 − X1; one of the two chirals carries η, the other does not; R-charge 1 − r for both chirals.
- The chiral-ring statement: the gauge-invariant monopole of A maps to the meson of B.

## Phase 1: the exact index identity
(a) Translate the S³ map into an index map: fugacities and R-charges on both sides (you may change gauge, for
    example shift the A-side flavour fugacities by a gauge transformation). Find the map for which I_A = I_B holds
    as an exact identity. Verify it numerically at at least 3 generic points (r, flavour fugacities, x), to
    better than 1e-12. If the naive translation of the S³ map fails, say so, show the failure, and explain what you
    had to change and why.
(b) Switch on a background flux n for A's topological symmetry. Find the B-side background fluxes and the
    prefactor that make the identity hold sector by sector in n. Verify numerically for several n.
(c) Check the chiral-ring statement in the index: identify the A monopole and the B meson as terms of the
    x-expansion at matching charges.

## Phase 2: the 2d limit of theory B
Take x = e^{−ε/2} → 1 (ε → 0) with:
- A's flavour (axial) fugacity y = q^{−iμ}, μ fixed;
- A's topological fugacity ω = e^{iθ}, θ fixed;
- A's background topological flux n ∝ 1/ε (n > 0), chosen so that a 2d FI parameter κ stays finite. You define
  κ. State your definition.
- r fixed.
Using your Phase 1 map, work out the corresponding limit of I_B. B is the side WITHOUT a good gauge-theory (GLSM)
limit. Expect the reduction to run over a continuum of gauge fluxes rather than a fixed flux. Find:
(a) the leading asymptotics of I_B as ε → 0, as an explicit finite-dimensional integral of special functions
    (state the overall power of ε and any ε-dependent prefactors);
(b) numerical evidence that I_B approaches your formula, with the error falling as ε decreases, at at least 2
    parameter points;
(c) the value of your limiting 2d expression (with prefactors stripped, and stated) at
    (r, μ, κ, θ) = (0.3, 0.25, −2.0, 0.8), (0.4, 0.6, −1.0, 2.0), (0.2, 0.1, −3.0, 1.3), (0.45, 0.8, −0.5, 2.7),
    to at least 10 digits;
(d) if you can, identify the 2d theory your expression is the S² partition function of, and the A-side 2d theory.

## What to report
For each phase: the final formulas, the scripts, the numbers, anything you had to guess or fix by hand, and every
place where you were not sure.
