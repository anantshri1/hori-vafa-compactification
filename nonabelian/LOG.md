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
