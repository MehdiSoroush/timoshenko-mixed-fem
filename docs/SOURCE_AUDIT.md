# Historical model correspondence and numerical reconstruction

**Historical reference:** A 2017 graduate-level FEM coursework report (64 PDF pages). The archival report is not included in this repository. PDF page references below identify locations in that source.

## Archival source correspondence

| Feature | Source location (PDF page, 1-indexed) | Status |
| --- | --- | --- |
| Mechanical problem: 2 m beam, left-half distributed load, central point load, support schematic | 9-10 | Configuration adopted |
| Material data and cross-section | 10, 32, 38 | Data adopted |
| Mixed beam formulation (`w` + bending moment `M`) | 17-27 | Formulation reviewed; not independently audited in full |
| Local 4x4 static stiffness matrix | 38 | **Corrected:** printed static listing has `-a` in final diagonal, inconsistent with PDF p.33 dynamic `+a`. The corrected `+a` converges and agrees with independent FEM; see independent verification |
| Static force assembly, constraints, linear solve | 39-40 | Reimplemented; some obvious numerical issues improved |
| Dynamic mass and time integrator | 33-37 | **Implemented; independent continuum bending-frequency and separate physical FE transient check now added; archival Abaqus datasets unavailable** |
| Abaqus comparisons and convergence figures | 42-63 | Historical reported results only; original numeric data unavailable |

## Numerical implementation choices

1. Separate `beam.width` from the source's element matrix variable `b`.
   Source PDF p. 38 uses `B=0.03`, `A=b*h` and `I=b*h^3/12` while it later
   defines `b` as an element stiffness coefficient; this is inconsistent
   in literal MATLAB variable case. Here dimensions and coefficients do
   not share names.
2. Use a sparse global stiffness matrix and sparse linear solver, not
   `inv(Kt)*Ft` (PDF p. 40).
3. Explicitly preserve the **three** constrained mixed DOFs in the
   fixed-pinned source configuration: `w(0)=0`, `w(L)=0`, `M(L)=0`.
4. For a force located between nodes, use linear interpolation to the
   neighboring deflection DOFs (source assumes odd `n`, center at a node).
5. For distributed loads ending inside an element, integrate the shape
   functions on the loaded subsegment, rather than assume load transitions
   must coincide with element edges.
6. Retain the **non-symmetric source element matrix** rather than
   silently replacing the author's mixed formulation with a standard
   displacement/rotation Timoshenko element. **Important correction:** the
   original printed *static* last diagonal is `-a`, while the dynamic
   listing uses `+a`. The printed `-a` form gives nonconvergent deflections; the
   implementation therefore uses the independently validated `+a`. The separate
   standard model is used for verification only (see `INDEPENDENT_VERIFICATION.md`).

## Analytical and numerical benchmarks

- **Analytical, independent:** Standard simply supported Timoshenko beam under
  a *full-length* uniform distributed load. The exact midspan deflection is
  `5 q L^4 / (384 EI) + q L^2/(8 k_s G A)`.
- **Analytical, independent:** Standard simply supported Timoshenko beam
  under a *central concentrated force*. Exact midspan deflection is
  `Q L^3/(48 EI) + Q L/(4 k_s G A)`.
- **Source study, numeric only:** Mesh refinement of the fixed-pinned beam
  subject to uniform `q=4800 N/m` on the left half and `Q=3000 N` at midspan.

These closed-form benchmarks verify the static implementation in special cases. Additional frequency and transient comparisons with an independently formulated physical FEM model are reported in [independent verification](INDEPENDENT_VERIFICATION.md). The historical Abaqus datasets are unavailable, and their exact reproduction has not been established.

## Archived static-deflection comparison

The original plot (PDF p. 42, printed figure 4-1) appears to have a peak
near 0.019 m, while the parameters and static matrix/BCs as transcribed
from PDF pp. 38-40 give a peak near 0.00966 m at 101 nodes.
The documented computational parameters do not reproduce the archived peak. Possible explanations include unrecorded runtime settings, a different load specification, or an error in the plotted data or description.

## Verification scope

The mixed dynamic mass matrix and historical time-integration recurrence have been implemented. The static and selected dynamic responses have been compared with analytical or independent FEM references. Stress recovery has been audited with separate definitions of outer-fiber and historical quarter-thickness normal stress.

The available numerical evidence covers the specified static, first-mode, and selected transient benchmarks. It does not establish equivalence with the unavailable archival Abaqus data or encompass all vibration modes and loading histories.

## Research provenance

This is a modern **Python reconstruction** of a 2017 MATLAB coursework model,
not the recovered MATLAB source code. The numerical tests, datasets, and figures are computed using the 2026 Python implementation; no experimental measurements are included as validation data.

