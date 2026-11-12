# Source audit and reconstruction notes (Phases 1–3)

**Primary source:** User-supplied `FEM_Project.pdf` (2017 coursework, 64 PDF pages).
The archival 2017 course report is **intentionally excluded** from this code-focused
repository at the author's request. Its source year is recorded without
publishing the original PDF.

## What we have reconstructed directly from the report

| Feature | Source location (PDF page, 1-indexed) | Status |
| --- | --- | --- |
| Mechanical problem: 2 m beam, left-half distributed load, central point load, support schematic | 9-10 | Configuration adopted |
| Material data and cross-section | 10, 32, 38 | Data adopted |
| Mixed beam formulation (`w` + bending moment `M`) | 17-27 | Formulation reviewed; not independently audited in full |
| Local 4x4 static stiffness matrix | 38 | **Corrected:** printed static listing has `-a` in final diagonal, inconsistent with PDF p.33 dynamic `+a`. The corrected `+a` converges and agrees with independent FEM; see Phase 3 audit |
| Static force assembly, constraints, linear solve | 39-40 | Reimplemented; some obvious numerical issues improved |
| Dynamic mass and time integrator | 33-37 | **Implemented; independent continuum bending-frequency and separate physical FE transient check now added; full historical Abaqus validation pending** |
| Abaqus comparisons and convergence figures | 42-63 | Historical reported results only; original numeric data unavailable |

## Explicit engineering decisions made in this reconstruction

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
   listing uses `+a`. Using `-a` gives nonconvergent deflections, so we use
   the independently validated `+a`. The separate standard model is used
   for verification only (see `PHASE3_VERIFICATION.md`).

## Verification actually performed

- **Analytical, independent:** Standard simply supported Timoshenko beam under
  a *full-length* uniform distributed load. The exact midspan deflection is
  `5 q L^4 / (384 EI) + q L^2/(8 k_s G A)`.
- **Analytical, independent:** Standard simply supported Timoshenko beam
  under a *central concentrated force*. Exact midspan deflection is
  `Q L^3/(48 EI) + Q L/(4 k_s G A)`.
- **Source study, numeric only:** Mesh refinement of the fixed-pinned beam
  subject to uniform `q=4800 N/m` on the left half and `Q=3000 N` at midspan.

These validation cases verify the solver in special circumstances, but
**do not independently validate** the source's original mixed dynamic model
or its own fixed-pinned dynamic comparisons.

## Unresolved source-to-figure discrepancy

The original plot (PDF p. 42, printed figure 4-1) appears to have a peak
near 0.019 m, while the parameters and static matrix/BCs as transcribed
from PDF pp. 38-40 give a peak near 0.00966 m at 101 nodes.
**Do not say that the source plot has been reproduced** until the difference
is explained. Candidate causes include a parameter used at run time that
was different from the documented values, plotting data not preserved,
force application, or a transcription/description error.

## Next phases before a public release

1. [Implemented Phase 2] Recreate the dynamic `Me` matrix and original Crank-Nicolson / three-level integration; see [dynamic audit](DYNAMIC_AUDIT.md).
2. [Phase 3 verified in selected cases] Independent modal time checks, exact continuum first bending frequency, and independent physical FE transient comparisons now exist; original Abaqus numerical comparison remains outstanding.
3. [Phase 3 complete] Check the historical `h/(4I)` vs outer-fiber `h/(2I)` stress location and generate stress plots, retaining both definitions.
4. Reconcile original static figure with reproduced numerical results.
5. Add appropriate Abaqus reference data **only if the underlying files
   can be supplied or the model is independently recreated**.
6. Add CI and an explicitly chosen open-source license prior to publishing.

## Provenance

This is a modern **Python reconstruction** of a 2017 MATLAB coursework model,
not the recovered MATLAB source code. Every new test and plot is new work.
No results in this project should be misrepresented as experimental data.

## Phase 3 update

A separate displacement–rotation solver, exact continuum frequency benchmark,
and source-vs-physical stress audit now exist. See [Phase 3 verification](PHASE3_VERIFICATION.md).
The uncorrected static source matrix and the historic figure cannot be called reproduced.
