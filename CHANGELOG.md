# Changes

## Documentation maintenance — November 2026

- Updated public-facing documentation and removed obsolete release-preparation material.
- Clarified verification scope and project history without modifying numerical code or results.

## 0.4.0 — November 2026 (English documentation and publication audit)

- Added a complete English technical report (Markdown and PDF), with mathematical methods, verification tables and figures, and provenance.
- Standardized public-facing descriptions of the 2017 coursework and November 2026 Python reconstruction.
- Added CITATION.cff, MIT license, and reproducibility documentation.
- Re-ran all 30 numerical tests; no numerical solver changes in this release.



## 0.3.0 — November 2026 (Phase 3)

- Added a physically independent, selectively integrated displacement–rotation
  Timoshenko finite element static and dynamic reference solver.
- Added an exact simply supported Timoshenko first-frequency benchmark with
  rotary inertia.
- Added static and dynamic independent comparisons and reproducible CSV/plots.
- Added stress recovery with clearly differentiated outer-fiber, historical,
  section-average, maximum rectangular, and effective Timoshenko definitions.
- Audited the contradictory printed static `-a` and dynamic `+a` coefficients.
- Expanded source transparency, mathematical formulation, and numerical verification documentation; archival 2017 coursework PDF remains excluded.
- Expanded automated numerical tests from 22 to 30.

## 0.2.0 — November 2026 (Phase 2)

- Source-derived mixed mass matrix and historical three-level time integrator.
- Newmark comparison; dynamic load profiles; modal numerical checks.

## 0.1.0 — November 2026 (Phase 1)

- Reconstructed static mixed FEM and analytical beam benchmarks.
