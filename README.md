# Mixed Finite Element Analysis of Timoshenko Beams

**Mehdi Soroush** · Research study: 2017 · Python reconstruction: 2026



This repository presents a reproducible mixed finite element implementation for static and transient Timoshenko beams. The mathematical formulation originated in a graduate computational mechanics study in **2017**; the executable Python reconstruction, independent numerical reference models, and verification studies were developed in **November 2026**.

The implementation combines a displacement–moment mixed formulation with a separately assembled displacement–rotation finite element model, analytical benchmarks, and reproducible static and dynamic calculations.

The [English technical report](docs/TECHNICAL_REPORT.md), also available as a [PDF](docs/Technical_Report.pdf), presents the mechanics, numerical method, validation data, figures, interpretation, and limitations. Additional documentation covers the [mathematical formulation](docs/MATHEMATICAL_FORMULATION.md), [historical model correspondence](docs/SOURCE_AUDIT.md), [dynamic formulation](docs/DYNAMIC_AUDIT.md), and [independent numerical verification](docs/INDEPENDENT_VERIFICATION.md).

## Implemented numerical methods

- **Static mixed formulation:** displacement `w` and bending moment `M` as alternating nodal DOFs; non-symmetric source-derived element matrix, consistent load assembly, sparse solve.
- **Dynamic mixed formulation:** source-derived non-symmetric 4×4 element mass matrix and assembled second-order system `M q_ddot + K q = f(t)`.
- **Weighted three-level integration:** Crank-Nicolson/trapezoidal startup and weighted-residual three-level recurrence with `gamma=1.5`, `lambda=0.8` (algorithmic damping).
- **Additional temporal reference:** Newmark average-acceleration integration for comparison with the source-derived recurrence.
- **Loading:** a point force at the beam midpoint, and three historical distributed load shapes (`constant`, `sin(t)`, `source-ramp`) active only during the first second, plus a separate unit-height triangular comparison profile.
- **Verification:** independent static closed-form benchmarks, mesh convergence, source-matrix tests, generalized FE eigenanalysis, exact semidiscrete modal time solution, and automated tests.
- **Independent verification:** a separate conventional displacement–rotation FEM with its own stiffness and physical mass matrices; Timoshenko continuum dispersion benchmark; independent dynamic response comparison; explicitly labeled stress recovery and source-code inconsistency audit.

## Default source study model

| Quantity | Value |
| --- | --- |
| Beam length | 2 m |
| Rectangular cross-section | `b=0.03 m`, `h=0.0431 m` |
| Young's modulus / Poisson ratio | `E=200 GPa`, `nu=0.3` |
| Density | 7800 kg/m³ |
| Shear correction | 5/6 |
| Distributed load | `P0=4800 N/m` over left half, for first 1 s in dynamic study |
| Concentrated force | `Q=3000 N` at midspan, throughout dynamic study |
| Mixed fixed-pinned constraints | `w(0)=w(L)=M(L)=0` |

The implementation assigns distinct identifiers to the rectangular-section width and the element coefficient denoted by the same symbol in the archived MATLAB listing.

## Installation and reproducibility

Python 3.11+ (no Abaqus, MATLAB, or proprietary solvers required):

```powershell
py -3.11 -m venv .venv
.venv\Scripts\Activate.ps1
python -m pip install -e ".[dev]"
python -m pytest -q
python examples/reproduce_original_static.py
python examples/verify_simply_supported.py
python examples/reproduce_dynamic_cases.py
python examples/verify_free_vibration.py
python examples/validate_independent_reference.py
python examples/recover_stresses.py
python examples/audit_historical_static.py
```

On Linux/macOS, the corresponding environment commands are `python3 -m venv .venv` and `source .venv/bin/activate`. Python 3.11 or newer is required; the project-local virtual environment isolates its dependencies.

## Numerical verification

| Claim | Level of evidence |
| --- | --- |
| Mixed static element stiffness | Uses the `+a` coefficient present in the dynamic listing, supported by mesh convergence and an independent FE reference |
| Static beam validation | Closed-form simply-supported uniform/point-load comparisons |
| Mixed dynamic element mass | Reproduced from original matrix specification |
| Dynamic time solver | Source recurrence reproduced and algebraically tested |
| Dynamic modal validation | Independent FE generalized eigenmode + exact semidiscrete solution, with new Newmark convergence |
| Independent continuum frequency | Exact simply-supported Timoshenko bending frequency 24.7214 Hz; mixed 41-node result 24.7341 Hz |
| Independent physical FE validation | Separate [w, rotation] elements: 0.0168% static midspan difference at 129 nodes; dynamic RMS 0.699% of reference peak over 0.2 s at 41 nodes |
| Stress-recovery audit | Explicit outer-fiber, quarter-thickness, section-average, maximum rectangular and effective Timoshenko shear definitions |
| Historical Abaqus comparisons | **Not independently reproduced**; underlying Abaqus datasets absent |
| Historical static deflection figure | Archived plot differs from calculations using the documented model inputs |
| Continuum-level dynamic verification | First-frequency and selected short-time responses verified; additional modes and loading histories are outside the reported cases |

**Stiffness-matrix convention.** The archived static MATLAB listing contains `−a` in the final stiffness diagonal, whereas the dynamic listing contains `+a`. The former is nonconvergent under the documented refinement study. The Python implementation uses the independently verified `+a` coefficient.

**Time-integration response.** The weighted three-level recurrence exhibits algorithmic damping in the reported numerical tests. The separate Newmark implementation provides a temporal comparison; see the [dynamic verification](docs/DYNAMIC_AUDIT.md).

## Main results

Run the scripts to regenerate the tracked CSV summaries and figures under `results/`:

- `results/figures/original_case_mesh_convergence.png`
- `results/figures/analytical_validation_convergence.png`
- `results/figures/dynamic_three_loading_cases.png`
- `results/figures/time_integrator_comparison.png`
- `results/figures/modal_time_integration_verification.png`
- `results/figures/first_frequency_mesh_convergence.png`
- `results/figures/independent_fem_validation.png`
- `results/figures/independent_frequency_verification.png`
- `results/figures/verified_stress_recovery.png`
- `results/figures/static_matrix_sign_audit.png`

## Repository layout

```text
src/timoshenko_mixed/       # Mixed static/dynamic FE kernels and benchmarks
examples/                  # Reproducible example scripts
results/data/              # Generated CSV datasets
results/figures/           # Generated plots
tests/                     # Numerical and regression tests
docs/                      # English report, mathematical notes, numerical audits
```

## Linear-elastic interpretation

Stress recovery is based on small-deformation, linear-elastic Timoshenko beam theory. The archived model specifies elastic properties without identifying a particular steel grade or yield strength. Plasticity, material failure, and structural strength assessment are outside the scope of the formulation. Evaluation conventions are presented in [independent verification](docs/INDEPENDENT_VERIFICATION.md) and the [mathematical formulation](docs/MATHEMATICAL_FORMULATION.md).

## Research provenance and license

The mixed finite element study originated in **2017**; its documented Python reconstruction and independent verification work date to **2026**. The archived static-deflection plot differs from the response computed using its stated parameters, and the original Abaqus numerical datasets are unavailable. The [historical model correspondence](docs/SOURCE_AUDIT.md) and [independent verification](docs/INDEPENDENT_VERIFICATION.md) describe these comparisons.

The software is distributed under the [MIT License](LICENSE); see [CITATION.cff](CITATION.cff). The 2026 computational results are distinct from archived Abaqus output or experimental measurements.
