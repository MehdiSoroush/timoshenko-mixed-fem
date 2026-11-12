# Mixed Finite Element Analysis of Timoshenko Beams

> **Status: Published Python reconstruction (v0.4.0).** Historical Abaqus data and a discrepant 2017 static plot remain unavailable/unresolved; this repository does **not** claim to reproduce those results exactly.

A reproducible Python finite element implementation built on a **2017 graduate mixed-FEM study**, with modern numerical reconstruction and independent verification in **November 2026**.

The original project analyzed a Timoshenko beam under static and dynamic loading. This repository focuses on *implementations and reproducible calculations*, not a scanned historical report.

**Start with the new English [technical report](docs/TECHNICAL_REPORT.md)** (also available as [PDF](docs/Technical_Report.pdf)). It presents the mechanics, numerical method, validation data, figures, interpretation, and limitations in one place. For implementation details, see [formulation](docs/MATHEMATICAL_FORMULATION.md), [source audit](docs/SOURCE_AUDIT.md), [dynamic audit](docs/DYNAMIC_AUDIT.md), and [independent verification](docs/PHASE3_VERIFICATION.md). The archival coursework report is not included.

## Implemented numerical methods

- **Static mixed formulation:** displacement `w` and bending moment `M` as alternating nodal DOFs; non-symmetric source-derived element matrix, consistent load assembly, sparse solve.
- **Dynamic mixed formulation:** source-derived non-symmetric 4×4 element mass matrix and assembled second-order system `M q_ddot + K q = f(t)`.
- **Original time scheme:** Crank-Nicolson/trapezoidal startup and weighted-residual three-level recurrence with `gamma=1.5`, `lambda=0.8` (algorithmic damping).
- **New comparison:** Newmark average-acceleration time integration, explicitly **not** claimed as the 2017 method.
- **Loading:** a point force at the beam midpoint, and three historical distributed load shapes (`constant`, `sin(t)`, `source-ramp`) active only during the first second, plus an explicitly **new** unit-height triangular pulse option.
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

The original MATLAB source used inconsistent variable naming for the section width and an element stiffness coefficient. This Python implementation separates these values; see the source audit.

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

On Linux/macOS, create the environment with `python3 -m venv .venv` and activate it with `source .venv/bin/activate`. On Windows choose an installed Python version (3.11 or newer). Existing environments need **not** be modified.

## Numerical verification at a glance

| Claim | Level of evidence |
| --- | --- |
| Mixed static element stiffness | Uses physically consistent `+a` term appearing in the source dynamic matrix; original printed *static* listing has a contradictory `−a`; see Phase 3 audit |
| Static beam validation | Closed-form simply-supported uniform/point-load comparisons |
| Mixed dynamic element mass | Reproduced from original matrix specification |
| Dynamic time solver | Source recurrence reproduced and algebraically tested |
| Dynamic modal validation | Independent FE generalized eigenmode + exact semidiscrete solution, with new Newmark convergence |
| Independent continuum frequency | Exact simply-supported Timoshenko bending frequency 24.7214 Hz; mixed 41-node result 24.7341 Hz |
| Independent physical FE validation | Separate [w, rotation] elements: 0.0168% static midspan difference at 129 nodes; dynamic RMS 0.699% of reference peak over 0.2 s at 41 nodes |
| Stress-recovery audit | Explicit outer-fiber, quarter-thickness, section-average, maximum rectangular and effective Timoshenko shear definitions |
| Historical Abaqus comparisons | **Not independently reproduced**; underlying Abaqus datasets absent |
| Historical static deflection figure | **Still discrepant** with documented inputs |
| Continuum-level dynamic validation | **Incomplete**; future work |

**Source correction note:** The historic static MATLAB listing prints `−a` where the dynamic listing has `+a`. The literal `−a` produces non-convergent solutions; this code uses the verified `+a`. This is a documented correction, **not an exact transcription** of the printed static listing.

**Numerical note:** The original three-level integration scheme exhibits pronounced algorithmic damping at practical step sizes. The new Newmark comparison helps isolate this numerical effect but is not proof that the entire historical physical model is correct. See [detailed audit](docs/DYNAMIC_AUDIT.md).

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

## Interpretation and responsible reuse

The calculated peak linear-elastic normal stresses may exceed the yield strength of some steels. The historic source does not identify a tested grade or its yield strength; the code does not model plasticity or failure. The plots must therefore be interpreted as **linear-elastic predictions**, not evidence of safe service behavior. See the [Phase 3 verification notes](docs/PHASE3_VERIFICATION.md) and [mathematical formulation](docs/MATHEMATICAL_FORMULATION.md).

## Project status and license

This public repository provides a Python reconstruction of a 2017 graduate-level computational mechanics study. Independent analytical and numerical benchmarks support the implementation, while the historical Abaqus datasets and a discrepancy in the original static-deflection figure remain unresolved. Technical limitations and outstanding validation tasks are recorded in the [source audit](docs/SOURCE_AUDIT.md) and [independent verification notes](docs/PHASE3_VERIFICATION.md).

The software is distributed under the [MIT License](LICENSE). The 2026 reconstruction and associated verification work should be distinguished from the original 2017 study when citing results.
