# Phase 3: independent continuum/reference FEM and stress verification

**Original project:** 2017 (Tir 1396); **reconstruction and new checks:** November 2026.
The archival 2017 coursework report is intentionally **not distributed** here.
These are newly computed checks, not archived Abaqus files, recovered numerical results, or experimental data.

## Source audit: an important static-matrix inconsistency

The 2017 manuscript writes the dynamic local matrix with `+a` in the last
stiffness diagonal (PDF page 33), but the *static MATLAB listing* writes `-a`
there (PDF page 38). These expressions are not identical. Literal use of the
printed `-a` static entry yields non-convergent results: for 55 nodes the
midspan deflection is about **−2.862 mm**, worsening to **−19.072 mm** at 101
nodes. The consistent `+a` matrix converges to about **+9.562 mm**, and agrees
with an entirely separate physical [deflection, rotation] FEM solver. The
Python implementation uses `+a`. This is a **documented correction/interpretation**,
not a literal copy of the source static line.

See `examples/audit_historical_static.py` and the generated
`results/data/static_matrix_sign_audit.csv`. The archived MATLAB listing also
contains a capitalization collision (`B` defines the section width before `b`
is used) and uses `inv(Kt)*Ft`, both addressed in the reconstruction.

**Still unresolved:** The static deflection plot in the original report (PDF
page 42) peaks near 19 mm, whereas consistent source parameters and two
independent current implementations yield roughly 9.66 mm **maximum** deflection
(9.56 mm at the midpoint). The discrepancy cannot be resolved
without the original runnable files, numerical Abaqus output, or additional
parameters. The `-a` typo does not satisfactorily explain this figure.

## Independent conventional FE model

A completely separate solver was added in `src/timoshenko_mixed/reference.py`:

- Nodal DOFs `[w_i, phi_i, w_j, phi_j]`, rather than the original mixed
  `[w_i, M_i, w_j, M_j]`.
- Consistent transverse and rotary mass (`rho*A` and `rho*I`).
- Standard bending energy `EI phi'^2` and shear energy
  `k_s*G*A*(w'−phi)^2`, with one-point shear integration to avoid locking.
- An independently assembled spatial load vector and standard Newmark
  average-acceleration temporal integration.
- Physical fixed–pinned conditions `w(0)=phi(0)=w(L)=0`, and simply-supported
  conditions `w(0)=w(L)=0`.

These are mathematically distinct spatial discretizations. Numerical agreement
is a meaningful check of the original mixed formulation rather than a
self-comparison of two time-stepping routines applied to the same matrices.

### Static convergence, original fixed–pinned loading

`Q=3000 N` at midpoint, `q=4800 N/m` on the left half, baseline material and
geometry. Midspan deflections (all newly computed):

| Nodes | Mixed [w,M], mm | Independent [w,phi], mm | Difference / reference |
| ---: | ---: | ---: | ---: |
| 9 | 9.581376 | 9.170208 | 4.4837% |
| 17 | 9.566786 | 9.463910 | 1.0870% |
| 33 | 9.563139 | 9.537415 | 0.2697% |
| 65 | 9.562227 | 9.555796 | 0.0673% |
| 129 | 9.561999 | 9.560391 | **0.0168%** |

These results support static convergence of both methods; they do not
establish agreement with the absent original Abaqus model.

### Exact simply supported continuum frequency

The frequency is obtained from the **Timoshenko continuum dispersion relation**
including rotary inertia, with wavenumber `k = pi/L` and bending branch:

`(kGA*k² − rho*A*omega²)*(EI*k² + kGA − rho*I*omega²) − (kGA*k)² = 0`

The analytical first bending frequency is **24.7213897 Hz**. The 41-node
original mixed generalized eigenproblem gives **24.7340902 Hz**, about
**0.0514% above** the continuum result. The separate physical [w,phi] FE model
gives **24.7404431 Hz** with 41 nodes, about **0.0771% above**. These
comparisons are more meaningful than a rough Euler–Bernoulli frequency check.

### Independent dynamic response comparison

For the fixed–pinned beam, 41 nodes, constant distributed load during the first
second (only the first 0.2 s evaluated), a 3000 N central force, zero initial
conditions, and a timestep of `dt=0.001 s`:

- Both solvers use **Newmark average acceleration** to isolate the spatial
  and mass-model differences; the original three-level time scheme is *not*
  being compared here.
- The midspan displacement RMS difference is **0.132761 mm**, or **0.699%**
  of the independent physical FEM response's peak magnitude.
- This checks the short-duration transient of the documented physical beam,
  not the full 2017 dynamic Abaqus history, which is unavailable.

The independent reference has its own stiffness, mass, force assembly, and
boundary elimination; only the shared scalar temporal load profile is reused.
The agreement is not proof that every high-frequency mode or all original
three-level algorithm outputs are free of issues. The source three-level method
was independently assessed for numerical damping in `DYNAMIC_AUDIT.md`.

## Stress recovery: explicit conventions and an important factor of two

The original stress-recovery listing (PDF pages 36, 40) computes
`M_mid*h/(4*I)` after multiplication by E. This value corresponds to the
bending stress at a position `|z|=h/4` for linear stress variation, not to the
surface stress at `|z|=h/2`. It is **one-half** the outer-fiber magnitude.
The source does not unambiguously say which through-thickness location was
intended; therefore both conventions are retained rather than declaring the
historic result definitively incorrect.

`stress.py` explicitly provides:

- `sigma_outer = M_mid*h/(2I)`, outer-fiber source-signed normal stress.
- `sigma_source_2017 = M_mid*h/(4I)`, literal historical quantity.
- `shear_force = (M_j−M_i)/element_length`.
- `tau_section_average = V/A`.
- `tau_rectangular_max = 3V/(2A)` for a rectangular section under classical
  parabolic Saint-Venant shear distribution.
- `tau_effective_timoshenko = V/(k_s*A)` corresponding to the 2017
  `G*epsilon_13` effective shear measure; **not** a unique section stress.

Signs: the mixed method's bending-moment DOF has the opposite sign to
`EI*d(phi)/dx` in the independent physical reference (at element midpoints).
Their magnitudes agree closely; the shear forces recovered from the nodal
moment gradient agree with the independent reference's `kGA*(w'−phi)`.
Sign conventions must be adjusted before comparing labeled tensile/compressive
faces to Abaqus outputs.

For the baseline fixed–pinned static case, at 101 nodes:

| Computed quantity | Maximum absolute value |
| --- | ---: |
| Outer-fiber bending normal stress | 259.585 MPa |
| 2017 source normal stress at quarter-thickness | 129.792 MPa |
| Average section shear stress | 4.864 MPa |

A representative 21-node, `dt=.001 s`, first-0.2-s constant-load Newmark
solution reaches about **467.993 MPa** peak *linear-elastic* outer-fiber bending
stress. These calculations presume linear elastic material behavior; the
report does not establish a material grade or yield stress. Therefore these
numbers must not be interpreted as evidence that the beam would remain elastic
under that loading in practice.

## Reproduce

```powershell
python -m pip install -e ".[dev]"
python -m pytest -q
python examples/audit_historical_static.py
python examples/validate_independent_reference.py
python examples/recover_stresses.py
```

Generated figures include:

- `results/figures/static_matrix_sign_audit.png`
- `results/figures/independent_fem_validation.png`
- `results/figures/independent_frequency_verification.png`
- `results/figures/verified_stress_recovery.png`

Supporting CSV files are in `results/data/`.

## Remaining work before describing this as a fully reproduced historical study

1. Reconcile the original 2017 static graph (apparently ~19 mm) with the
   actual documented input parameters and now cross-verified ~9.56 mm midpoint.
2. If historical Abaqus input/output data can be recovered, compare like-for-like
   geometry, formulation, mass, loading, supports, and time sampling.
3. Verify additional dynamic benchmark cases and source three-level numerical
   damping across a range of steps and modes. Compare physical shear/bending
   definitions with intended plotting locations.

The available evidence supports describing this repository as a **verified Python reconstruction of a 2017 computational mechanics study**, with explicitly documented modifications and independent reference calculations. It does not establish exact reproduction of the historical Abaqus results.
