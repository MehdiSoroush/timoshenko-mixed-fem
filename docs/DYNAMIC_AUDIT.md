# Phase 2: Dynamic mixed FEM reconstruction and verification

**Original project:** 2017 (Tir 1396, Iranian calendar).  
**Current implementation:** Python reconstruction, November 2026.  
**Provenance:** Archival 2017 coursework report remains private; it is **not distributed in this repository**.

## What the report actually specifies

The 2017 report's *PDF pages* (1-indexed) contain:

| PDF pages | Source content | Status of reconstruction |
| --- | --- | --- |
| 9-10 | Beam geometry, three time-varying distributed load cases | Implemented |
| 26-27 | Crank-Nicolson startup, three-point weighted-residual time recurrence, \(\gamma,\lambda\) table | Implemented |
| 32-33 | Material quantities, element stiffness and mixed mass matrices | Implemented directly |
| 34 | Boundary and loading assembly | Implemented, with more robust DOF elimination |
| 35 | MATLAB Crank-Nicolson startup, three-point recurrence | Implemented directly in algebraic form |
| 36-37 | Postprocessing formulas and plotting | Deflection and moments in Phase 2; physically labeled element stress checks now implemented in Phase 3 |
| 45-63 | Historical time convergence, numerical deflections, moment plots, Abaqus comparisons | Figures reviewed; **not numerically reproduced/validated** |

## Source-derived equation

The mixed semidiscrete unknowns alternate transverse displacement and bending moment:

`q = [w_0, M_0, w_1, M_1, ..., w_N, M_N]`

The original dynamic system is

`M q_ddot + K q = f(t)`.

The mixed element mass matrix (PDF p. 33) is nonsymmetric. With `c=rho*A`, `d=rho/E`, `g=rho/(k_s*G*E*A)`:

```
Me = le * [[c/3, d/3, c/6, d/6],
           [0,   g/3, 0,   g/6],
           [c/6, d/6, c/3, d/3],
           [0,   g/6, 0,   g/3]]
```

The implementation does **not** silently replace it with a symmetric beam mass matrix. The first physical modes of the generalized `(K, M)` matrix pencil are verified with a separate generalized-eigenvalue computation, but that is *not* sufficient to independently establish continuum accuracy of the full dynamic model.

### Initial step

For initial displacement `u0` and velocity `v0`, the historical Crank-Nicolson startup can be written (zero damping):

`(M + dt² K / 4) u1 = (M - dt² K / 4) u0 + dt M v0 + dt² (f0 + f1) / 4`

The source program uses zero initial displacement and velocity. The Python implementation also accepts consistent nonzero initial states for modal verification.

### Subsequent steps

The report chooses `gamma=3/2`, `lambda=4/5`. Defining

```
alpha = 1/2 - 2*lambda + gamma
beta  = 1/2 + lambda - gamma
A = M + lambda*dt²*K
B = 2*M - alpha*dt²*K
C = -M - beta*dt²*K
```

the resulting recurrence is

`A u[j+1] = B u[j] + C u[j-1] + dt²*(beta*f[j-1] + alpha*f[j] + lambda*f[j+1])`.

At each time step, the sparse linear system is solved using a reused LU factorization; no explicit matrix inverse is formed.

## Load-profile caveats

The central `Q=3000 N` point force is held on for the full simulation. The distributed load of `P0=4800 N/m` over the **left half** is active only during the original `L_D=1 s` interval:

1. `constant`: `P(t)/P0 = 1` while active.
2. `sine`: `P(t)/P0 = sin(t)` while active, using the original report's **1 rad/s** time dependence rather than an inferred frequency.
3. `source-ramp`: linearly grows and decays with an actual peak **0.5**, not 1.0. The historical discrete formula uses `floor(1/dt)` and slightly changes the profile for arbitrary step sizes. This reconstruction uses its continuous-time counterpart with exact peak 0.5 at the half-second mark when `duration=1`.

A *new*, explicitly labeled `unit-triangle` profile peaks at 1.0 and is available only for demonstration. Do not conflate it with the 2017 coursework code.

## Verification results and what they do **not** establish

- The element mass coefficients have an exact-code unit test against the printed source.
- The translational mass submatrix sums to `rho*A*L` as expected for a uniform field.
- All three historical forcing functions produce finite results and preserve displacement/moment boundary constraints.
- The source three-level equation and the Newmark equations satisfy small numerical residuals at each step (these are **algebraic checks**, not evidence of physical validation).
- An *independently computed FE generalized eigenvector* gives an exact scalar modal reference `q(t) = phi*cos(omega*t)` for **the discretized problem**. A newly added Newmark average-acceleration integrator demonstrates second-order temporal convergence against this solution.
- **Significant algorithmic damping** is present in the original three-level scheme for the first mode at the selected timesteps. It may strongly influence the original dynamic graphs. The new Newmark integrator is a comparison, **not** a replacement claimed to be present in the historical work.
- The first simply-supported natural frequency approaches the expected order of magnitude from the Euler-Bernoulli slender-beam expression. That expression is only an approximate reference for the Timoshenko beam.

**Updated after Phase 3:** The independent physical [w, rotation] FE transient and exact continuum first bending frequency are now verified in selected cases; physically labeled dynamic stress recovery is implemented. Important remaining caveats: no original Abaqus data files; original dynamic plots not quantitatively digitized; additional physical loading/time benchmarks and original three-level damping studies are warranted; original static figure discrepancy remains unresolved. See [Phase 3 verification](PHASE3_VERIFICATION.md).

## Running

```bash
python examples/reproduce_dynamic_cases.py
python examples/verify_free_vibration.py
python -m pytest -q
```

Generated figures and CSV summaries are **new results of this reconstruction**, not original experimental data or recovered graphs.
