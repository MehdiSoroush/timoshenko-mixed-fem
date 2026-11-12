# Mathematical formulation and reproducibility guide

This English description accompanies a **2026 Python reconstruction** of a
**2017 Timoshenko-beam mixed FEM graduate course project**. The archival 2017 source document is not distributed. Source-specific uncertainties are documented in
`SOURCE_AUDIT.md` and `PHASE3_VERIFICATION.md`.

## Beam model and variables

Length `L`, cross-sectional area `A=b*h`, rectangular section second moment
`I=b*h^3/12`, Young's modulus `E`, shear modulus `G=E/[2(1+nu)]`, shear
factor `ks=5/6`, density `rho`.

The original mixed model solves for transverse displacement `w(x,t)` and a
bending-moment DOF `M(x,t)`. The separate validation model uses `w(x,t)` and
cross-section rotation `phi(x,t)`.

For the reference solver, the shear strain is `gamma = dw/dx - phi`, the
transverse shear force is `V = ks*G*A*gamma`, and the conventional bending
moment is `M_physical = E*I*d(phi)/dx`. The original mixed DOF is found to have
the **opposite sign** to `M_physical` under the implementation's load/rotation
sign conventions. Thus `M_mixed ≈ -E*I*phi'`, and `V ≈ d(M_mixed)/dx`.

The static continuum equations with these conventions include
`-d²(M_mixed)/dx² = q(x)` and
`w'' + M_mixed/(EI) - M_mixed''/(ks*G*A) = 0`, away from concentrated forces.
A point force creates a jump in shear `V` across its application point.

## Source-derived mixed static element

For a two-node element of length `le` with ordering
`[w_i, M_i, w_j, M_j]`, define

```text
a = le**2/(3*EI) + 1/(ks*G*A)
beta = le**2/(6*EI) - 1/(ks*G*A)

Ke = 1/le * [[ 0,  1,    0, -1],
             [-1,  a,    1, beta],
             [ 0, -1,    0,  1],
             [ 1, beta, -1,  a]]
```

This mixed local matrix is nonsymmetric. The 2017 *static* printed listing
shows `-a` in the last position; the 2017 *dynamic* listing has `+a`. The code
uses the physically convergent **corrected `+a`** and independently verifies it.
The printed `-a` fails numerical convergence and must **not** be silently
presented as the implemented source matrix.

The distributed transverse loading is assembled by exact integration of
linear shape functions over the loaded portions of elements. The point load
is assigned by nodal interpolation when necessary.

**Original support conditions:** a clamped left and pinned right end,
`w(0)=0`, `phi(0)=0`, `w(L)=0`, `M_mixed(L)=0`. In the mixed weak formulation,
`phi(0)=0` is a natural boundary condition and no explicit rotation DOF exists.
The constrained mixed DOFs are therefore `w(0)`, `w(L)`, and `M(L)`.
For the simply-supported benchmark: `w(0)=M(0)=w(L)=M(L)=0`.

## Source-derived mixed mass matrix

The 2017 source defines `c=rho*A`, `d=rho/E`, and
`g=rho/(ks*G*E*A)` and uses

```text
Me = le * [[c/3, d/3, c/6, d/6],
           [  0, g/3,   0, g/6],
           [c/6, d/6, c/3, d/3],
           [  0, g/6,   0, g/3]]
```

The semidiscrete mixed system is `Mglobal * u_ddot + Kglobal * u = f(t)`.
Both assembled mixed matrices are generally nonsymmetric. The original
weighted-residual three-level recurrence and its Crank–Nicolson startup
are documented in `DYNAMIC_AUDIT.md`. An additional Newmark average-
acceleration integrator is marked as **new 2026 verification work**, not
part of the original coursework.

## Independent physical displacement–rotation FEM

The independent 2026 reference uses
`[w_i, phi_i, w_j, phi_j]` per two-node element. Define

```text
B_shear = [-1/le, -1/2,  1/le, -1/2]
B_curv  = [    0, -1/le,     0,  1/le]
Ke_ref  = (ks*G*A)*le*(B_shear.T @ B_shear)
          + (E*I)*le*(B_curv.T @ B_curv)
```

The shear term uses **one-point integration**, preventing the spurious
shear-locking response of fully integrated linear Timoshenko elements in
slender beams. Bending is exact for linear rotations.

The physical mass includes consistent transverse mass `rho*A` on the `w`
DOFs and rotary mass `rho*I` on the `phi` DOFs; for each pair of nodes, both
use `le/6*[[2,1],[1,2]]` times the corresponding mass density coefficient.
No source-derived mixed matrices are reused in the physical reference model.

## Independent checks and formulas

For a simply-supported beam of length `L`, the midspan deflections are:

- Full-length uniform load `q`:
  `5*q*L**4/(384*EI) + q*L**2/(8*ks*G*A)`.
- Center point load `Q`:
  `Q*L**3/(48*EI) + Q*L/(4*ks*G*A)`.

For simply-supported free vibration with rotary inertia and `k=n*pi/L`,
the lower bending-branch frequency is obtained by solving

```text
[(ks*G*A)*k² - (rho*A)*omega²]
  *[(E*I)*k² + (ks*G*A) - (rho*I)*omega²]
  - [(ks*G*A)*k]² = 0
```

The repository uses a numerically stable quadratic formula for the smaller
positive `omega²`. See `analytical.py` and the independent tests.

## Stress recovery and interpretation

At an element midpoint, approximate `M_mid=(M_i+M_j)/2` and the source-signed
shear `V=(M_j-M_i)/le`. Linear bending gives `sigma_x=M_mid*z/I` in the
source's sign convention; the sign must be inverted to compare with a physical
moment definition `EI*phi'` depending on face orientation.

- Surface (`z=h/2`) normal stress: `M_mid*h/(2I)`.
- Historical normal stress formula: `M_mid*h/(4I)`, corresponding to `z=h/4`.
- Section-average transverse shear stress: `V/A`.
- Rectangular-section peak (classical parabolic distribution): `1.5*V/A`.
- Historical effective shear measure: `V/(ks*A)`.

Because shear correction converts a nonuniform 3D cross-sectional field into
an equivalent beam energy, the historical effective shear measure must **not**
be described as the actual maximum transverse stress.

## Limitations

All computed examples use small-strain *linear elasticity*. A predicted peak
stress exceeding the yield stress of a particular grade would invalidate a
purely linear-elastic interpretation, but the historical source does not name
a tested steel grade or its yield strength. No experimental strength claims
are made. The original Abaqus results and source MATLAB runtime are not
available for definitive historical reproduction; independent analytical and
physical finite-element validations are used instead.
