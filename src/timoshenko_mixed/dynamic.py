"""Source-derived *dynamic* mixed beam finite elements, with an audited time scheme.

Historical source: 2017 coursework report, PDF pp. 32--37 (1-indexed).
This module is a new Python reconstruction, not recovered MATLAB code.

The original equations are M u_ddot + K u = f(t), with mixed unknowns
[w_0, M_0, w_1, M_1, ...]. Both M and K are non-symmetric in this
formulation. The reported 3-level, weighted-residual recurrence is preserved
with gamma=3/2 and lambda=4/5, with a Crank--Nicolson startup (equivalent
to the initial displacement step of trapezoidal/Newmark average acceleration).
The scheme has algorithmic high-frequency damping; it is not energy preserving.

NOTE: The historical report includes no raw Abaqus numerical data and does not
provide a quantified dynamic accuracy assessment. See docs/DYNAMIC_AUDIT.md.
"""
from __future__ import annotations

from dataclasses import dataclass
from typing import Literal

import numpy as np
from scipy.sparse import lil_matrix
from scipy.sparse.linalg import splu

from .model import Beam, _assemble_loads, assemble_stiffness

DynamicLoad = Literal["constant", "sine", "source-ramp", "unit-triangle"]
Boundary = Literal["fixed-pinned", "pinned-pinned"]


@dataclass
class DynamicSolution:
    nodes: np.ndarray
    times: np.ndarray
    deflection: np.ndarray  # (n_times, n_nodes), in metres
    bending_moment: np.ndarray  # (n_times, n_nodes), in N*m
    boundary: str
    load_case: str
    load_duration: float
    distributed_intensity: float
    point_force: float
    startup_residual: float
    maximum_step_residual: float
    _all_dofs: np.ndarray  # (n_times, 2*n_nodes)
    integrator: str

    def midpoint_deflection(self) -> np.ndarray:
        """Interpolate deflection at the beam midpoint at each time step."""
        j = np.searchsorted(self.nodes, self.nodes[-1] / 2, side="right")
        j = min(max(j, 1), len(self.nodes) - 1)
        x0, x1 = self.nodes[j-1], self.nodes[j]
        p = ((self.nodes[-1] / 2) - x0) / (x1-x0)
        return (1-p)*self.deflection[:, j-1] + p*self.deflection[:, j]


def element_mass(beam: Beam, element_length: float) -> np.ndarray:
    """Reproduce the non-symmetric 4x4 mixed mass matrix in source PDF p. 33.

    The mass coefficient g is rho/(k_s*G*E*A), precisely as printed. It has
    a very different magnitude from c=rho*A, so nonsymmetric scaling and
    conditioning should be monitored carefully. Do not symmetrize this matrix.
    """
    if element_length <= 0:
        raise ValueError("element_length must be positive")
    c = beam.density * beam.area
    d = beam.density / beam.young_modulus
    g = beam.density / (beam.kGA * beam.young_modulus)
    return element_length * np.array([
        [c/3, d/3, c/6, d/6],
        [0,   g/3, 0,   g/6],
        [c/6, d/6, c/3, d/3],
        [0,   g/6, 0,   g/3],
    ], dtype=float)


def assemble_mass(beam: Beam, n_nodes: int):
    if n_nodes < 3:
        raise ValueError("Use at least three nodes")
    nodes = np.linspace(0, beam.length, n_nodes)
    M = lil_matrix((2*n_nodes, 2*n_nodes), dtype=float)
    for e in range(n_nodes-1):
        dofs = np.arange(2*e, 2*e+4)
        Me = element_mass(beam, nodes[e+1]-nodes[e])
        for r in range(4):
            for c in range(4):
                M[dofs[r], dofs[c]] += Me[r, c]
    return nodes, M.tocsr()


def load_multiplier(t: float, case: DynamicLoad, duration: float = 1.0) -> float:
    """Historical distributed-load time profiles (point load stays on).

    `sine` uses sin(t) with t in seconds, as written in the original source.
    `source-ramp` peaks at 0.5 for duration=1, reflecting the reported code;
    `unit-triangle` is an explicitly NEW, peak-1 alternative. The source's
    integer-floor discretization is approximated by continuous piecewise time.
    """
    if duration <= 0:
        raise ValueError("duration must be positive")
    if case not in ("constant", "sine", "source-ramp", "unit-triangle"):
        raise ValueError(f"Unsupported dynamic load case: {case}")
    if t < 0 or t > duration:
        return 0.0
    if case == "constant":
        return 1.0
    if case == "sine":
        return float(np.sin(t))
    half_triangle = min(t/duration, 1.0-t/duration)
    return half_triangle if case == "source-ramp" else 2.0*half_triangle


def _free_dofs(n_nodes: int, boundary: Boundary) -> np.ndarray:
    if boundary == "fixed-pinned":
        fixed = [0, 2*n_nodes-2, 2*n_nodes-1]
    elif boundary == "pinned-pinned":
        fixed = [0, 1, 2*n_nodes-2, 2*n_nodes-1]
    else:
        raise ValueError(f"Unsupported boundary: {boundary}")
    mask = np.ones(2*n_nodes, dtype=bool)
    mask[fixed] = False
    return np.flatnonzero(mask)


def solve_dynamic(
    n_nodes: int = 21,
    *,
    beam: Beam | None = None,
    boundary: Boundary = "fixed-pinned",
    load_case: DynamicLoad = "constant",
    duration: float = 1.0,
    t_end: float = 2.0,
    dt: float = 0.005,
    distributed_intensity: float = 4800.0,
    distributed_end: float | None = None,
    point_force: float = 3000.0,
    point_position: float | None = None,
    initial_dofs: np.ndarray | None = None,
    initial_velocity: np.ndarray | None = None,
    gamma: float = 1.5,
    lam: float = 0.8,
    integrator: Literal["source-three-level", "newmark-average"] = "source-three-level",
) -> DynamicSolution:
    """Solve the original mixed semidiscrete dynamics using its 3-level method.

    Source convention: distributed loading on the left half, applied only
    for the first `duration=1 s`; a central point load Q remains applied.
    Zero displacement and velocity are assumed unless supplied. Initial DOF
    arrays use alternating [w,M] ordering and must satisfy constraints.

    The source's constant, sin(t), and half-height triangle are all available.
    A fully normalized triangle is an explicitly *new* comparison case.
    Constraints are imposed by free-DOF elimination rather than row stamping.
    """
    b = beam if beam is not None else Beam()
    if not np.isfinite(t_end) or not np.isfinite(dt) or t_end <= 0 or dt <= 0:
        raise ValueError("t_end and dt must be finite and positive")
    n_steps = int(round(t_end/dt))
    if n_steps < 2 or not np.isclose(n_steps*dt, t_end, rtol=1e-10, atol=1e-12):
        raise ValueError("t_end/dt must be an integer >= 2")
    if duration <= 0 or not np.isfinite(duration):
        raise ValueError("duration must be finite and positive")
    if integrator not in ("source-three-level", "newmark-average"):
        raise ValueError(f"Unsupported integrator: {integrator}")
    # The source scheme has only one reported (gamma, lambda) parameter pair.
    if integrator == "source-three-level" and not (np.isclose(gamma, 1.5) and np.isclose(lam, 0.8)):
        raise ValueError("The source recurrence uses gamma=1.5, lambda=0.8")
    x, K = assemble_stiffness(b, n_nodes)
    xm, M = assemble_mass(b, n_nodes)
    np.testing.assert_allclose(x, xm)
    free = _free_dofs(n_nodes, boundary)
    Kf = K[free, :][:, free].tocsc()
    Mf = M[free, :][:, free].tocsc()
    tf = b.length/2 if distributed_end is None else distributed_end
    xp = b.length/2 if point_position is None else point_position
    unit_q = _assemble_loads(x, 1.0, tf, 0.0, xp)[free]
    point = _assemble_loads(x, 0.0, tf, point_force, xp)[free]
    times = np.linspace(0.0, t_end, n_steps+1)
    scales = np.asarray([load_multiplier(t, load_case, duration) for t in times])
    forces = point[None, :] + (distributed_intensity * scales[:, None]) * unit_q[None, :]
    u = np.zeros((n_steps+1, len(free)), dtype=float)
    v0 = np.zeros(2*n_nodes) if initial_velocity is None else np.asarray(initial_velocity, dtype=float)
    u0 = np.zeros(2*n_nodes) if initial_dofs is None else np.asarray(initial_dofs, dtype=float)
    if v0.shape != (2*n_nodes,) or u0.shape != (2*n_nodes,):
        raise ValueError("Initial DOF arrays must have length 2*n_nodes")
    if not np.isfinite(v0).all() or not np.isfinite(u0).all():
        raise ValueError("Initial DOFs must be finite")
    fixed_mask = np.ones(2*n_nodes, dtype=bool)
    fixed_mask[free] = False
    if np.any(u0[fixed_mask] != 0) or np.any(v0[fixed_mask] != 0):
        raise ValueError("Initial DOFs must satisfy constrained boundary conditions")
    u[0] = u0[free]
    vf0 = v0[free]
    dt2 = dt*dt

    if integrator == "source-three-level":
        # Equation (48), p.26, after multiplying both sides by dt/2.
        # The source's Crank--Nicolson startup uses u0, v0 and force at
        # t0 and t1. Unlike a first-order rule, it is consistent with
        # second-order dynamics when initial velocity is included.
        starter = (Mf + .25*dt2*Kf).tocsc()
        starter_lu = splu(starter)
        rhs1 = (Mf - .25*dt2*Kf) @ u[0] + dt*(Mf @ vf0) + (.25*dt2)*(forces[0]+forces[1])
        u[1] = starter_lu.solve(rhs1)
        startup_residual = float(np.linalg.norm(starter@u[1]-rhs1, ord=np.inf))

        # Equation (49), source PDF p.26; implementation on PDF p.35.
        alpha = .5 - 2*lam + gamma
        beta = .5 + lam - gamma
        A = (Mf + lam*dt2*Kf).tocsc()
        B = 2*Mf - alpha*dt2*Kf
        C = -Mf - beta*dt2*Kf
        lu = splu(A)
        max_res = 0.0
        for j in range(1, n_steps):
            rhs = (B@u[j] + C@u[j-1] + dt2*(beta*forces[j-1] + alpha*forces[j] + lam*forces[j+1]))
            u[j+1] = lu.solve(rhs)
            res = np.linalg.norm(A@u[j+1]-rhs, ord=np.inf)
            max_res = max(max_res, float(res))
    else:
        # NEW: Standard Newmark average acceleration (beta=1/4, gamma=1/2)
        # used as an independent time-integration comparison. For an
        # undamped stable semidiscrete system it is second-order accurate
        # and avoids the original algorithm's artificial modal damping.
        beta_nm, gamma_nm = 0.25, 0.5
        acceleration = np.zeros_like(u)
        velocity = np.zeros_like(u)
        velocity[0] = vf0
        mass_lu = splu(Mf)
        acceleration[0] = mass_lu.solve(forces[0] - Kf@u[0])
        A = (Mf + beta_nm*dt2*Kf).tocsc()
        lu = splu(A)
        max_res = 0.0
        for j in range(n_steps):
            u_predict = u[j] + dt*velocity[j] + dt2*(0.5-beta_nm)*acceleration[j]
            v_predict = velocity[j] + dt*(1-gamma_nm)*acceleration[j]
            rhs = forces[j+1] - Kf@u_predict
            a_next = lu.solve(rhs)
            acceleration[j+1] = a_next
            u[j+1] = u_predict + beta_nm*dt2*a_next
            velocity[j+1] = v_predict + gamma_nm*dt*a_next
            max_res = max(max_res, float(np.linalg.norm(Mf@a_next + Kf@u[j+1] - forces[j+1], ord=np.inf)))
        startup_residual = float(np.linalg.norm(Mf@acceleration[0] + Kf@u[0] - forces[0], ord=np.inf))
    full = np.zeros((n_steps+1, 2*n_nodes))
    full[:, free] = u
    return DynamicSolution(
        nodes=x, times=times, deflection=full[:, ::2],
        bending_moment=full[:, 1::2], boundary=boundary, load_case=load_case,
        load_duration=duration, distributed_intensity=distributed_intensity,
        point_force=point_force, startup_residual=startup_residual,
        maximum_step_residual=max_res, _all_dofs=full, integrator=integrator,
    )
