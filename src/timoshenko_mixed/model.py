"""Static mixed (deflection, bending moment) 1D finite-element solver.

The mixed static stiffness matrix follows the documented source formulation.
The printed static MATLAB listing (PDF page 38) has `-a` in its bottom-right
diagonal, whereas the dynamic listing (PDF page 33) has `+a`. The literal
`-a` fails mesh convergence; `+a` matches an independent physical FEM solver.
The independently verified coefficient is used in the 2026 implementation.

Global degrees of freedom: [w_0, M_0, w_1, M_1, ..., w_(n-1), M_(n-1)].
Displacements and moments are mixed primary unknowns.

Loads and supports are expressed explicitly. No commercial FEM solver is
used. Source correspondence is described in docs/SOURCE_AUDIT.md.
"""
from __future__ import annotations

from dataclasses import dataclass
from typing import Literal
import numpy as np
from scipy.sparse import lil_matrix
from scipy.sparse.linalg import spsolve


@dataclass(frozen=True)
class Beam:
    length: float = 2.0
    width: float = 0.03
    height: float = 0.0431
    young_modulus: float = 200e9
    poisson_ratio: float = 0.3
    density: float = 7800.0
    shear_correction: float = 5.0 / 6.0

    def __post_init__(self) -> None:
        if any(x <= 0 for x in (
            self.length, self.width, self.height, self.young_modulus,
            self.density, self.shear_correction
        )):
            raise ValueError("Beam lengths and physical parameters must be positive")
        if not (-1.0 < self.poisson_ratio < 0.5):
            raise ValueError("Poisson ratio must lie between -1 and 0.5")

    @property
    def area(self) -> float:
        return self.width * self.height

    @property
    def second_moment(self) -> float:
        return self.width * self.height**3 / 12.0

    @property
    def shear_modulus(self) -> float:
        return self.young_modulus / (2.0 * (1.0 + self.poisson_ratio))

    @property
    def EI(self) -> float:
        return self.young_modulus * self.second_moment

    @property
    def kGA(self) -> float:
        return self.shear_correction * self.shear_modulus * self.area


@dataclass
class StaticSolution:
    nodes: np.ndarray
    deflection: np.ndarray
    bending_moment: np.ndarray
    free_residual_l2: float
    boundary: str
    distributed_intensity: float
    distributed_end: float
    point_force: float
    point_position: float


def element_stiffness(beam: Beam, element_length: float) -> np.ndarray:
    """Return the corrected 4 x 4 source-aligned mixed stiffness matrix.

    Local DOFs: [w_i, M_i, w_j, M_j]. Matrix is non-symmetric by formulation.
    Bottom-right `+a` is a correction of the printed static `-a` inconsistency
    supported by the historical dynamic listing and independent verification.
    """
    if element_length <= 0:
        raise ValueError("element_length must be positive")
    a = element_length**2 / (3.0 * beam.EI) + 1.0 / beam.kGA
    b = element_length**2 / (6.0 * beam.EI) - 1.0 / beam.kGA
    return np.array([
        [0.0,  1.0, 0.0, -1.0],
        [-1.0, a,  1.0, b],
        [0.0, -1.0, 0.0,  1.0],
        [1.0, b, -1.0,  a],
    ], dtype=float) / element_length


def assemble_stiffness(beam: Beam, n_nodes: int):
    """Assemble a sparse, global static mixed matrix."""
    if n_nodes < 3:
        raise ValueError("Use at least 3 nodes")
    nodes = np.linspace(0, beam.length, n_nodes)
    K = lil_matrix((2 * n_nodes, 2 * n_nodes), dtype=float)
    for e in range(n_nodes - 1):
        dofs = np.arange(2 * e, 2 * e + 4)
        Ke = element_stiffness(beam, nodes[e + 1] - nodes[e])
        for i in range(4):
            for j in range(4):
                K[dofs[i], dofs[j]] += Ke[i, j]
    return nodes, K.tocsr()


def _assemble_loads(nodes: np.ndarray, q: float, q_end: float,
                    point_force: float, point_x: float) -> np.ndarray:
    """Consistent linear distributed-load vector + interpolated point load.

    A partial element is integrated exactly for constant q. This is more
    general than the report's half-beam loop, and works for even node counts.
    """
    n = len(nodes)
    L = float(nodes[-1])
    if not 0 <= q_end <= L:
        raise ValueError("distributed_end must lie in [0, length]")
    if not 0 <= point_x <= L:
        raise ValueError("point_position must lie in [0, length]")
    f = np.zeros(2 * n, dtype=float)
    for e in range(n - 1):
        x0 = float(nodes[e]); x1 = float(nodes[e + 1])
        hi = min(x1, q_end)
        if hi <= x0 or q == 0:
            continue
        le = x1 - x0
        # Integral N1 dx, integral N2 dx on [x0, hi]
        int_n2 = (hi - x0)**2 / (2.0 * le)
        int_n1 = (hi - x0) - int_n2
        f[2 * e] += q * int_n1
        f[2 * (e + 1)] += q * int_n2
    if point_force != 0:
        j = min(int(np.searchsorted(nodes, point_x, side="right")), n - 1)
        if j == 0:
            f[0] += point_force
        else:
            x0, x1 = nodes[j - 1], nodes[j]
            frac = (point_x - x0) / (x1 - x0)
            f[2 * (j - 1)] += point_force * (1.0 - frac)
            f[2 * j] += point_force * frac
    return f


def solve_static(
    n_nodes: int = 55,
    beam: Beam | None = None,
    boundary: Literal["fixed-pinned", "pinned-pinned"] = "fixed-pinned",
    distributed_intensity: float = 4800.0,
    distributed_end: float | None = None,
    point_force: float = 3000.0,
    point_position: float | None = None,
) -> StaticSolution:
    """Solve the static mixed system.

    Source study: fixed-pinned, q=4800 N/m on [0, L/2], central Q=3000 N.
    Fixed-pinned mixed DOFs: w(0)=0, w(L)=0, M(L)=0.
    Pinned-pinned mixed DOFs: w(0)=M(0)=w(L)=M(L)=0.
    The first case's zero rotation at x=0 is a natural condition associated
    with the weak form; it is not explicitly a rotation DOF.
    """
    b = beam or Beam()
    q_end = b.length / 2 if distributed_end is None else distributed_end
    x_point = b.length / 2 if point_position is None else point_position
    nodes, K = assemble_stiffness(b, n_nodes)
    f = _assemble_loads(nodes, distributed_intensity, q_end, point_force, x_point)

    if boundary == "fixed-pinned":
        constrained = np.array([0, 2 * n_nodes - 2, 2 * n_nodes - 1], dtype=int)
    elif boundary == "pinned-pinned":
        constrained = np.array([0, 1, 2 * n_nodes - 2, 2 * n_nodes - 1], dtype=int)
    else:
        raise ValueError("Unsupported boundary condition")
    free_mask = np.ones(2 * n_nodes, dtype=bool)
    free_mask[constrained] = False
    free = np.flatnonzero(free_mask)
    d = np.zeros(2 * n_nodes, dtype=float)
    d[free] = spsolve(K[free, :][:, free], f[free])
    residual = K @ d - f
    return StaticSolution(
        nodes=nodes, deflection=d[::2], bending_moment=d[1::2],
        free_residual_l2=float(np.linalg.norm(residual[free])),
        boundary=boundary, distributed_intensity=distributed_intensity,
        distributed_end=q_end, point_force=point_force, point_position=x_point,
    )


def midpoint_deflection(solution: StaticSolution) -> float:
    return float(np.interp(solution.nodes[-1] / 2.0, solution.nodes, solution.deflection))
