"""Independent displacement--rotation Timoshenko beam reference solver.

This is NEW 2026 verification work, deliberately unlike the source-derived
mixed [w, M] element.  Two-node [w, phi] elements use exact bending
integration and *one-point* shear integration to avoid shear locking.

The consistent physical mass includes transverse inertia rho*A and rotary
inertia rho*I.  The historical mixed dynamic matrix uses a different set of
variables and must not be assumed equivalent without numerical comparison.
"""
from __future__ import annotations
from dataclasses import dataclass
from typing import Literal
import numpy as np
from scipy.linalg import eigh
from scipy.sparse import lil_matrix
from scipy.sparse.linalg import splu, spsolve
from .model import Beam

ReferenceBoundary = Literal["fixed-pinned", "pinned-pinned"]


@dataclass
class ReferenceStaticSolution:
    nodes: np.ndarray
    deflection: np.ndarray
    rotation: np.ndarray
    free_residual_l2: float


@dataclass
class ReferenceDynamicSolution:
    nodes: np.ndarray
    times: np.ndarray
    deflection: np.ndarray
    rotation: np.ndarray

    def midpoint_deflection(self) -> np.ndarray:
        return np.asarray([np.interp(self.nodes[-1]/2, self.nodes, w)
                           for w in self.deflection])


def reference_element_matrices(beam: Beam, length: float):
    """Return independent conventional symmetric stiffness and physical mass.

    DOFs [w_i, phi_i, w_j, phi_j]; engineering shear strain w' - phi.
    Shear term is integrated at element midpoint; bending term is exact.
    """
    if length <= 0:
        raise ValueError("Element length must be positive")
    b_shear = np.array([-1/length, -0.5, 1/length, -0.5])
    b_curv = np.array([0, -1/length, 0, 1/length])
    k = (beam.kGA * length * np.outer(b_shear, b_shear)
         + beam.EI * length * np.outer(b_curv, b_curv))
    m = np.zeros((4, 4))
    m[np.ix_([0, 2], [0, 2])] = beam.density*beam.area*length/6 * np.array([[2,1],[1,2]])
    m[np.ix_([1, 3], [1, 3])] = beam.density*beam.second_moment*length/6 * np.array([[2,1],[1,2]])
    return k, m


def assemble_reference(beam: Beam, n_nodes: int):
    if n_nodes < 3:
        raise ValueError("At least three nodes are required")
    x = np.linspace(0, beam.length, n_nodes)
    k = lil_matrix((2*n_nodes, 2*n_nodes))
    m = lil_matrix((2*n_nodes, 2*n_nodes))
    for e in range(n_nodes-1):
        ke, me = reference_element_matrices(beam, x[e+1]-x[e])
        for a in range(4):
            for b in range(4):
                k[2*e+a, 2*e+b] += ke[a,b]
                m[2*e+a, 2*e+b] += me[a,b]
    return x, k.tocsr(), m.tocsr()


def reference_free_dofs(n_nodes: int, boundary: ReferenceBoundary):
    if boundary == "fixed-pinned":
        constrained = (0, 1, 2*n_nodes-2)
    elif boundary == "pinned-pinned":
        constrained = (0, 2*n_nodes-2)
    else:
        raise ValueError("Unknown boundary")
    mask = np.ones(2*n_nodes,dtype=bool)
    mask[list(constrained)] = False
    return np.flatnonzero(mask)


def reference_load(nodes: np.ndarray, q: float, q_end: float, point: float, x_point: float):
    """Independent consistent force assembly: [vertical force, 0 moment] per node."""
    f = np.zeros(2*len(nodes))
    if not (0 <= q_end <= nodes[-1] and 0 <= x_point <= nodes[-1]):
        raise ValueError("Load positions must lie within the beam")
    for i, (xl, xr) in enumerate(zip(nodes[:-1], nodes[1:])):
        stop = min(xr,q_end)
        if stop > xl:
            span = stop-xl
            f[2*i] += q*(span-span**2/(2*(xr-xl)))
            f[2*(i+1)] += q*span**2/(2*(xr-xl))
    if point:
        i = min(np.searchsorted(nodes,x_point,side="right"),len(nodes)-1)
        if i == 0:
            f[0] += point
        else:
            fraction=(x_point-nodes[i-1])/(nodes[i]-nodes[i-1])
            f[2*(i-1)] += point*(1-fraction)
            f[2*i] += point*fraction
    return f


def solve_reference_static(n_nodes: int = 81, *, beam: Beam | None = None,
                           boundary: ReferenceBoundary="fixed-pinned",
                           distributed_intensity:float=4800.,
                           distributed_end:float|None=None,
                           point_force:float=3000.,
                           point_position:float|None=None) -> ReferenceStaticSolution:
    b=beam or Beam()
    x,K,_=assemble_reference(b,n_nodes)
    f=reference_load(x,distributed_intensity,
                     b.length/2 if distributed_end is None else distributed_end,
                     point_force,b.length/2 if point_position is None else point_position)
    free=reference_free_dofs(n_nodes,boundary)
    u=np.zeros(2*n_nodes)
    u[free]=spsolve(K[free,:][:,free],f[free])
    residual=K@u-f
    return ReferenceStaticSolution(x,u[::2],u[1::2],float(np.linalg.norm(residual[free])))


def reference_first_frequencies(beam: Beam, n_nodes: int, *, boundary: ReferenceBoundary="pinned-pinned", count:int=4) -> np.ndarray:
    x,K,M=assemble_reference(beam,n_nodes)
    free=reference_free_dofs(n_nodes,boundary)
    vals=eigh(K[free,:][:,free].toarray(),M[free,:][:,free].toarray(),
              subset_by_index=(0,count-1),eigvals_only=True)
    if np.any(vals <= 0):
        raise ValueError("Nonpositive eigenfrequency squared")
    return np.sqrt(vals)/(2*np.pi)


def solve_reference_dynamic(n_nodes:int=41, *, beam:Beam|None=None,
                            boundary:ReferenceBoundary="fixed-pinned",
                            t_end:float=.25,dt:float=.0005,
                            load_case:str="constant", duration:float=1.,
                            distributed_intensity:float=4800.,
                            distributed_end:float|None=None,
                            point_force:float=3000., point_position:float|None=None) -> ReferenceDynamicSolution:
    """Average-acceleration Newmark response of independent physical FE model."""
    # Import only the shared *load time profile*, never the mixed matrices.
    from .dynamic import load_multiplier
    b=beam or Beam()
    if not np.isfinite(dt) or not np.isfinite(t_end) or dt<=0 or t_end<=0:
        raise ValueError("Time and timestep must be finite and positive")
    n_steps=int(round(t_end/dt))
    if n_steps < 2 or not np.isclose(n_steps*dt,t_end,atol=1e-12,rtol=1e-10):
        raise ValueError("t_end/dt must be an integer >= 2")
    nodes,K,M=assemble_reference(b,n_nodes)
    free=reference_free_dofs(n_nodes,boundary)
    k=K[free,:][:,free].tocsc(); m=M[free,:][:,free].tocsc()
    xp=b.length/2 if point_position is None else point_position
    qend=b.length/2 if distributed_end is None else distributed_end
    unit=reference_load(nodes,1.,qend,0.,xp)[free]
    point=reference_load(nodes,0.,qend,point_force,xp)[free]
    times=np.linspace(0,t_end,n_steps+1)
    forces=point[None,:]+distributed_intensity*np.asarray([load_multiplier(t,load_case,duration) for t in times])[:,None]*unit[None,:]
    u=np.zeros((n_steps+1,len(free)))
    v=np.zeros_like(u); a=np.zeros_like(u)
    a[0]=splu(m).solve(forces[0])
    lu=splu(m+0.25*dt**2*k)
    for j in range(n_steps):
        up=u[j]+dt*v[j]+dt**2*.25*a[j]
        vp=v[j]+dt*.5*a[j]
        a[j+1]=lu.solve(forces[j+1]-k@up)
        u[j+1]=up+0.25*dt**2*a[j+1]
        v[j+1]=vp+dt*.5*a[j+1]
    full=np.zeros((n_steps+1,2*n_nodes))
    full[:,free]=u
    return ReferenceDynamicSolution(nodes,times,full[:,::2],full[:,1::2])
