"""Beam stress recovery with transparent physical and historical conventions.

NEW verification work (2026), not the original MATLAB source. The source
uses h/(4I) times the midpoint moment for sigma_11, while the full outer-
fibre bending stress has magnitude h/(2I) times the midpoint moment.

Two distinct shear conventions are exposed: section-average V/A, peak
rectangular-section shear 3V/(2A), and the source's effective G*gamma=V/(ks A).
"""
from __future__ import annotations
from dataclasses import dataclass
import numpy as np
from .model import Beam


@dataclass
class ElementStress:
    x_mid: np.ndarray
    moment_mid: np.ndarray
    shear_force: np.ndarray
    sigma_outer: np.ndarray
    sigma_source_2017: np.ndarray
    tau_section_average: np.ndarray
    tau_rectangular_max: np.ndarray
    tau_effective_timoshenko: np.ndarray


def recover_element_stress(nodes:np.ndarray,moments:np.ndarray,beam:Beam|None=None)->ElementStress:
    """Recover piecewise bending stress and shear from nodal mixed moments.

    `moments` has last axis = number of nodes; any preceding dimensions
    (including time) are preserved.  Sign follows the original moment DOF.
    """
    b=beam or Beam()
    nodes=np.asarray(nodes,dtype=float)
    moments=np.asarray(moments,dtype=float)
    if nodes.ndim!=1 or len(nodes)<2 or np.any(np.diff(nodes)<=0):
        raise ValueError("nodes must be strictly increasing 1D coordinates")
    if moments.shape[-1]!=len(nodes) or not np.isfinite(moments).all():
        raise ValueError("moment array must have a finite value at each node")
    mm=.5*(moments[...,1:]+moments[...,:-1])
    shear=np.diff(moments,axis=-1)/np.diff(nodes)
    return ElementStress(
        x_mid=.5*(nodes[:-1]+nodes[1:]),
        moment_mid=mm,
        shear_force=shear,
        sigma_outer=mm*b.height/(2*b.second_moment),
        sigma_source_2017=mm*b.height/(4*b.second_moment),
        tau_section_average=shear/b.area,
        tau_rectangular_max=1.5*shear/b.area,
        tau_effective_timoshenko=shear/(b.shear_correction*b.area),
    )
