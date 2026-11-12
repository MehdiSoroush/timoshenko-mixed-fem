"""Independent closed-form static validation for standard simply supported beams.

For a prismatic linearly elastic Timoshenko beam under uniform distributed
load q or a point force Q at midspan, with zero deflection at the ends and
free end rotations. These are not analytical solutions to the original
fixed-pinned / half-length load case.
"""
from .model import Beam


def simply_supported_uniform_midspan(beam: Beam, q: float) -> float:
    return 5.0 * q * beam.length**4 / (384.0 * beam.EI) + q * beam.length**2 / (8.0 * beam.kGA)


def simply_supported_center_point_midspan(beam: Beam, point_force: float) -> float:
    return point_force * beam.length**3 / (48.0 * beam.EI) + point_force * beam.length / (4.0 * beam.kGA)


def simply_supported_timoshenko_frequency(beam: Beam, mode: int = 1) -> float:
    """Exact uniform-beam Timoshenko simply-supported bending frequency (Hz).

    Independent continuum dispersion relation, including rotary inertia:
      (kGA k² - rho A omega²)(EI k² + kGA - rho I omega²) - (kGA k)² = 0
    with k = mode*pi/L. Returns the lower (bending) branch.
    """
    import numpy as np
    if mode < 1:
        raise ValueError("mode must be >= 1")
    k = mode*np.pi/beam.length
    rA=beam.density*beam.area
    rI=beam.density*beam.second_moment
    D=beam.kGA*k**2*rI+rA*(beam.EI*k**2+beam.kGA)
    C=beam.kGA*beam.EI*k**4
    disc=D**2-4*rA*rI*C
    if disc<0 and abs(disc)/D**2 < 1e-12:
        disc=0
    if disc<0:
        raise ValueError("Negative dispersion discriminant")
    omega2=2*C/(D+np.sqrt(disc))
    return float(np.sqrt(omega2)/(2*np.pi))
