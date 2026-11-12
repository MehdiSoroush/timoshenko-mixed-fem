"""Independent 2026 displacement/rotation, continuum, and stress checks."""
import numpy as np
import pytest
from timoshenko_mixed.model import Beam, solve_static, midpoint_deflection
from timoshenko_mixed.dynamic import solve_dynamic
from timoshenko_mixed.analytical import (simply_supported_uniform_midspan,
    simply_supported_timoshenko_frequency)
from timoshenko_mixed.reference import (
    reference_element_matrices, solve_reference_static,
    reference_first_frequencies, solve_reference_dynamic)
from timoshenko_mixed.stress import recover_element_stress


def test_independent_element_is_symmetric_with_positive_nonrigid_energy():
    b=Beam();K,M=reference_element_matrices(b,.1)
    np.testing.assert_allclose(K,K.T,atol=1e-7)
    np.testing.assert_allclose(M,M.T,atol=1e-12)
    assert np.linalg.eigvalsh(M).min()>0
    # Two independent rigid-body modes (vertical translation and rigid rotation)
    np.testing.assert_allclose(K@np.array([1.,0.,1.,0.]),0.,atol=1e-7)
    np.testing.assert_allclose(K@np.array([0.,1.,.1,1.]),0.,atol=1e-7)


def test_reference_uniform_approaches_exact_static_solution():
    b=Beam(); exact=simply_supported_uniform_midspan(b,4800.)
    errors=[]
    for n in (17,33,65):
        r=solve_reference_static(n,beam=b,boundary='pinned-pinned',
               distributed_intensity=4800,distributed_end=2,point_force=0)
        errors.append(abs(r.deflection[n//2]-exact))
    assert errors[2] < errors[1]/3.9
    assert errors[1] < errors[0]/3.9
    assert errors[-1]/exact < .001


def test_independent_static_agrees_with_mixed_fixed_pinned():
    n=129;b=Beam()
    m=solve_static(n,b)
    r=solve_reference_static(n,beam=b)
    # Distinct formulations converge to the same physical solution.
    assert abs(midpoint_deflection(m)-r.deflection[n//2])/abs(r.deflection[n//2])<.0003


def test_independent_timoshenko_continuum_frequency():
    b=Beam();exact=simply_supported_timoshenko_frequency(b)
    computed=reference_first_frequencies(b,41,count=1)[0]
    assert abs(computed/exact - 1)<.001
    assert exact == pytest.approx(24.721389716118285,rel=1e-8)
    with pytest.raises(ValueError):
        simply_supported_timoshenko_frequency(b,mode=0)


def test_stress_outer_fiber_is_twice_original_report_formula():
    b=Beam()
    nodes=np.array([0.,1.,2.]); moments=np.array([0.,100.,0.])
    e=recover_element_stress(nodes,moments,b)
    np.testing.assert_allclose(e.sigma_outer,2*e.sigma_source_2017)
    np.testing.assert_allclose(e.tau_section_average,e.shear_force/b.area)
    np.testing.assert_allclose(e.tau_rectangular_max,1.5*e.tau_section_average)
    np.testing.assert_allclose(e.tau_effective_timoshenko,e.tau_section_average/b.shear_correction)
    with pytest.raises(ValueError):
        recover_element_stress(np.array([0.,0.]),moments[:2],b)


def test_stress_accepts_dynamic_time_history():
    s=solve_dynamic(9,dt=.005,t_end=.02)
    stress=recover_element_stress(s.nodes,s.bending_moment)
    assert stress.sigma_outer.shape == (5,8)
    assert np.isfinite(stress.tau_rectangular_max).all()


def test_independent_dynamic_matches_low_frequency_fixed_pinned_response():
    # Both solvers have distinct physical DOFs/mass matrices. Agreement is
    # assessed for the macroscopic midspan response, not algebraic residual.
    common=dict(n_nodes=41,boundary='fixed-pinned',dt=.001,t_end=.15,
                load_case='constant')
    mixed=solve_dynamic(**common,integrator='newmark-average')
    ref=solve_reference_dynamic(**common)
    wm=mixed.midpoint_deflection();wr=ref.midpoint_deflection()
    rmse=float(np.sqrt(np.mean((wm-wr)**2)))
    assert rmse/np.max(abs(wr))<.04


def test_source_moment_sign_and_independent_shear_force():
    """The source M DOF is opposite to EI*d(phi)/dx in the reference convention."""
    b=Beam();n=129
    m=solve_static(n,beam=b)
    r=solve_reference_static(n,beam=b)
    le=np.diff(r.nodes)
    physical_moment=b.EI*np.diff(r.rotation)/le
    source_element_moment=.5*(m.bending_moment[1:]+m.bending_moment[:-1])
    physical_shear=b.kGA*(np.diff(r.deflection)/le-.5*(r.rotation[1:]+r.rotation[:-1]))
    source_shear=np.diff(m.bending_moment)/le
    assert np.linalg.norm(physical_moment+source_element_moment)/np.linalg.norm(physical_moment)<.0001
    assert np.linalg.norm(physical_shear-source_shear)/np.linalg.norm(physical_shear)<.0001
