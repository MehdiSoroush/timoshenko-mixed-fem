"""Regression and temporal-verification tests for the dynamic reconstruction.

Verification is limited to the discretized matrix equations plus an
independent single-mode ODE reference. These are not Abaqus comparisons.
"""
import numpy as np
import pytest
from scipy.linalg import eig

from timoshenko_mixed import Beam, assemble_mass, element_mass, load_multiplier, solve_dynamic
from timoshenko_mixed.model import assemble_stiffness
from timoshenko_mixed.dynamic import _free_dofs


def test_source_element_mass_coefficients():
    b = Beam()
    le = .4
    c, d, g = b.density*b.area, b.density/b.young_modulus, b.density/(b.kGA*b.young_modulus)
    expected = le*np.array([
        [c/3, d/3, c/6, d/6], [0, g/3, 0, g/6],
        [c/6, d/6, c/3, d/3], [0, g/6, 0, g/3],
    ])
    np.testing.assert_allclose(element_mass(b, le), expected)
    # Nonsymmetry is intentional (the source matrix is a mixed formulation).
    assert not np.allclose(expected / np.abs(expected).max(), expected.T / np.abs(expected).max(), atol=1e-12)


def test_mass_matrix_translational_mass_check():
    b = Beam()
    x, M = assemble_mass(b, 11)
    w = M[::2, ::2].toarray()
    assert M.shape == (22,22)
    assert np.isclose(np.sum(w), b.density*b.area*b.length)


@pytest.mark.parametrize('case,expected', [
    ('constant', 1.0), ('sine', np.sin(.5)),
    ('source-ramp', .5), ('unit-triangle', 1.0),
])
def test_midpulse_load_profiles(case,expected):
    assert load_multiplier(.5, case) == pytest.approx(expected)
    assert load_multiplier(1.01, case) == 0.0


def test_source_ramp_is_not_a_unit_triangle():
    assert load_multiplier(.25,'source-ramp') == pytest.approx(.25)
    assert load_multiplier(.75,'source-ramp') == pytest.approx(.25)
    assert load_multiplier(.75,'unit-triangle') == pytest.approx(.5)


@pytest.mark.parametrize('method', ['source-three-level','newmark-average'])
def test_all_three_source_loads_finite_and_enforce_supports(method):
    for shape in ('constant','sine','source-ramp'):
        sol=solve_dynamic(9, t_end=.10, dt=.005, load_case=shape, integrator=method)
        assert np.isfinite(sol._all_dofs).all()
        assert np.max(abs(sol.deflection[:,0])) == 0
        assert np.max(abs(sol.deflection[:,-1])) == 0
        assert np.max(abs(sol.bending_moment[:,-1])) == 0
        assert sol.startup_residual < 1e-5
        assert sol.maximum_step_residual < 1e-5


def test_unforced_zero_state_remains_zero():
    for method in ('source-three-level','newmark-average'):
        sol=solve_dynamic(9, dt=.005, t_end=.1, point_force=0,
                          distributed_intensity=0, integrator=method)
        assert np.all(sol._all_dofs == 0)


def _first_mode(n_nodes=9):
    """Independent generalized eigenproblem for a semidiscrete modal reference."""
    b=Beam(); x,K=assemble_stiffness(b,n_nodes); _,M=assemble_mass(b,n_nodes)
    free=_free_dofs(n_nodes,'pinned-pinned')
    # This reduced matrix pencil contains very high-frequency modes. The
    # smallest physically relevant positive eigenvalue is well separated.
    v, vectors=eig(K.toarray()[np.ix_(free,free)],M.toarray()[np.ix_(free,free)])
    finite=np.flatnonzero(np.isfinite(v) & (v.real > 0) & (np.abs(v.imag) < 1e-5*v.real))
    idx=finite[np.argmin(v.real[finite])]
    omega=float(np.sqrt(v[idx].real))
    initial=np.zeros(2*n_nodes)
    initial[free]=vectors[:,idx].real
    initial*=0.001/max(abs(initial[::2]))
    return b, initial, omega


def test_first_natural_frequency_converges_with_mesh():
    w=[]
    for n in (7, 11, 17):
        _, _, omega=_first_mode(n)
        w.append(omega)
    assert 20 < w[-1]/(2*np.pi) < 30  # expected simply-supported order of magnitude
    assert abs(w[2]-w[1]) < abs(w[1]-w[0])


def test_newmark_modal_response_is_second_order_accurate():
    _, x0, omega=_first_mode(9)
    errors=[]
    for dt in (.002, .001):
        sol=solve_dynamic(9, boundary='pinned-pinned', dt=dt, t_end=.04,
                          distributed_intensity=0, point_force=0,
                          initial_dofs=x0, integrator='newmark-average')
        exact=x0[2*(9//2)]*np.cos(omega*sol.times)
        errors.append(np.max(np.abs(sol.midpoint_deflection()-exact)) / abs(exact[0]))
    assert errors[1] < errors[0]/3.5
    assert errors[1] < .015


def test_source_three_level_damps_free_modal_oscillation():
    _, x0, _=_first_mode(9)
    s=solve_dynamic(9, boundary='pinned-pinned', dt=.001, t_end=.04,
                    distributed_intensity=0, point_force=0, initial_dofs=x0)
    n=solve_dynamic(9, boundary='pinned-pinned', dt=.001, t_end=.04,
                    distributed_intensity=0, point_force=0, initial_dofs=x0,
                    integrator='newmark-average')
    # Characterize numerical damping of the weighted-residual recurrence.
    assert abs(s.midpoint_deflection()[-1]) < .75*abs(n.midpoint_deflection()[-1])


def test_invalid_time_grid_and_initial_dofs_rejected():
    with pytest.raises(ValueError,match='integer'):
        solve_dynamic(9,t_end=.1,dt=.03)
    with pytest.raises(ValueError,match='length'):
        solve_dynamic(9,t_end=.1,dt=.01,initial_dofs=np.ones(8))
    with pytest.raises(ValueError,match='constrained'):
        solve_dynamic(9,t_end=.1,dt=.01,initial_dofs=np.ones(18))
    with pytest.raises(ValueError,match='gamma'):
        solve_dynamic(9,t_end=.1,dt=.01,gamma=0.5)
    with pytest.raises(ValueError,match='Unsupported'):
        load_multiplier(.5,'invalid')
