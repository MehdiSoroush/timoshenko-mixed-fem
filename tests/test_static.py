"""Tests for source-derived static mixed FEM and independent analytical checks."""
import numpy as np
import pytest
from timoshenko_mixed import (
    Beam, midpoint_deflection, solve_static,
    simply_supported_uniform_midspan, simply_supported_center_point_midspan,
)
from timoshenko_mixed.model import element_stiffness, assemble_stiffness


def test_corrected_local_stiffness_matches_dynamic_listing_structure():
    b = Beam(); le = 0.5
    a = le**2 / (3*b.EI) + 1/b.kGA
    beta = le**2 / (6*b.EI) - 1/b.kGA
    expected = np.array([[0,1,0,-1],[-1,a,1,beta],[0,-1,0,1],[1,beta,-1,a]])/le
    np.testing.assert_allclose(element_stiffness(b, le), expected)


def test_assembled_matrix_has_expected_shape():
    x, K = assemble_stiffness(Beam(), 5)
    assert x.size == 5
    assert K.shape == (10, 10)


def test_source_case_boundary_values_and_small_free_residual():
    sol = solve_static(55)
    assert abs(sol.deflection[0]) < 1e-13
    assert abs(sol.deflection[-1]) < 1e-13
    assert abs(sol.bending_moment[-1]) < 1e-10
    assert sol.free_residual_l2 < 1e-6


def test_original_case_midspan_converges():
    values = [midpoint_deflection(solve_static(n)) for n in (9, 17, 33, 65)]
    assert all(x > 0 for x in values)
    differences = np.abs(np.diff(values))
    assert differences[-1] < differences[0]


def test_uniform_load_matches_independent_timoshenko_formula():
    b = Beam()
    exact = simply_supported_uniform_midspan(b, 4800.0)
    calc = midpoint_deflection(solve_static(
        65, b, boundary="pinned-pinned", distributed_intensity=4800,
        distributed_end=b.length, point_force=0))
    assert abs(calc / exact - 1) < 0.001


def test_midpoint_force_matches_independent_timoshenko_formula():
    b = Beam()
    exact = simply_supported_center_point_midspan(b, 3000.0)
    calc = midpoint_deflection(solve_static(
        65, b, boundary="pinned-pinned", distributed_intensity=0,
        point_force=3000))
    assert abs(calc / exact - 1) < 0.001


def test_even_mesh_interpolates_midpoint_load():
    sol = solve_static(10)
    assert np.isfinite(sol.deflection).all()
    assert sol.deflection.max() > 0


def test_invalid_inputs_are_rejected():
    with pytest.raises(ValueError):
        solve_static(2)
    with pytest.raises(ValueError):
        Beam(young_modulus=-1.0)
    with pytest.raises(ValueError):
        solve_static(5, distributed_end=5.0)
