# tests/test_active_sensing.py

import numpy as np
import pytest
from soundfielduq.active_sensing import ActiveSensorPlacement


def test_greedy_sensor_selection():
    """Engine must pick the candidate with highest uncertainty when unconstrained."""
    candidates = np.array([[0.0, 0.0], [1.0, 1.0], [2.0, 2.0]])
    uncertainties = np.array([0.2, 0.9, 0.4])

    placement = ActiveSensorPlacement(candidate_grid=candidates, min_distance=0.1)
    selected_idx = placement.select_next_sensor(uncertainties)

    assert selected_idx == 1  # Index 1 has maximum uncertainty (0.9)


def test_minimum_distance_constraint():
    """Engine must skip candidate with max uncertainty if within min_distance of existing sensor."""
    candidates = np.array([[0.0, 0.0], [0.05, 0.05], [2.0, 2.0]])
    uncertainties = np.array([0.8, 0.95, 0.5])

    # min_distance = 0.5m; Candidate 1 (0.05, 0.05) is too close to Candidate 0 (0.0, 0.0)
    placement = ActiveSensorPlacement(candidate_grid=candidates, min_distance=0.5)

    idx1 = placement.select_next_sensor(uncertainties)  # Picks index 1 (max uncertainty = 0.95)
    idx2 = placement.select_next_sensor(uncertainties)  # Should skip index 0 and pick index 2

    assert idx1 == 1
    assert idx2 == 2


def test_get_selected_coordinates():
    """Verify returned coordinate array shapes."""
    candidates = np.random.rand(10, 2)
    placement = ActiveSensorPlacement(candidate_grid=candidates, min_distance=0.1)

    placement.select_next_sensor(np.random.rand(10))
    placement.select_next_sensor(np.random.rand(10))

    selected_coords = placement.get_selected_coordinates()
    assert selected_coords.shape == (2, 2)
