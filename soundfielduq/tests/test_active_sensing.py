import numpy as np
import pytest
from soundfielduq.active_sensing import ActiveSensorPlacement


def test_active_sensing_grid_validation():
    with pytest.raises(ValueError):
        ActiveSensorPlacement(candidate_grid=np.ones((10, 4)))


def test_active_sensing_distance_constraints():
    candidates = np.array([[0.0, 0.0], [0.05, 0.05], [1.0, 1.0]])
    placement = ActiveSensorPlacement(candidate_grid=candidates, min_distance=0.5)

    uncertainties = np.array([0.7, 0.95, 0.4])
    idx1 = placement.select_next_sensor(uncertainties)
    idx2 = placement.select_next_sensor(uncertainties)

    assert idx1 == 1  # Maximum uncertainty picked first
    assert idx2 == 2  # Skip candidate 0 because it's < 0.5m from candidate 1
