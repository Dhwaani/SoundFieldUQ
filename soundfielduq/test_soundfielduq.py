import numpy as np
import pytest
from soundfielduq.aci import RollingACI
from soundfielduq.active_sensing import ActiveSensorPlacement
from soundfielduq.conformal import DynamicConformalPredictor
from soundfielduq.fast_split import FastSplitConformal


@pytest.fixture
def rng():
    return np.random.default_rng(seed=42)


def test_aci_exact_quantile():
    aci = RollingACI(target_alpha=0.10)
    scores = np.array([0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0])
    q = aci.compute_quantile(scores)
    assert np.isclose(q, 1.0)  # (10+1)*(0.9)/10 = 0.99 -> 1.0 quantile


def test_aci_invalid_inputs():
    with pytest.raises(ValueError):
        RollingACI(target_alpha=1.5)
    with pytest.raises(ValueError):
        RollingACI(lr=-0.1)


def test_fast_split_calibration_and_coverage(rng):
    engine = FastSplitConformal(target_alpha=0.10)
    y_true = rng.normal(0, 1, 1000)
    y_pred = y_true + rng.normal(0, 0.1, 1000)

    engine.calibrate(y_true, y_pred)
    assert engine.q_threshold > 0.0

    test_frame = rng.normal(0, 1, 100)
    lower, upper = engine.predict_frame(test_frame)
    assert lower.shape == (100,)
    assert np.all(lower < upper)


def test_active_sensing_dynamic_update():
    grid = np.array([[0.0, 0.0], [0.1, 0.1], [1.0, 1.0], [2.0, 2.0]])
    placement = ActiveSensorPlacement(candidate_grid=grid, min_distance=0.5)

    unc = np.array([0.9, 0.95, 0.4, 0.3])
    idx1 = placement.select_next_sensor(unc)
    assert idx1 == 1  # Pick 0.95 first

    # Dynamic step: Next step masks out index 0 (within 0.5m of index 1)
    idx2 = placement.select_next_sensor(unc)
    assert idx2 == 2  # Pick index 2 (1.0, 1.0)


def test_missing_candidates_raises_error():
    grid = np.array([[0.0, 0.0], [0.1, 0.1]])
    placement = ActiveSensorPlacement(candidate_grid=grid, min_distance=0.5)
    placement.select_next_sensor(np.array([0.5, 0.8]))

    with pytest.raises(RuntimeError):
        placement.select_next_sensor(np.array([0.5, 0.8]))
