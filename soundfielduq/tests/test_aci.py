import numpy as np
import pytest
from soundfielduq.aci import RollingACI
from soundfielduq.conformal import DynamicConformalPredictor


def test_aci_initialization_validation():
    with pytest.raises(ValueError):
        RollingACI(target_alpha=1.5)
    with pytest.raises(ValueError):
        RollingACI(lr=-0.01)


def test_aci_alpha_adaptation_direction():
    aci = RollingACI(target_alpha=0.10, lr=0.10)

    # Out of bound measurement -> alpha increases
    new_alpha = aci.update_and_predict(y_true=10.0, lower_bound=-1.0, upper_bound=1.0)
    assert new_alpha > 0.10

    # In-bound measurements -> alpha decreases back
    for _ in range(15):
        new_alpha = aci.update_and_predict(y_true=0.0, lower_bound=-1.0, upper_bound=1.0)
    assert new_alpha < 0.10


def test_dynamic_conformal_pipeline():
    aci = RollingACI(target_alpha=0.10)
    predictor = DynamicConformalPredictor(aci)

    y_true = np.linspace(1.0, 5.0, 50)
    y_pred = y_true + np.random.normal(0, 0.1, 50)
    predictor.calibrate(y_true, y_pred)

    lower, upper = predictor.predict_next(3.0)
    assert lower < 3.0 < upper

    predictor.update(y_true_t=3.0, y_pred_t=3.0, lower=lower, upper=upper)
    assert len(predictor.calibration_scores) <= aci.window_size
