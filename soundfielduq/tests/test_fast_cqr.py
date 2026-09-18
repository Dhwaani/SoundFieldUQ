import numpy as np
import pytest
from soundfielduq.fast_cqr import FastSplitConformal


def test_fast_cqr_uncalibrated_raises():
    engine = FastSplitConformal()
    with pytest.raises(RuntimeError):
        engine.predict_frame_fast(np.array([1.0, 2.0]))


def test_fast_cqr_vectorized_bounds():
    y_true = np.random.normal(0, 1, (100, 10))
    y_pred = y_true + np.random.normal(0, 0.05, (100, 10))

    engine = FastSplitConformal(target_alpha=0.10)
    engine.calibrate(y_true, y_pred)

    frame = np.ones(10)
    lower, upper = engine.predict_frame_fast(frame)

    assert lower.shape == (10,)
    assert np.all(lower < frame)
    assert np.all(upper > frame)
