"""Dynamic conformal prediction pipeline wrapper."""

from typing import Union
import numpy as np
from soundfielduq.aci import RollingACI


class DynamicConformalPredictor:
    """Wrapper pipeline coordinating online ACI state and score window buffers."""

    def __init__(self, aci_engine: RollingACI) -> None:
        self.aci: RollingACI = aci_engine
        self.calibration_scores: list[float] = []

    def calibrate(self, y_true: np.ndarray, y_pred: np.ndarray) -> None:
        """Populates initial non-conformity calibration memory buffer.

        Parameters
        ----------
        y_true : np.ndarray
            Reference ground truth signal vector or matrix.
        y_pred : np.ndarray
            Predicted baseline pressure values matching y_true dimensions.
        """
        y_true = np.asarray(y_true, dtype=np.float64)
        y_pred = np.asarray(y_pred, dtype=np.float64)

        if y_true.shape != y_pred.shape:
            raise ValueError(
                f"Shape mismatch: y_true {y_true.shape} vs y_pred {y_pred.shape}"
            )

        scores = np.abs(y_true - y_pred).ravel()
        self.calibration_scores = list(scores[-self.aci.window_size :])

    def predict_next(self, y_pred_t: float) -> tuple[float, float]:
        """Calculates prediction interval bounds for step t."""
        if not self.calibration_scores:
            raise RuntimeError("Predictor is uncalibrated. Run calibrate() first.")

        q = self.aci.compute_quantile(np.array(self.calibration_scores))
        return float(y_pred_t - q), float(y_pred_t + q)

    def update(
        self, y_true_t: float, y_pred_t: float, lower: float, upper: float
    ) -> None:
        """Ingests true ground truth observation to update memory buffer and ACI state."""
        score_t = float(abs(y_true_t - y_pred_t))
        self.calibration_scores.append(score_t)

        if len(self.calibration_scores) > self.aci.window_size:
            self.calibration_scores.pop(0)

        self.aci.update_and_predict(y_true_t, lower, upper)
