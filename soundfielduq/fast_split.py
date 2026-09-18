"""Low-latency Split-Conformal Prediction Engine for sound field frame processing."""

import numpy as np


class FastSplitConformal:
    """Symmetric absolute-residual split-conformal prediction engine.

    Notes
    -----
    - Calibration complexity: O(N log N) due to quantile sorting, where N is calibration size.
    - Frame prediction complexity: O(M) vector addition, linear in the target frame size M.
    """

    def __init__(self, target_alpha: float = 0.10) -> None:
        if not 0.0 < target_alpha < 1.0:
            raise ValueError(f"target_alpha must be in (0, 1), got {target_alpha}")
        self.target_alpha: float = target_alpha
        self.q_threshold: float | None = None

    def calibrate(self, y_true: np.ndarray, y_pred: np.ndarray) -> None:
        """Computes non-conformity threshold over calibration set.

        Parameters
        ----------
        y_true : np.ndarray
            Reference ground truth values (N samples).
        y_pred : np.ndarray
            Model predicted values (N samples).
        """
        y_true = np.asarray(y_true, dtype=np.float64)
        y_pred = np.asarray(y_pred, dtype=np.float64)

        if y_true.shape != y_pred.shape:
            raise ValueError(f"Shape mismatch: {y_true.shape} vs {y_pred.shape}")

        scores = np.abs(y_true - y_pred).ravel()
        n = scores.size
        if n == 0:
            raise ValueError("Calibration dataset cannot be empty.")

        q_level = np.ceil((n + 1) * (1.0 - self.target_alpha)) / n
        q_level = float(np.clip(q_level, 0.0, 1.0))
        self.q_threshold = float(np.quantile(scores, q_level, method="higher"))

    def predict_frame(
        self, y_pred_frame: np.ndarray
    ) -> tuple[np.ndarray, np.ndarray]:
        """Calculates point-wise conformal prediction intervals for an audio frame.

        Parameters
        ----------
        y_pred_frame : np.ndarray
            Predicted frame values of size M.

        Returns
        -------
        tuple[np.ndarray, np.ndarray]
            (lower_bounds, upper_bounds) both of shape matching input frame.
        """
        if self.q_threshold is None:
            raise RuntimeError("Engine is uncalibrated. Run calibrate() first.")

        y_pred_frame = np.asarray(y_pred_frame, dtype=np.float64)
        return y_pred_frame - self.q_threshold, y_pred_frame + self.q_threshold
