"""Adaptive Conformal Inference (ACI) for dynamic sound fields."""

from typing import Union
import numpy as np


class RollingACI:
    """Adaptive Conformal Inference engine using online learning gradient updates.

    Dynamically adjusts target miscoverage rate (alpha) in response to empirical
    coverage errors during non-stationary acoustic domain shifts.
    """

    def __init__(
        self,
        target_alpha: float = 0.10,
        lr: float = 0.05,
        window_size: int = 100,
        alpha_bounds: tuple[float, float] = (0.001, 0.999),
    ) -> None:
        if not 0.0 < target_alpha < 1.0:
            raise ValueError(f"target_alpha must be in (0, 1), got {target_alpha}")
        if lr <= 0.0:
            raise ValueError(f"Learning rate lr must be positive, got {lr}")
        if window_size < 1:
            raise ValueError(f"window_size must be >= 1, got {window_size}")

        self.target_alpha: float = target_alpha
        self.lr: float = lr
        self.window_size: int = window_size
        self.alpha_bounds: tuple[float, float] = alpha_bounds
        self.alpha_t: float = target_alpha

    def update_and_predict(
        self, y_true: float, lower_bound: float, upper_bound: float
    ) -> float:
        """Evaluates empirical coverage of step t and updates target alpha.

        Parameters
        ----------
        y_true : float
            Ground truth acoustic pressure measurement.
        lower_bound : float
            Lower bound of predicted conformal interval.
        upper_bound : float
            Upper bound of predicted conformal interval.

        Returns
        -------
        float
            Updated alpha parameter for subsequent frame interval estimation.
        """
        covered = lower_bound <= y_true <= upper_bound
        err = 0.0 if covered else 1.0

        # Gradient update step (Angelopoulos et al., 2021)
        self.alpha_t += self.lr * (self.target_alpha - err)
        self.alpha_t = float(
            np.clip(self.alpha_t, self.alpha_bounds[0], self.alpha_bounds[1])
        )
        return self.alpha_t

    def compute_quantile(self, scores: np.ndarray) -> float:
        """Computes score quantile corresponding to current adaptive alpha_t.

        Parameters
        ----------
        scores : np.ndarray
            1D array of non-conformity calibration residual scores.

        Returns
        -------
        float
            Calibrated quantile threshold value.
        """
        scores = np.asarray(scores, dtype=np.float64).ravel()
        n = scores.size
        if n == 0:
            return 0.0

        q_level = np.ceil((n + 1) * (1.0 - self.alpha_t)) / n
        q_level = float(np.clip(q_level, 0.0, 1.0))
        return float(np.quantile(scores, q_level, method="higher"))
