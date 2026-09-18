"""Sequential Active Sensor Placement with Dynamic Variance/Uncertainty Update."""

import numpy as np


class ActiveSensorPlacement:
    """Sequential greedy microphone placement engine with dynamic score updating."""

    def __init__(
        self, candidate_grid: np.ndarray, min_distance: float = 0.20
    ) -> None:
        grid = np.asarray(candidate_grid, dtype=np.float64)
        if grid.ndim != 2 or grid.shape[1] not in (2, 3):
            raise ValueError(f"candidate_grid must be (N, 2) or (N, 3), got {grid.shape}")
        if min_distance < 0.0:
            raise ValueError(f"min_distance must be non-negative, got {min_distance}")

        self.candidate_grid: np.ndarray = grid
        self.min_distance: float = min_distance
        self.selected_indices: list[int] = []

    def select_next_sensor(
        self, spatial_covariance_or_unc: np.ndarray, update_fn=None
    ) -> int:
        """Selects next optimal sensor and updates spatial uncertainty state.

        Parameters
        ----------
        spatial_covariance_or_unc : np.ndarray
            Current spatial uncertainty vector (N,) or covariance matrix (N, N).
        update_fn : callable, optional
            Function `fn(current_unc, selected_idx) -> updated_unc` to recompute
            uncertainty field after selection.

        Returns
        -------
        int
            Index of chosen sensor candidate.
        """
        unc = np.asarray(spatial_covariance_or_unc, dtype=np.float64)
        if unc.ndim == 2:
            unc = np.diag(unc)  # Extract marginal variances if covariance given

        unc = unc.ravel()
        if unc.size != self.candidate_grid.shape[0]:
            raise ValueError("Uncertainty dimension must match candidate count.")

        # Zero out candidates violating spatial proximity constraints
        masked_unc = unc.copy()
        for selected in self.selected_indices:
            distances = np.linalg.norm(
                self.candidate_grid - self.candidate_grid[selected], axis=1
            )
            masked_unc[distances < self.min_distance] = -np.inf

        best_idx = int(np.argmax(masked_unc))
        if masked_unc[best_idx] == -np.inf:
            raise RuntimeError("No candidate locations available matching distance constraints.")

        self.selected_indices.append(best_idx)
        return best_idx

    def get_selected_coordinates(self) -> np.ndarray:
        """Returns physical coordinates of selected sensors."""
        if not self.selected_indices:
            return np.empty((0, self.candidate_grid.shape[1]), dtype=np.float64)
        return self.candidate_grid[self.selected_indices]
