# soundfielduq/active_sensing.py

import numpy as np


class ActiveSensorPlacement:
    """
    Greedy sensor placement engine using Uncertainty Quantification (UQ).
    Selects optimal sensor locations to minimize sound field reconstruction uncertainty.
    """

    def __init__(self, candidate_grid: np.ndarray, min_distance: float = 0.2):
        """
        Args:
            candidate_grid: (N, 2) or (N, 3) spatial coordinates of possible microphone locations.
            min_distance: Minimum allowable physical distance between sensors (in meters).
        """
        self.candidate_grid = candidate_grid
        self.min_distance = min_distance
        self.selected_indices = []

    def select_next_sensor(self, spatial_uncertainty: np.ndarray) -> int:
        """
        Selects the candidate position with maximum uncertainty, respecting physical constraints.

        Args:
            spatial_uncertainty: (N,) array of prediction interval widths at candidate positions.

        Returns:
            index of the selected candidate sensor position.
        """
        # Sort candidate indices by descending uncertainty
        ranked_indices = np.argsort(spatial_uncertainty)[::-1]

        for idx in ranked_indices:
            candidate_pos = self.candidate_grid[idx]

            # Check distance against already selected sensors
            if len(self.selected_indices) > 0:
                selected_positions = self.candidate_grid[self.selected_indices]
                distances = np.linalg.norm(selected_positions - candidate_pos, axis=1)

                if np.any(distances < self.min_distance):
                    continue  # Skip if too close to an existing sensor

            self.selected_indices.append(idx)
            return idx

        raise RuntimeError("No candidate locations available matching distance constraints.")

    def get_selected_coordinates(self) -> np.ndarray:
        """Returns physical coordinates of all selected sensors."""
        return self.candidate_grid[self.selected_indices]
