# examples/demo_active_sensing.py

import numpy as np
from soundfielduq.active_sensing import ActiveSensorPlacement


def simulate_sound_field(grid_x: np.ndarray, grid_y: np.ndarray, source_pos: tuple = (1.5, 2.0), k: float = 5.0):
    """Simulates Helmholtz pressure wave field from a point source."""
    dx = grid_x - source_pos[0]
    dy = grid_y - source_pos[1]
    dist = np.sqrt(dx**2 + dy**2) + 1e-6
    pressure = np.cos(k * dist) / dist
    return pressure


def run_active_sensing_demo():
    # Define a 2D room spatial grid (5m x 5m, resolution 0.1m)
    x = np.linspace(0, 5, 50)
    y = np.linspace(0, 5, 50)
    xx, yy = np.meshgrid(x, y)
    grid_coords = np.column_stack([xx.ravel(), yy.ravel()])

    # True synthetic sound field
    true_field = simulate_sound_field(xx, yy).ravel()

    # Initialize Active Sensing Engine
    placement = ActiveSensorPlacement(candidate_grid=grid_coords, min_distance=0.5)

    num_sensors = 8
    print(f"Placing {num_sensors} microphones actively based on UQ...")

    for step in range(num_sensors):
        if step == 0:
            # Place initial sensor arbitrarily or at domain center
            uncertainty_map = np.ones(len(grid_coords))
        else:
            # Synthetic spatial uncertainty proxy: distance to nearest chosen sensor * field energy
            chosen_coords = placement.get_selected_coordinates()
            distances = np.min(
                np.linalg.norm(grid_coords[:, None, :] - chosen_coords[None, :, :], axis=2),
                axis=1,
            )
            uncertainty_map = distances * (np.abs(true_field) + 0.1)

        idx = placement.select_next_sensor(uncertainty_map)
        pos = grid_coords[idx]
        print(f" Sensor {step + 1}: Positioned at ({pos[0]:.2f}m, {pos[1]:.2f}m) | Uncertainty: {uncertainty_map[idx]:.4f}")


if __name__ == "__main__":
    run_active_sensing_demo()
