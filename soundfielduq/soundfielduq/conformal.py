# soundfielduq/conformal.py

import numpy as np

class DynamicConformalPredictor:
    """
    Wraps dynamic prediction intervals around base sound field reconstruction models.
    """
    def __init__(self, aci_engine):
        self.aci = aci_engine
        self.calibration_scores = []

    def calibrate(self, y_true: np.ndarray, y_pred: np.ndarray):
        """Initial calibration on offline/stationary dataset."""
        self.calibration_scores = list(np.abs(y_true - y_pred))

    def predict_next(self, y_pred_t: float) -> tuple[float, float]:
        """Calculates prediction interval for time step t using current alpha."""
        q = self.aci.compute_quantile(np.array(self.calibration_scores))
        lower = y_pred_t - q
        upper = y_pred_t + q
        return lower, upper

    def update(self, y_true_t: float, y_pred_t: float, lower: float, upper: float):
        """Updates calibration scores and ACI engine upon receiving true sensor measurement."""
        score_t = abs(y_true_t - y_pred_t)
        
        # Maintain rolling window of recent non-conformity scores
        self.calibration_scores.append(score_t)
        if len(self.calibration_scores) > self.aci.window_size:
            self.calibration_scores.pop(0)
            
        # Update dynamic alpha
        self.aci.update_and_predict(y_true_t, lower, upper)
