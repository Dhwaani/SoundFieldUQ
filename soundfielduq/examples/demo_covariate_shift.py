# examples/demo_covariate_shift.py

import numpy as np
from soundfielduq.aci import RollingACI
from soundfielduq.conformal import DynamicConformalPredictor

def run_covariate_shift_experiment():
    np.random.seed(42)
    n_steps = 300
    target_alpha = 0.10  # 90% coverage
    
    # Generate Synthetic Sound Field Prediction Errors
    # Phase 1 (t=0..150): Low noise room acoustics
    # Phase 2 (t=150..300): Covariate Shift (Increased room reverberation / sensor noise)
    true_signal = np.sin(np.linspace(0, 10 * np.pi, n_steps))
    base_predictions = true_signal + np.random.normal(0, 0.1, n_steps)
    
    # Introduce sudden acoustic domain shift at t = 150
    shifted_errors = np.random.normal(0, 0.35, n_steps - 150)
    base_predictions[150:] = true_signal[150:] + shifted_errors

    # Setup Predictors
    aci_engine = RollingACI(target_alpha=target_alpha, lr=0.05, window_size=50)
    predictor = DynamicConformalPredictor(aci_engine)
    
    # Pre-calibrate on initial stationary data
    predictor.calibrate(true_signal[:30], base_predictions[:30])

    # Online Tracking
    aci_coverage = []
    
    for t in range(30, n_steps):
        y_true = true_signal[t]
        y_pred = base_predictions[t]
        
        # Predict & Evaluate
        lower, upper = predictor.predict_next(y_pred)
        covered = lower <= y_true <= upper
        aci_coverage.append(covered)
        
        # Online Update
        predictor.update(y_true, y_pred, lower, upper)

    post_shift_coverage = np.mean(aci_coverage[120:]) # Coverage after shift
    print(f"Target Coverage: {1 - target_alpha:.2%}")
    print(f"Post-Shift Empirical Coverage (ACI): {post_shift_coverage:.2%}")

if __name__ == "__main__":
    run_covariate_shift_experiment()
