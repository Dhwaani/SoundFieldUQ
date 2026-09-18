# tests/test_aci.py

import numpy as np
import pytest
from soundfielduq.aci import RollingACI

def test_aci_alpha_increases_on_miscoverage():
    """Alpha should increase (tighten intervals) when miscoverage occurs repeatedly."""
    aci = RollingACI(target_alpha=0.10, lr=0.05)
    initial_alpha = aci.alpha_t
    
    # Simulate repeated miscoverage (true value falls outside bounds)
    for _ in range(5):
        aci.update_and_predict(y_true=10.0, lower_bound=-1.0, upper_bound=1.0)
        
    assert aci.alpha_t > initial_alpha

def test_aci_alpha_decreases_on_coverage():
    """Alpha should decrease (widen intervals) when bounds consistently cover."""
    aci = RollingACI(target_alpha=0.10, lr=0.05)
    initial_alpha = aci.alpha_t
    
    # Simulate continuous valid coverage
    for _ in range(5):
        aci.update_and_predict(y_true=0.0, lower_bound=-1.0, upper_bound=1.0)
        
    assert aci.alpha_t < initial_alpha

def test_alpha_clipping_bounds():
    """Alpha should remain strictly within valid probability bounds [0.001, 0.999]."""
    aci = RollingACI(target_alpha=0.10, lr=0.5)
    
    # Force extreme miscoverage updates
    for _ in range(20):
        aci.update_and_predict(y_true=100.0, lower_bound=0.0, upper_bound=1.0)
    assert aci.alpha_t <= 0.999
    
    # Force extreme valid coverage updates
    for _ in range(20):
        aci.update_and_predict(y_true=0.5, lower_bound=0.0, upper_bound=1.0)
    assert aci.alpha_t >= 0.001
