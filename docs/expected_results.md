# Expected results

These benchmark values serve as reference checkpoints to verify your implementation. They were generated using an independent Python implementation before writing the MATLAB code to provide a cross-language validation check.

**Note on random seeds**: MATLAB and NumPy use different pseudo-random number generators, so outputs won't match to the last decimal place. Agreement within roughly $\pm 0.02$ on coverage and $\pm 1\text{ dB}$ on NMSE indicates a successful port.

## Sanity checks (should be near-exact)

| Quantity | Expected |
|---|---|
| Helmholtz kernel min eigenvalue, k = 3.7 … 37 | ≥ −1e−9 (PSD) |
| Helmholtz residual of synthetic field, 300 Hz | ~3e−5 relative |
| Helmholtz residual of synthetic field, 1000 Hz | ~3e−4 relative |
| KRR mean vs GP mean, λ = σ² | identical to 1e−8 |

## Experiment 01 — accuracy (800 Hz, 1 m aperture)

NMSE drops rapidly once the aperture is sampled below half a wavelength. The Helmholtz kernel should consistently outperform the squared-exponential (RBF) baseline across all microphone counts.

| Microphones | Helmholtz (dB) | Squared-exponential (dB) |
|---|---|---|
| 40 | ≈ −2.8 | ≈ +1.7 |
| 80 | ≈ −7.8 | ≈ −0.5 |
| 160 | ≈ −19.3 | ≈ −6.9 |

If the RBF control *wins*, the length scale is mismatched — check `ell = pi/k`.

## Experiment 02 — marginal calibration (120 mics, 90% target)

All three methods land close to nominal. **This is the expected, honest result.** A fitted Helmholtz GP is already well calibrated for interpolation.

| Method | Coverage at 90% nominal |
|---|---|
| GP posterior | ≈ 0.89 |
| Split conformal | ≈ 0.90 |
| Normalized conformal | ≈ 0.90 |

*Troubleshooting*: If the GP reports ≈ 1.000 with enormous radii, the signal variance is not being fitted — check `fitGP`'s closed-form `s`.

## Experiment 03 — conditional coverage (the headline)

Marginal coverage near 0.90 for all three, but the spatial breakdown diverges sharply. Coverage by distance from the aperture centre:

| Distance band | GP | Split conformal | Normalized |
|---|---|---|---|
| 0.00 – 0.35 m | ≈ 0.89 | ≈ **1.00** | ≈ 0.90 |
| 0.35 – 0.55 m | ≈ 0.88 | ≈ 0.98 | ≈ 0.91 |
| 0.55 – 0.75 m | ≈ 0.89 | ≈ **0.75** | ≈ 0.89 |

*Key Pattern*: Split conformal's spread across bins should be **≳ 0.20**; the GP's and normalized conformal's should be **≲ 0.03**. That contrast is the paper.

## Experiment 04 — aperture shift

Calibration inside the aperture, test region scaled outward. Split conformal collapses; the other two hold.

| Test region / aperture | GP | Split conformal | Normalized |
|---|---|---|---|
| 1.0 | ≈ 0.89 | ≈ 0.90 | ≈ 0.90 |
| 2.0 | ≈ 0.89 | ≈ **0.28** | ≈ 0.89 |

A split-conformal coverage of ~0.28 against a nominal 0.90, with no warning emitted, is the single most striking number in the project.

## Experiment 05 — frequency sweep

NMSE degrades with frequency as the fixed layout under-samples the field. The question is whether the intervals *notice*. They should: coverage stays near
nominal even as accuracy collapses, because the GP variance grows. If coverage falls while NMSE worsens, the intervals are over-confident exactly where a practitioner is most likely to be misled — report that prominently if it happens on your run.

## If something is off

| Symptom | Likely cause |
|---|---|
| GP coverage ≈ 1.00, huge radii | signal variance not fitted |
| All coverages ≈ 0.95 at 90% nominal | Rayleigh factor of 2 wrong in `gpRadius` |
| Split conformal slightly under nominal | interpolated quantile instead of order statistic |
| `chol` fails | duplicate microphone positions, or k too high for the layout |
| NMSE flat in microphone count | field and kernel wavenumbers disagree |


---

# Part II — the contribution (experiments 06–08)

## Experiment 06 — weighted conformal under designed density shift

Calibration microphones follow a truncated Gaussian of width `sigma_p` centred in the aperture; query positions are uniform over the same cube.
Nominal coverage 0.90.

| `sigma_p` | ESS / n | split CP | normalized CP | **weighted CP** |
|---|---|---|---|---|
| 0.18 | ≈ 0.12 | ≈ 0.58 | ≈ 0.93 | ≈ 0.93 |
| 0.25 | ≈ 0.37 | ≈ 0.73 | ≈ 0.92 | ≈ 0.92 |
| 0.35 | ≈ 0.77 | ≈ 0.83 | ≈ 0.91 | ≈ 0.91 |
| 0.50 | ≈ 0.94 | ≈ 0.86 | ≈ 0.90 | ≈ 0.90 |
| 0.80 | ≈ 0.99 | ≈ 0.89 | ≈ 0.90 | ≈ 0.90 |

**The pattern that matters:** split conformal degrades monotonically as the calibration density diverges from the query density (0.89 → 0.58), while weighted conformal stays at or above nominal throughout. At a representative mid setting the headline numbers are **0.765 unweighted vs 0.901 weighted**.

Slight over-coverage at very low ESS (≈ 0.93 rather than 0.90) is expected and correct: with few effective points the discrete weighted quantile rounds up.

If weighted conformal tracks split conformal instead of holding at nominal, the likelihood ratio is being computed the wrong way round — check the argument order in `sfuq.geom.likelihoodRatio(Pcal, Qquery, X)`.

## Experiment 07 — beyond the calibration support

Calibration uniform over the aperture; query region scaled outward.

| query/calib size | split CP coverage | weighted CP (certified pts) | abstention |
|---|---|---|---|
| 1.0 | ≈ 0.90 | ≈ 0.90 | ≈ 0.00 |
| 1.5 | ≈ 0.55 | ≈ 0.90 | ≈ 0.70 |
| 2.0 | ≈ 0.28 | ≈ 0.90 | ≈ 0.87 |

**The single most striking number in the project** is split conformal at 0.283 coverage against 0.900 nominal, with a 0.000 abstention rate — it never
declines to answer and is wrong 72% of the time. Weighted conformal abstains on 87% of those points and is correct on the remainder.

The certifiable fraction should track the volume ratio: at size ratio *r* it is
≈ 1/r³ (0.30 at r = 1.5, 0.125 at r = 2.0).

## Experiment 08 — placement

Three calibration layouts, identical microphone counts and model.

| design | ESS / n | weighted coverage | split coverage | width (rel.) |
|---|---|---|---|---|
| compact | ≈ 0.04–0.12 | ≈ 0.93 | ≈ 0.58 | ≈ 0.95–1.05 |
| uniform | ≈ 1.00 | ≈ 0.90 | ≈ 0.90 | ≈ 1.00 |
| matched | ≈ 1.00 | ≈ 0.90 | ≈ 0.90 | 1.00 |

Coverage is essentially identical across layouts —
weighted conformal is valid for all of them, which is the point of a distribution-free guarantee. What separates them is effective sample size, which collapses by an order of magnitude for a compact array. The cost of a convenient layout is paid in efficiency, not in coverage, and is invisible unless ESS is reported.

`matchedDesign` must give ESS/n = 1 to within 1e-10 — it draws from the query density, so every weight is exactly 1. If it does not, the ratio is not
cancelling and something is wrong in `likelihoodRatio`.

## New sanity checks

| Quantity | Expected |
|---|---|
| Weighted CP with all weights equal | identical to split conformal, to 1e-12 |
| `likelihoodRatio(P, P, X)` | exactly 1 for all X |
| `likelihoodRatio` for nested uniform boxes | vol(P)/vol(Q) inside, Inf outside |
| ESS of n equal weights | exactly n |
| ESS with any infinite weight | 0 |
| Dominating query weight (w ≫ Σw_cal) | radius = Inf, not the max score |

## If something is off (Part II)

| Symptom | Likely cause |
|---|---|
| Weighted coverage tracks split coverage | densities passed in the wrong order |
| Weighted coverage ≈ 1.00 everywhere | most radii are Inf — check abstention rate, not just coverage |
| ESS/n > 1 | weights not derived from a true density ratio |
| `infiniteCalWeight` error | Pcal and Qquery swapped in the call |
| Abstention 0 when query region is larger | query density support not actually exceeding calibration support |
