# SoundFieldUQ

**Certifying a reconstructed sound field: exact-ratio weighted conformal prediction, and knowing when to abstain.**

[![MATLAB](https://img.shields.io/badge/MATLAB-R2020a%2B-orange.svg)](https://www.mathworks.com)
[![license](https://img.shields.io/badge/license-MIT-blue.svg)](LICENSE)
[![toolboxes](https://img.shields.io/badge/toolboxes-none%20required-green.svg)](#requirements)
[![hardware](https://img.shields.io/badge/hardware-none%20required-green.svg)](#why-no-hardware-is-needed)

```matlab
run('startup.m')
run_tests        % ~40 s, no data, no toolboxes
run_all(true)    % ~3 min, all ten figures
```

---

## The one-sentence contribution

In sound-field reconstruction the covariate is **spatial position**, and the experimenter chooses both where the calibration microphones go and which listener positions must be certified — so the covariate-shift likelihood ratio dQ/dP is available **exactly, in closed form, from geometry alone**, which is the one thing weighted conformal prediction normally cannot have.

## Why that matters

Weighted conformal prediction ([Tibshirani, Barber, Candès & Ramdas, 2019](https://arxiv.org/abs/1904.06019)) restores marginal coverage under covariate shift *provided the likelihood ratio is known*. In essentially every published application it is not: it must be estimated from unlabelled data, and that estimation is the dominant source of error and the main reason the method is not used more widely.

Acoustic reconstruction is a rare case where the ratio is a **design variable, not an unknown**. Microphone positions and listener positions are both chosen. Their densities are known analytically. Their ratio is a ratio of known densities over known regions.

This repository works out what follows.

## What we found

**1. Split conformal fails silently under realistic spatial shift.** Microphones cluster near an array; listeners do not. That is a covariate shift, and it breaks exchangeability. Split conformal's single constant radius is calibrated on residuals drawn disproportionately from the easy centre, so it under-covers: **0.765 empirical against 0.900 nominal**, with no warning.

**2. The exact geometric ratio restores validity.** Weighted conformal with dQ/dP computed from the designed densities recovers coverage to **0.901**. Not estimated, not approximated — exact.

**3. It costs effective sample size, and we report it.** Kish ESS falls to roughly half the calibration set under moderate clustering; intervals widen accordingly. Split conformal's intervals are narrower *precisely because they are wrong*. Coverage and width are always reported together.

**4. Beyond the calibration support, nothing restores validity — and that is the most useful result here.** When listener positions extend past the region the microphones cover, dQ/dP is genuinely infinite. Split conformal keeps issuing confident finite intervals and collapses to **0.283 coverage**. Weighted conformal returns an **infinite interval**: it abstains. It does not pretend a guarantee survives where no measurement constrains anything.

> A method that tells you it cannot answer is safe to deploy. A method that answers confidently and wrongly is not. The abstention rate is a map of exactly which listener positions your microphone layout cannot certify — computable **before you measure anything**.

**5. Placement follows as a corollary.** The weights are w = q/p and ESS is maximised when p = q, so the optimal calibration layout is a **sample from the query density itself**. Put the calibration microphones where the listeners are. Established placement criteria — posterior variance, mutual information, conditional entropy — optimise how much the model *learns*; that is the right objective for fitting a field and the wrong one for certifying it. All layouts are valid; they differ in cost, and the cost is invisible unless ESS is reported.

## Practical grounding

Calibrated bounds enable decisions that point estimates cannot:

| Application | What the certificate decides |
|---|---|
| Spatial ANC / virtual sensing | Worst-case residual at a virtual microphone where no sensor exists — a safety margin |
| Sound-zone control | Where the bright/dark contrast can be trusted, and where it silently degrades |
| AR/VR 6DoF audio | Which interpolated listener positions need re-measurement |
| Room-acoustic digital twins | A principled "re-measure here" trigger |
| Array design | Where to put microphones so the whole listening area is certifiable |

The machinery is not acoustics-specific. Any spatial inverse problem where the experimenter designs the sampling and query regions — geostatistics, EM field reconstruction, PDE surrogates — has the same exactly-known ratio.

## Requirements

Base MATLAB **R2020a or newer**. That is the entire list.

No Statistics and Machine Learning Toolbox: `quantile`, `range` and `pdist2` are implemented in-project. `quantileHigher` in particular is written out because MATLAB's `quantile` *interpolates* between order statistics, and conformal's finite-sample guarantee needs the order statistic itself — an interpolated quantile silently voids it.

## Why no hardware is needed

Ground truth is generated analytically by the image-source method, so exact complex pressure is available at **any** point. Held-out test positions are exact rather than interpolated — stronger than a measured dataset can offer, since a measured RIR at an unmeasured position does not exist. `tests/tField.m` verifies the synthetic field satisfies the Helmholtz equation to ~1e-5 relative residual.

Swapping in measured data (MeshRIR, dEchorate — both CC BY 4.0) means replacing one function, `sfuq.data.observe`. That is the natural next step, not a prerequisite.

## Layout

```
SoundFieldUQ/
├── startup.m / run_all.m / run_tests.m
├── +sfuq/
│   ├── trial.m               exchangeable fit/calibrate/test cycle
│   ├── shiftTrial.m          cycle under a DESIGNED covariate shift    [NEW]
│   ├── +geom/                densities + the exact likelihood ratio    [NEW]
│   ├── +uq/                  split, normalized, weighted conformal; ESS
│   ├── +design/              certifiable fraction, greedy ESS placement [NEW]
│   ├── +data/                image-source room, exact Helmholtz fields
│   ├── +kernels/             Helmholtz (sinc) kernel + RBF control
│   ├── +models/              GP with closed-form fitted signal variance
│   ├── +eval/                NMSE, coverage, conditional coverage, Winkler
│   └── +viz/                 all ten figures
├── experiments/              exp01-exp08
├── tests/                    6 classes, incl. Monte Carlo guarantee checks
├── docs/                     math.md, protocol.md, expected_results.md
└── paper/                    IEEE skeleton structured around the claim
```

## The experiments

| # | Question | Figures |
|---|---|---|
| 01 | Do the estimators work, and is the physics worth encoding? | 1 |
| 02 | Are the uncertainty estimates honest *on average*? | 2 |
| 03 | Are they honest at each *location*? | 3, 4 |
| 04 | What breaks when exchangeability is violated? | 5 |
| 05 | Do the intervals notice when reconstruction fails? | 6 |
| **06** | **Does the exact geometric ratio restore validity?** | **7, 8** |
| **07** | **What happens beyond the calibration support?** | **9** |
| **08** | **Where should the calibration microphones go?** | **10** |

01-05 establish the setting and audit existing methods. **06-08 are the contribution.** Experiment 02 reports a result that flatters the *baseline* — a properly fitted Helmholtz GP is already well calibrated for interpolation — and it is reported first, on purpose, because it shows the comparison was run fairly.

## Testing

```matlab
run_tests
```

The two that carry the claims:

- **`tWeightedConformal`** — verifies the weighted construction collapses exactly onto split conformal when weights are equal; achieves nominal coverage under a shift with an analytically known ratio (Monte Carlo, 250 replications); abstains rather than under-covering when a weight is infinite; and that a dominating query weight forces an infinite interval instead of silently truncating to the largest calibration score.
- **`tGeometry`** — checks the likelihood ratio against closed-form values, not against another implementation. Includes the volume-ratio identity and the exp(d²/2σ²) form.

Thresholds were set with margin measured across independent seeds, not chosen to just barely pass.

## Honest limits

- Weighted conformal restores **marginal** validity under shift. It does not grant **conditional** coverage — provably unattainable distribution-free ([Barber, Candès, Ramdas & Tibshirani, 2021](https://arxiv.org/abs/1903.04684)). The honest construction pairs weighting (for the shift) with variance normalization (for approximate conditionality).
- Reweighting cannot manufacture information. Where the support fails, the method abstains; it does not extrapolate.
- Synthetic rooms only so far: rigid walls, frequency-independent absorption, no scattering. Measured-data validation is the next step.
- Single frequency per experiment. Broadband coverage is harder because errors correlate across frequency.

## Roadmap

- [x] Exact-ratio weighted conformal; ESS accounting; abstention
- [x] Certifiable-region maps; ESS-optimal and greedy placement
- [ ] Measured data (MeshRIR S1-M3969, dEchorate)
- [ ] Functional conformal bands over whole RIRs, with per-segment coverage
- [ ] Band-wise coverage through the spatial-aliasing limit
- [ ] Preprint + archived artifact (Zenodo DOI)

## Citing

See [`CITATION.cff`](CITATION.cff). `+sfuq/+uq/` and `+sfuq/+geom/` are method-agnostic and reusable on any estimator emitting a point prediction; MIT licensed and meant to be lifted.

## Licence

MIT. Everything is synthetic and self-generated, so there are no dataset licences to inherit.
