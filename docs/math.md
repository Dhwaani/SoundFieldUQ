# The mathematics

Everything in the project reduces to four objects. This file states them
precisely so the code can be read against the maths. 

## 1. The Helmholtz kernel

A sound field $p(\mathbf{r})$ in a source-free region at wavenumber
$k=\omega/c$ satisfies the homogeneous Helmholtz equation

$$\nabla^2 p + k^2 p = 0.$$

Every such field can be written as a superposition of plane waves whose
wavevectors lie on the sphere of radius $k$. Placing a uniform measure on that
sphere and integrating gives the reproducing kernel

$$\kappa(\mathbf{r},\mathbf{r}') = j_0\!\big(k\|\mathbf{r}-\mathbf{r}'\|\big)
= \frac{\sin\big(k\|\mathbf{r}-\mathbf{r}'\|\big)}{k\|\mathbf{r}-\mathbf{r}'\|}.$$

Functions in the associated RKHS satisfy the Helmholtz equation **by
construction**. That is what makes it a far stronger prior than a generic
squared-exponential kernel at the same number of measurements, and
`exp01_accuracy` quantifies the gap.

Note $\kappa(\mathbf{r},\mathbf{r})=1$, so the kernel is normalised and the
signal power is carried separately (§3).

## 2. The synthetic field

For a shoebox room the image-source expansion gives, exactly,

$$p(\mathbf{r}) = \sum_i \beta^{n_i}\,
\frac{e^{-\mathrm{j}k\|\mathbf{r}-\mathbf{r}_i\|}}{4\pi\|\mathbf{r}-\mathbf{r}_i\|},$$

where $\mathbf{r}_i$ are mirrored source positions and $n_i$ the number of
reflections. Each term is a free-field Green's function and therefore solves
the Helmholtz equation away from its own singularity; so does the sum.

Working at a single wavenumber rather than synthesising a time-domain impulse
response avoids fractional-delay interpolation entirely, and gives ground
truth at **any** point — which is what makes exact held-out evaluation
possible without hardware.

## 3. Gaussian process with a fitted scale

Parameterise the covariance as

$$\mathbf{C} = s\big(\mathbf{K} + \sigma_{\mathrm{rel}}^2\mathbf{I}\big),$$

so a single factor $s$ multiplies signal and noise together. The type-II
maximum-likelihood estimate is then available in closed form:

$$\hat{s} = \frac{\mathbf{y}^{H}
\big(\mathbf{K}+\sigma_{\mathrm{rel}}^2\mathbf{I}\big)^{-1}\mathbf{y}}{2n},$$

the factor $2$ because each complex observation carries two real degrees of
freedom. No numerical optimisation is needed.

Posterior mean and per-component variance:

$$\mu(\mathbf{r}) = \mathbf{k}_*^{\!\top}\mathbf{C}_0^{-1}\mathbf{y},
\qquad
v(\mathbf{r}) = \hat{s}\big(1 - \mathbf{k}_*^{\!\top}\mathbf{C}_0^{-1}\mathbf{k}_*\big).$$

$\hat{s}$ cancels in the mean, which is why kernel ridge regression with
$\lambda=\sigma_{\mathrm{rel}}^2$ gives an identical point estimate — pinned by
`tModels/krrMeanMatchesGpMean`.

**Fitting $s$ is not a detail.** Leaving it at 1 while the field magnitude is
order $10^{-2}$ inflates every interval by more than an order of magnitude and
produces the false impression that GP uncertainty is worthless. Most of the
"GP is miscalibrated" folklore is this mistake.

## 4. Prediction regions

For complex pressure the natural region is a **disk**: $\{z : |z-\mu| \le R\}$.

**GP disk.** If the real and imaginary errors are each $\mathcal{N}(0,v)$ then
$|e|^2/v \sim \chi^2_2$, whose CDF is $1-e^{-x/2}$. Hence

$$R_{\mathrm{GP}} = \sqrt{-2v\log\alpha}.$$

Equivalently $|e|$ is Rayleigh with scale $\sqrt{v}$. The factor of two is the
easiest thing in the project to get wrong; `tConformal` pins it by simulation.

**Split conformal.** With calibration scores $s_i = |y_i - \mu(\mathbf{x}_i)|$
on a set disjoint from the training set,

$$R_{\mathrm{CP}} = s_{(\lceil (n+1)(1-\alpha)\rceil)},$$

the order statistic — *not* an interpolated quantile. Then, for exchangeable
calibration and test points,

$$\mathbb{P}\big(|y-\mu(\mathbf{x})| \le R_{\mathrm{CP}}\big) \ge 1-\alpha$$

for **any** model, kernel and data distribution, at finite $n$. The $(n+1)$
accounts for the unseen test point joining the exchangeable set.

**Normalized conformal.** Take scores $s_i = |y_i-\mu(\mathbf{x}_i)|/\sigma(\mathbf{x}_i)$
and set

$$R_{\mathrm{nCP}}(\mathbf{x}) = \hat{q}\,\sigma(\mathbf{x}).$$

Marginal validity is retained (conformal is agnostic to the score); the radius
now follows the model's own sense of difficulty. With $\sigma$ from the
Helmholtz GP, the normalizer is physics-informed.

## 6. Weighted conformal with an exactly known ratio (the contribution)

### The shift

Calibration points are drawn from a density $p$ over a region; the positions to
be certified are drawn from a density $q$. In a real deployment these differ:
microphones sit where the array is, listeners sit where the listeners are.
Exchangeability fails and split conformal is no longer valid.

### Why acoustics is special

Weighted conformal prediction (Tibshirani, Barber, Candès & Ramdas, 2019)
restores marginal coverage *provided the likelihood ratio*

$$w(\mathbf{r}) = \frac{\mathrm{d}Q}{\mathrm{d}P}(\mathbf{r}) = \frac{q(\mathbf{r})}{p(\mathbf{r})}$$

*is known.* In almost every application it is not, and must be estimated from
unlabelled data — the dominant error source.

Here the covariate is position and **both densities are designed**. For two
uniform boxes, $w = \mathrm{vol}(P)/\mathrm{vol}(Q)$ on the overlap. For a
truncated-Gaussian calibration density against a uniform query density on the
same box,

$$w(\mathbf{r}) \;\propto\; \exp\!\left(\frac{\|\mathbf{r}-\mathbf{c}\|^2}{2\sigma^2}\right),$$

the normalising constants cancelling in the ratio. No estimation, no
estimation error.

### The construction

For calibration scores $s_i$ with weights $w_i = w(\mathbf{r}_i)$ and a query
point with weight $w(\mathbf{r})$, each score carries normalised mass

$$p_i(\mathbf{r}) = \frac{w_i}{\sum_j w_j + w(\mathbf{r})},
\qquad
p_\infty(\mathbf{r}) = \frac{w(\mathbf{r})}{\sum_j w_j + w(\mathbf{r})},$$

and the radius is the $(1-\alpha)$ quantile of $\sum_i p_i \delta_{s_i} + p_\infty \delta_{+\infty}$.
The atom at $+\infty$ plays the role the $(n+1)$ plays in the unweighted case:
it accounts for the unseen query point joining the weighted exchangeable set,
which is what makes the guarantee finite-sample.

Setting all $w_i = 1$ recovers split conformal exactly — pinned by
`tWeightedConformal/reducesToSplitConformalWhenWeightsAreEqual`.

### Absolute continuity, and abstention

The guarantee requires $Q \ll P$: wherever $q > 0$, we need $p > 0$. If a
listener position lies outside the region the calibration microphones cover,
then $w = +\infty$ there and **no reweighting can help** — there is simply no
measurement constraining that point.

The correct behaviour is then $p_\infty > \alpha$, the quantile is $+\infty$,
and the interval is infinite. The method **abstains**.

This is the project's most useful finding. Split conformal in the same
situation returns a confident finite interval and collapses to $\approx 0.28$
empirical coverage against $0.90$ nominal, *silently*. An explicit "I cannot
certify this point" is worth far more than a confident wrong answer, and the
abstention map is computable from geometry before any measurement is taken.

### The price: effective sample size

Weighting costs efficiency. With Kish effective sample size

$$\mathrm{ESS} = \frac{\left(\sum_i w_i\right)^2}{\sum_i w_i^2},$$

the quantile is effectively estimated from $\mathrm{ESS}$ points rather than
$n$. As $p$ and $q$ diverge, ESS falls and intervals widen even though coverage
stays correct. Reporting ESS alongside coverage is what stops the method being
presented as free.

### Placement follows

Since $w = q/p$ and ESS is maximised exactly when $p = q$:

Put the calibration microphones where the listeners are. Existing placement
criteria — posterior variance, mutual information, conditional entropy —
optimise how much the model *learns*, which is the right objective for fitting
a field and the wrong one for certifying it. The two optima differ, and
`exp08_placement.m` measures the gap: all layouts are valid, but a compact
array pays for its convenience in width, not in coverage, which is exactly why
the cost is easy to miss.

## 7. What is guaranteed and what is not

| Property | GP | Split CP | Normalized CP | **Weighted CP** |
|---|---|---|---|---|
| Marginal coverage, no shift | if model correct | **guaranteed** | **guaranteed** | **guaranteed** |
| Marginal coverage, under shift | approx. | **fails** | approx. | **guaranteed** (exact $w$) |
| Conditional coverage | approx. | no | approx. | approx. |
| Beyond calibration support | degrades | **fails silently** | degrades | **abstains** |
| Needs a variance estimate | yes | no | yes | no |
| Needs a known $\mathrm{d}Q/\mathrm{d}P$ | no | no | no | yes — and here it is exact |

Conditional validity is provably unattainable in full generality without
distributional assumptions (Vovk; Lei & Wasserman; Barber, Candès, Ramdas &
Tibshirani 2021). Weighting fixes the *shift*, not the *conditionality*; the
two are addressed by different mechanisms and the paper must say so.

## 8. Legacy comparison (unweighted methods only)

### 5. What is guaranteed and what is not

| Property | GP | Split conformal | Normalized conformal |
|---|---|---|---|
| Marginal coverage | if model correct | **guaranteed** | **guaranteed** |
| Conditional (spatial) coverage | approx., empirically good | **no** — constant radius | approx., empirically good |
| Robust to covariate shift | degrades gracefully | **fails, silently** | degrades gracefully |
| Needs a variance estimate | yes | no | yes |

Conditional validity is provably unattainable in full generality without distributional assumptions (Vovk; Lei & Wasserman; Barber et al. on the limits of conditional conformal inference). Normalized conformal is a practical approximation, not a solution, and the paper must say so.

## 9. Threats to validity

- **Synthetic-only.** Rigid walls, frequency-independent absorption, no
  scattering or diffraction from furnishings. Real rooms are messier; the
  measured-data extension is the check.
- **Single frequency per experiment.** Broadband coverage is a separate and
  harder question, since errors correlate across frequency.
- **Noise is circular complex Gaussian.** Real measurement chains are not.
  Conformal is distribution-free and should be indifferent; the GP is not, and
  that is a point in conformal's favour that these experiments do not yet
  exercise.
- **Distance from the region centre** is a proxy for "distance from the
  supporting microphones". A local-density covariate would be sharper.
- **The designed densities are idealisations.** Real microphone positions carry
  placement error, so the "exact" ratio is exact only up to how precisely the
  positions are known. Quantifying sensitivity of coverage to position error is
  open work.
- **ESS is a summary, not a guarantee.** A low ESS warns that intervals are
  estimated from few effective points; it does not by itself bound the error in
  the quantile.
