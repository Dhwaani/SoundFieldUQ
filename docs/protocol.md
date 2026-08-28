# Experimental protocol

Written before the results and committed ahead of them. If a decision here
turns out to be wrong it changes in its own commit, saying so, and every
affected number is regenerated. A protocol that moves after the numbers arrive
is not measuring anything.

## The claim under test

> For sparse sound-field interpolation, spatial uncertainty estimates from standard methods aren't equally reliable across space. Specifically, split conformal prediction—despite having the strongest theoretical guarantee—is the least trustworthy conditionally.

What is **not** claimed: that any method reconstructs more accurately; that
these results transfer to measured rooms without checking; that normalized
conformal is valid under covariate shift.

## Data discipline

Three disjoint, independently drawn sets per trial:

| Set | Used for | Must not be used for |
|---|---|---|
| microphones | fitting the GP | anything else |
| calibration | conformal quantiles | fitting |
| test | evaluation | fitting or calibration |

Reusing training points for calibration makes residuals optimistically small
and the radii too narrow. The guarantee is voided, and the numbers still look
plausible. This is the single easiest way to produce a wrong result here.

Accuracy (NMSE) is computed against the **noiseless** ground truth; coverage
against the **noisy observation**, since that is what a measurement would
actually return. Both are available because the field is analytic.

## Fixed conditions

- Room 6 × 5 × 3 m, wall reflection coefficient 0.7, image-source order 3.
- Measurement aperture: 1 m cube at the room centre.
- Measurement noise: circular complex Gaussian at 2% of field RMS.
- Nominal coverage: 90% (α = 0.10) unless a sweep says otherwise.
- Trials per condition: 8–12, seeds fixed and derived from the condition index
  so every number is reproducible.

## Reporting rules

- Coverage is always reported **with** interval width or interval score.
  Coverage alone is gamed by making regions enormous; width alone by making
  them zero. The Winkler interval score resolves this and is reported.
- Conditional coverage bins with fewer than 20 points report `NaN`, never a
  coverage estimated from a handful of samples.
- Marginal and conditional coverage are always reported together. Reporting
  only the marginal number would hide the paper's own main finding.
- Negative results are reported first. Experiment 02 shows the GP is already
  well calibrated marginally; that goes in the abstract, not a footnote.

## Sequence

1. `exp01` — accuracy. If NMSE does not fall with microphone count, or the
   Helmholtz kernel does not beat a generic RBF, stop: something is wrong
   upstream and no downstream calibration analysis would mean anything.
2. `exp02` — marginal calibration. Establishes the comparison is fair.
3. `exp03` — conditional coverage. The result.
4. `exp04` — aperture shift. Test boundary failure modes.
5. `exp05` — frequency sweep. Separates expected physical degradation from
   genuine over-confidence.
6. `exp06` — weighted conformal under a designed density shift. The main result.
7. `exp07` — beyond the calibration support. Document failure limits and expected abstention behavior.
8. `exp08` — placement. The corollary, and the practically useful output.

## Known threats

Listed in `math.md` §6.
