# Contributing

## Ground rules

1. **Never commit a number you did not compute.** Placeholder values must say so.
2. **The protocol precedes the result.** Changing `docs/protocol.md` gets its own
   commit explaining why, and every affected number is regenerated.
3. **Abstentions are never scored as successes.** An infinite interval covers by
   definition. Weighted conformal is scored only where it issued a finite
   interval, with the abstention rate reported alongside — always.
4. **Coverage is never reported without width and effective sample size.**
   Coverage alone is gamed by widening; width alone by narrowing.
5. **A failed reproduction is a finding.** Record it rather than hiding it.

## Before opening a PR

```matlab
run_tests
```

All six classes must pass. The suite runs on base MATLAB with no toolboxes and
no downloads; if it starts needing either, a dependency has leaked in and the
fix belongs in the code, not the test.

## Adding a density

Implement `.sample(n, stream)`, `.pdf(X)` and `.support(X)` following
`+sfuq/+geom/uniformBox.m`. The pdf may be unnormalised **only if** it will
always be paired against another density on the same support, so the constants
cancel in the ratio — set `.normalised = false` and say so in the header, as
`gaussianBox` does.

Then add a case to `tGeometry` checking the resulting likelihood ratio against
a closed form. The exactness of that ratio is the project's claim to novelty;
it is checked against mathematics, never against another implementation.

## Commit messages

State what changed and why it is correct, not which file you touched.
