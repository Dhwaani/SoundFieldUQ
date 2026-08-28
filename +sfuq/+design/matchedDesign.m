function X = matchedDesign(Qquery, n, stream)
%MATCHEDDESIGN  The provably ESS-optimal calibration layout: sample from Q.
%
%   X = MATCHEDDESIGN(QQUERY, N) returns N calibration positions drawn from the
%   query density itself.
%
%   This is the continuum optimum. The conformal weights are w = q/p, so when
%   p = q every weight is one, the Kish effective sample size equals n exactly,
%   and the weighted quantile is estimated as efficiently as an unweighted one.
%   No layout can do better on that criterion.
%
%   It is included as the reference against which practical layouts are scored
%   in EXP08_PLACEMENT: a compact array (calibration clustered at the centre)
%   loses most of its effective sample size, while a matched layout retains all
%   of it. The gap between the two is the cost of putting the microphones where
%   they are convenient rather than where the listeners are.
%
%   Note the deliberate implication: if the listening area is where you must
%   certify, then that is where the calibration microphones belong, even though
%   for FITTING the field a compact array may be perfectly adequate. Fitting and
%   certifying have different optimal designs, and conflating them is the
%   mistake this function exists to make visible.

    if nargin < 3, stream = []; end
    X = Qquery.sample(n, stream);
end
