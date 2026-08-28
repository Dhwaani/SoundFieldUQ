function [w, isCertifiable] = likelihoodRatio(Pcal, Qquery, X)
%LIKELIHOODRATIO  Exact covariate-shift weight w(x) = dQ/dP at points X.
%
%   [W, CERT] = LIKELIHOODRATIO(PCAL, QQUERY, X) returns the ratio of the query
%   density to the calibration density at each row of X, together with a
%   logical flag marking the points where that ratio is finite.
%
%   This function is the heart of the project's contribution.
%
%   Weighted conformal prediction (Tibshirani, Barber, Candes & Ramdas, 2019)
%   restores marginal coverage under covariate shift provided the likelihood
%   ratio dQ/dP is known. In essentially every published application it is NOT
%   known and must be estimated from unlabelled data, which is the dominant
%   source of error and the reason the method is not used more widely.
%
%   Sound-field reconstruction is different. The covariate is spatial position,
%   and BOTH distributions are chosen by the experimenter: where the
%   calibration microphones go, and which listener positions are to be
%   certified. The ratio is therefore available in closed form from geometry
%   alone, with no estimation and no estimation error.
%
%   The absolute-continuity condition, and why it is reported not hidden
%   -------------------------------------------------------------------
%   Weighted conformal requires Q << P: wherever the query density is positive,
%   the calibration density must be too. If a listener position lies outside
%   the region the calibration microphones cover, the ratio there is genuinely
%   infinite and NO reweighting can rescue it. The correct behaviour is to
%   return Inf, which makes the downstream prediction interval infinite -- an
%   explicit abstention.
%
%   That is a feature. Plain split conformal in the same situation returns a
%   confident, finite, and badly wrong interval; EXP07_SUPPORT_FAILURE measures
%   it collapsing to roughly 0.28 empirical coverage against a nominal 0.90,
%   silently. An honest "I cannot certify this point" is worth far more to a
%   practitioner than a confident lie, and CERT is exactly the map of which
%   listener positions the current microphone layout can and cannot certify.
%
%   Outputs
%     W     N-by-1 weights, Inf where the calibration density vanishes
%     CERT  N-by-1 logical, true where W is finite (the certifiable points)
%
%   See also SFUQ.UQ.WEIGHTEDCONFORMAL, SFUQ.DESIGN.CERTIFIABLEFRACTION.

    p = Pcal.pdf(X);
    q = Qquery.pdf(X);

    w = zeros(size(p));

    live = p > 0;
    w(live) = q(live) ./ p(live);

    % q > 0 where p == 0: genuinely outside the certifiable region.
    w(~live & q > 0) = Inf;

    % q == 0 and p == 0: the point is irrelevant to the query distribution.
    % Weight zero, and it simply drops out of the weighted quantile.
    w(~live & q <= 0) = 0;

    isCertifiable = isfinite(w);
end
