function [Rq, qHat] = normalizedConformal(residCal, sigmaCal, sigmaQ, alpha)
%NORMALIZEDCONFORMAL  Locally adaptive conformal prediction radii.
%
%   [RQ, QHAT] = NORMALIZEDCONFORMAL(RESIDCAL, SIGMACAL, SIGMAQ, ALPHA)
%   returns a per-query-point radius RQ = QHAT * SIGMAQ, where QHAT is the
%   conformal quantile of the NORMALIZED scores |y - mu(x)| / sigma(x).
%
%   Inputs
%     residCal  absolute residuals on the calibration set
%     sigmaCal  predictive standard deviation at the calibration points
%     sigmaQ    predictive standard deviation at the query points
%     alpha     target miscoverage
%
%   The idea
%   --------
%   Split conformal gives a guarantee but a constant radius. A model that
%   reports its own uncertainty gives a shape but no guarantee. Dividing the
%   residual by the model's own sigma before taking the conformal quantile
%   keeps the distribution-free marginal guarantee while letting the radius
%   follow the model's sense of where the problem is hard.
%
%   In this project sigma comes from the Helmholtz GP posterior, so the
%   normalizer is physics-informed: it grows away from the microphones and
%   beyond the aperture. The combination is the method this project
%   recommends, and the experiments show it is the only one of the three that
%   is simultaneously (i) marginally valid, (ii) approximately conditionally
%   valid across space, and (iii) robust when calibration and test points come
%   from different regions.
%
%   Note that (iii) is empirical, not guaranteed. Normalization does not
%   repair a violation of exchangeability; it merely degrades far more
%   gracefully, because the normalizer itself is aware that extrapolated
%   points are harder. Do not overclaim this in the paper.
%
%   See also SFUQ.UQ.SPLITCONFORMAL.

    residCal = abs(residCal(:));
    sigmaCal = sigmaCal(:);
    sigmaQ   = sigmaQ(:);

    if numel(residCal) ~= numel(sigmaCal)
        error('sfuq:normalizedConformal:sizeMismatch', ...
              'residCal and sigmaCal must have the same length.');
    end

    % Floor the normalizer. A sigma of exactly zero (a query point sitting on
    % a microphone) would make the score infinite and drag the quantile to the
    % maximum, destroying the radii everywhere else.
    floorVal = 1e-12 * max(sigmaCal);
    sigmaCal = max(sigmaCal, floorVal);
    sigmaQ   = max(sigmaQ,   floorVal);

    scores = residCal ./ sigmaCal;

    n = numel(scores);
    level = min(ceil((n + 1) * (1 - alpha)) / n, 1);
    qHat = sfuq.uq.quantileHigher(scores, level);

    Rq = qHat * sigmaQ;
end
