function [R, info] = weightedConformal(residCal, wCal, wQuery, alpha)
%WEIGHTEDCONFORMAL  Covariate-shift conformal prediction with EXACT weights.
%
%   [R, INFO] = WEIGHTEDCONFORMAL(RESIDCAL, WCAL, WQUERY, ALPHA) returns a
%   per-query-point prediction radius R that retains finite-sample marginal
%   coverage 1-ALPHA when the calibration and query covariate distributions
%   differ, provided the weights are the true likelihood ratio dQ/dP.
%
%   Inputs
%     residCal  n-by-1 absolute residuals on a calibration set disjoint from
%               the training set
%     wCal      n-by-1 likelihood ratios at the calibration points
%     wQuery    m-by-1 likelihood ratios at the query points
%     alpha     target miscoverage, e.g. 0.10
%
%   Outputs
%     R     m-by-1 radii; Inf marks a query point that cannot be certified
%     INFO  struct with .ess, .essFraction, .abstentionRate, .qHat, .nFinite
%
%   The construction
%   ----------------
%   Following Tibshirani, Barber, Candes & Ramdas (NeurIPS 2019), each
%   calibration score carries normalised mass
%
%       p_i(x) = w_i / ( sum_j w_j + w(x) ),
%
%   together with a point mass at +infinity of
%
%       p_inf(x) = w(x) / ( sum_j w_j + w(x) ).
%
%   The radius is the (1-alpha) quantile of that weighted distribution. The
%   extra atom at infinity is what makes the guarantee finite-sample rather
%   than asymptotic: it plays the same role as the (n+1) in unweighted split
%   conformal, accounting for the unseen query point joining the weighted
%   exchangeable set. Note that the radius depends on the query point through
%   w(x) in the denominator, so unlike split conformal it is not constant.
%
%   When the atom at infinity alone exceeds alpha the quantile IS infinity, and
%   the method abstains. That happens exactly when a query point is far more
%   likely under Q than anything the calibration set represents -- which is the
%   honest answer, not a failure.
%
%   Efficiency, reported rather than buried
%   ---------------------------------------
%   Weighting costs effective sample size. INFO.ESS is the Kish effective
%   sample size of the calibration weights; as the calibration and query
%   densities diverge it falls, the quantile is estimated from fewer effective
%   points, and the intervals widen. Reporting it alongside coverage is what
%   stops the method being presented as a free lunch, and it is also the
%   objective that SFUQ.DESIGN.GREEDYESS optimises when choosing microphone
%   positions.
%
%   See also SFUQ.GEOM.LIKELIHOODRATIO, SFUQ.UQ.SPLITCONFORMAL,
%   SFUQ.DESIGN.GREEDYESS.

    residCal = abs(residCal(:));
    wCal     = wCal(:);
    wQuery   = wQuery(:);

    validateattributes(alpha, {'numeric'}, {'scalar', '>', 0, '<', 1});
    if numel(residCal) ~= numel(wCal)
        error('sfuq:weightedConformal:sizeMismatch', ...
              'residCal (%d) and wCal (%d) must have the same length.', ...
              numel(residCal), numel(wCal));
    end
    if any(~isfinite(wCal))
        error('sfuq:weightedConformal:infiniteCalWeight', ...
              ['A calibration point has infinite weight, meaning the query ' ...
               'density is positive where the calibration density is zero AT ' ...
               'A CALIBRATION POINT. That is contradictory; check that the ' ...
               'densities are the right way round.']);
    end

    % Zero-weight calibration points contribute nothing; drop them so they do
    % not distort the effective sample size.
    keep     = wCal > 0;
    residCal = residCal(keep);
    wCal     = wCal(keep);

    if isempty(wCal)
        R = inf(size(wQuery));
        info = struct('ess', 0, 'essFraction', 0, 'abstentionRate', 1, ...
                      'qHat', Inf, 'nFinite', 0, 'nCal', 0);
        return
    end

    [scores, ord] = sort(residCal);
    wSorted = wCal(ord);
    cumW    = cumsum(wSorted);
    total   = cumW(end);

    % A query point needs cumulative calibration mass >= (1-alpha) of the
    % total INCLUDING its own atom. Rearranged, the threshold on the raw
    % cumulative weight is (1-alpha) * (total + w(x)).
    thr = (1 - alpha) * (total + wQuery);

    R = inf(size(wQuery));
    finiteQ = isfinite(wQuery);

    if any(finiteQ)
        idx = firstIndexAtLeast(cumW, thr(finiteQ));
        r = inf(sum(finiteQ), 1);
        ok = idx <= numel(scores);
        r(ok) = scores(idx(ok));
        R(finiteQ) = r;
    end

    info = struct();
    info.nCal           = numel(wCal);
    info.ess            = total^2 / sum(wSorted.^2);
    info.essFraction    = info.ess / numel(wCal);
    info.abstentionRate = mean(~isfinite(R));
    info.nFinite        = sum(isfinite(R));
    info.qHat           = median(R(isfinite(R)));
    if isempty(info.qHat), info.qHat = Inf; end
end

function idx = firstIndexAtLeast(cumW, thr)
%   For each threshold, the index of the first cumulative weight that reaches
%   it. cumW is ascending, so this is a search; returns numel(cumW)+1 when no
%   index qualifies, which the caller maps to an infinite radius.
%
%   Chunked to bound peak memory: the naive outer comparison is
%   numel(cumW)-by-numel(thr), which is fine for a few hundred thousand
%   elements and ruinous beyond that.

    n = numel(cumW);
    m = numel(thr);
    idx = zeros(m, 1);

    chunk = max(1, floor(4e6 / max(n, 1)));
    for a = 1:chunk:m
        b = min(a + chunk - 1, m);
        idx(a:b) = sum(cumW(:) < thr(a:b).', 1).' + 1;
    end
end
