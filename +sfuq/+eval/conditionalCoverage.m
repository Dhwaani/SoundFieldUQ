function T = conditionalCoverage(residuals, radii, covariate, edges)
%CONDITIONALCOVERAGE  Coverage broken down by a covariate, e.g. distance.
%
%   T = CONDITIONALCOVERAGE(RESIDUALS, RADII, COVARIATE, EDGES) bins the test
%   points by COVARIATE and reports coverage within each bin.
%
%   Returned struct array T has fields:
%     .lo, .hi     bin edges
%     .n           points in the bin
%     .coverage    empirical coverage in the bin
%     .meanRadius  mean prediction radius in the bin
%
%   Why this is the important measurement
%   -------------------------------------
%   Marginal coverage is an average, and averages hide structure. A method can
%   hit 90% overall while covering 100% near the microphones and 60% at the
%   edge of the aperture. For a practitioner deciding whether to trust a
%   reconstruction at a particular seat in a room, the conditional number is
%   the one that matters, and it is precisely the one that split conformal
%   does not control.
%
%   Bins with fewer than 20 points report NaN rather than a meaningless
%   coverage estimated from a handful of samples.

    residuals = abs(residuals(:));
    radii = radii(:);
    if isscalar(radii)
        radii = repmat(radii, numel(residuals), 1);
    end
    covariate = covariate(:);

    nb = numel(edges) - 1;
    T = repmat(struct('lo', NaN, 'hi', NaN, 'n', 0, ...
                      'coverage', NaN, 'meanRadius', NaN), nb, 1);

    for b = 1:nb
        if b < nb
            sel = covariate >= edges(b) & covariate < edges(b + 1);
        else
            sel = covariate >= edges(b) & covariate <= edges(b + 1);
        end

        T(b).lo = edges(b);
        T(b).hi = edges(b + 1);
        T(b).n  = sum(sel);

        if T(b).n >= 20
            T(b).coverage   = mean(residuals(sel) <= radii(sel));
            T(b).meanRadius = mean(radii(sel));
        end
    end
end
