function c = coverage(residuals, radii)
%COVERAGE  Empirical coverage of a set of prediction regions.
%
%   C = COVERAGE(RESIDUALS, RADII) returns the fraction of points whose
%   absolute residual falls inside the corresponding radius. RADII may be a
%   scalar (a constant-width region) or a vector the same length as RESIDUALS.
%
%   This is the quantity every uncertainty claim in the project is judged
%   against. A method advertising 90% intervals should score 0.90 here; the
%   distance from 0.90 is the calibration error.

    residuals = abs(residuals(:));
    radii = radii(:);
    if isscalar(radii)
        radii = repmat(radii, numel(residuals), 1);
    elseif numel(radii) ~= numel(residuals)
        error('sfuq:coverage:sizeMismatch', ...
              'radii must be scalar or the same length as residuals.');
    end

    c = mean(residuals <= radii);
end
