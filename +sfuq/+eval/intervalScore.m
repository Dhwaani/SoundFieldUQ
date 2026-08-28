function s = intervalScore(residuals, radii, alpha)
%INTERVALSCORE  Winkler interval score: sharpness penalised by miscoverage.
%
%   S = INTERVALSCORE(RESIDUALS, RADII, ALPHA) returns the mean score
%
%       2*R + (2/alpha) * max(|e| - R, 0)
%
%   Lower is better.
%
%   Coverage alone can be gamed by making every interval enormous, and width
%   alone by making them all zero. The interval score is the standard proper
%   scoring rule that resolves this: it rewards narrow regions but charges a
%   steep penalty, scaled by 1/alpha, whenever the truth falls outside. A
%   method should be judged on coverage AND this score together, never on
%   width in isolation.
%
%   The factor two makes this the two-sided form, consistent with reporting a
%   radius rather than a half-width.

    residuals = abs(residuals(:));
    radii = radii(:);
    if isscalar(radii)
        radii = repmat(radii, numel(residuals), 1);
    end

    penalty = (2 / alpha) * max(residuals - radii, 0);
    s = mean(2 * radii + penalty);
end
