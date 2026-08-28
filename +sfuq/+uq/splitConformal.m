function R = splitConformal(residCal, alpha)
%SPLITCONFORMAL  Distribution-free prediction radius by split conformal.
%
%   R = SPLITCONFORMAL(RESIDCAL, ALPHA) returns the scalar radius R such that
%   the prediction region {z : |z - mu(x)| <= R} covers the truth with
%   probability at least 1 - ALPHA.
%
%   Inputs
%     residCal  vector of absolute residuals |y - mu(x)| on a CALIBRATION set
%               that was NOT used to fit the model
%     alpha     target miscoverage, e.g. 0.10 for 90% coverage
%
%   Guarantee and its price
%   -----------------------
%   Provided the calibration and test points are exchangeable, this radius
%   satisfies P(|y - mu(x)| <= R) >= 1 - alpha for ANY model, ANY kernel and
%   ANY data distribution, at finite sample size. Nothing needs to be Gaussian
%   and the model does not need to be correct.
%
%   The price is that the guarantee is MARGINAL, averaged over test points. A
%   constant radius spends the same uncertainty budget everywhere, so it
%   over-covers where the field is easy and under-covers where it is hard.
%   Experiment 03 measures exactly that, and it is the reason
%   SFUQ.UQ.NORMALIZEDCONFORMAL exists.
%
%   Exchangeability is a real assumption, not a formality. Experiment 04
%   breaks it deliberately by calibrating inside the microphone aperture and
%   testing outside it, where coverage collapses.
%
%   For complex pressure the score |y - mu| is the modulus, so the region is a
%   disk in the complex plane. Conformal prediction is agnostic to the choice
%   of score, so this is legitimate; it is stated explicitly because a
%   reviewer will ask.
%
%   See also SFUQ.UQ.NORMALIZEDCONFORMAL, SFUQ.UQ.QUANTILEHIGHER.

    residCal = abs(residCal(:));
    n = numel(residCal);

    validateattributes(alpha, {'numeric'}, {'scalar', '>', 0, '<', 1});

    % The (n+1) is what makes the guarantee finite-sample rather than
    % asymptotic: it accounts for the unseen test point joining the
    % exchangeable set.
    level = min(ceil((n + 1) * (1 - alpha)) / n, 1);

    if (n + 1) * (1 - alpha) > n
        warning('sfuq:splitConformal:tooFewPoints', ...
                ['Calibration set of %d points is too small for %.1f%% ' ...
                 'coverage; need at least %d. The radius falls back to the ' ...
                 'maximum residual and the guarantee is vacuous.'], ...
                n, 100 * (1 - alpha), ceil(1 / alpha) - 1);
    end

    R = sfuq.uq.quantileHigher(residCal, level);
end
