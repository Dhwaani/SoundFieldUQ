function [mu, variance] = predictGP(M, Xq)
%PREDICTGP  Posterior mean and variance of a fitted sound-field GP.
%
%   [MU, VAR] = PREDICTGP(M, XQ) returns the posterior mean MU (complex) and
%   the posterior variance VAR (real, non-negative) at query points XQ.
%
%   VAR is the variance of EACH of the real and imaginary parts, which is the
%   convention SFUQ.UQ.GPRADIUS assumes when converting it into a predictive
%   disk radius. Stating that convention explicitly matters: getting the
%   factor of two wrong here shifts every coverage number by several points
%   and is invisible unless the calibration is checked against a known case,
%   which TCONFORMAL does.
%
%   The variance is computed as s - ||L'\ks||^2 rather than by forming an
%   explicit inverse; the triangular solve is both faster and better
%   conditioned.

    Ks = M.kernelFcn(Xq, M.X);          % Nq-by-n
    mu = Ks * M.alpha;

    if nargout > 1
        V = M.L' \ Ks.';                % n-by-Nq
        variance = M.s * (1 - sum(V .^ 2, 1).');

        % The subtraction can go slightly negative for query points that
        % coincide with training points. Clamp rather than let a negative
        % variance propagate into a complex radius.
        variance = max(variance, eps);
    end
end
