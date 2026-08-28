function R = gpRadius(variance, alpha)
%GPRADIUS  Predictive disk radius implied by a GP posterior variance.
%
%   R = GPRADIUS(VARIANCE, ALPHA) converts the per-component posterior
%   variance of a complex GP into the radius of the central (1-ALPHA)
%   predictive disk.
%
%   Derivation
%   ----------
%   Let the error e = y - mu have independent real and imaginary parts, each
%   N(0, v). Then |e|^2 / v is chi-squared with two degrees of freedom, whose
%   CDF is 1 - exp(-x/2). Setting that equal to 1 - alpha gives
%   x = -2*log(alpha), hence
%
%       R = sqrt( -2 * v * log(alpha) ).
%
%   Equivalently |e| is Rayleigh distributed with scale sqrt(v).
%
%   This is the GP's own claim about its uncertainty, made under the
%   assumption that its model is correct. The point of the experiments is to
%   check that claim empirically rather than take it on trust, so this
%   function is the "Bayesian baseline" against which the distribution-free
%   methods are measured.
%
%   Getting the factor of two wrong here is the easiest mistake in the whole
%   project and would bias every coverage number; TCONFORMAL pins it with a
%   Monte Carlo test against a known Gaussian.

    validateattributes(alpha, {'numeric'}, {'scalar', '>', 0, '<', 1});
    R = sqrt(-2 * variance(:) * log(alpha));
end
