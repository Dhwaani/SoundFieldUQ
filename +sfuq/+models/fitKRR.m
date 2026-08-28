function M = fitKRR(X, y, k, lambda, kernelFcn)
%FITKRR  Helmholtz kernel ridge regression (point estimate only).
%
%   M = FITKRR(X, Y, K, LAMBDA) fits the regularised interpolant
%
%       min_f  sum_i |f(x_i) - y_i|^2 + lambda * ||f||_H^2
%
%   whose solution is f(x) = k(x,X) * (K + lambda I)^{-1} y.
%
%   This is the classical wave-based baseline (Ueno, Koyama and Saruwatari's
%   kernel interpolation of sound fields). Its posterior mean coincides with
%   the GP mean when lambda equals the noise variance, so it is included
%   mainly to make that equivalence explicit and to serve as the point
%   estimator that the distribution-free wrappers are applied to when no
%   variance estimate is available at all.
%
%   See also SFUQ.MODELS.PREDICTKRR, SFUQ.MODELS.FITGP.

    if nargin < 5 || isempty(kernelFcn)
        kernelFcn = @(a, b) sfuq.kernels.helmholtz(a, b, k);
    end

    y = y(:);
    n = numel(y);
    K = kernelFcn(X, X) + lambda * eye(n);
    K = (K + K.') / 2;

    M = struct('X', X, 'k', k, 'lambda', lambda, 'kernelFcn', kernelFcn, ...
               'alpha', K \ y, 'n', n);
end
