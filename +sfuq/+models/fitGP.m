function M = fitGP(X, y, k, sigmaRel, kernelFcn)
%FITGP  Fit a Gaussian process sound-field model with a Helmholtz prior.
%
%   M = FITGP(X, Y, K, SIGMAREL) fits a GP to complex pressure observations Y
%   at positions X, using the Helmholtz kernel at wavenumber K and a relative
%   noise standard deviation SIGMAREL.
%
%   M = FITGP(X, Y, K, SIGMAREL, KERNELFCN) uses a custom kernel, called as
%   KERNELFCN(A, B). Pass @(a,b) sfuq.kernels.gaussian(a, b, ell) for the
%   physics-free ablation.
%
%   The signal variance is fitted, not guessed
%   ------------------------------------------
%   The covariance is parameterised as
%
%       C = s * ( Kernel(X,X) + sigmaRel^2 * I )
%
%   so that the single scale factor s multiplies both signal and noise. Under
%   that parameterisation the type-II maximum-likelihood estimate of s is
%   available in closed form,
%
%       s_hat = y' * (Kernel + sigmaRel^2 I)^{-1} * y  /  (2 * n),
%
%   with the factor 2 because each complex observation carries two real
%   degrees of freedom. No numerical optimisation is required.
%
%   This matters more than it might appear. An unfitted signal variance is the
%   single largest source of mis-calibrated GP intervals: leaving s at 1 while
%   the field magnitude is order 0.01 inflates every predictive interval by
%   more than an order of magnitude, and produces the misleading impression
%   that GP uncertainty is useless. Fitting s is what gives the GP a fair
%   hearing in the comparison against conformal methods.
%
%   Output struct M
%     .X, .k, .sigmaRel, .kernelFcn   model definition
%     .L          Cholesky factor (upper triangular, C0 = L'*L)
%     .alpha      weights for the posterior mean
%     .s          fitted signal variance
%     .n          number of training points
%
%   See also SFUQ.MODELS.PREDICTGP.

    if nargin < 5 || isempty(kernelFcn)
        kernelFcn = @(a, b) sfuq.kernels.helmholtz(a, b, k);
    end

    y = y(:);
    n = numel(y);
    if size(X, 1) ~= n
        error('sfuq:fitGP:sizeMismatch', ...
              'X has %d rows but y has %d elements.', size(X, 1), n);
    end

    C0 = kernelFcn(X, X) + sigmaRel ^ 2 * eye(n);

    % Symmetrise before factorising: the kernel is mathematically symmetric
    % but the distance expansion introduces asymmetry at the 1e-16 level,
    % which CHOL rejects.
    C0 = (C0 + C0.') / 2;

    [L, flag] = chol(C0);
    if flag ~= 0
        % Jitter escalation. Reaching the largest jitter means the microphone
        % layout is effectively rank deficient at this wavenumber, which is a
        % real modelling failure, not a numerical hiccup, so it raises.
        jitter = 1e-12 * trace(C0) / n;
        for attempt = 1:8
            [L, flag] = chol(C0 + jitter * eye(n));
            if flag == 0
                break
            end
            jitter = jitter * 10;
        end
        if flag ~= 0
            error('sfuq:fitGP:notPositiveDefinite', ...
                  ['Kernel matrix is not positive definite even with ' ...
                   'jitter %.2e. Check for duplicate microphone positions.'], ...
                  jitter);
        end
    end

    alpha = L \ (L' \ y);

    % Closed-form type-II ML scale. real() guards against a residual
    % imaginary part of order 1e-18 from the complex quadratic form.
    s = real(y' * alpha) / (2 * n);

    M = struct('X', X, 'k', k, 'sigmaRel', sigmaRel, 'kernelFcn', kernelFcn, ...
               'L', L, 'alpha', alpha, 's', s, 'n', n);
end
