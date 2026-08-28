function q = quantileHigher(x, p)
%QUANTILEHIGHER  Empirical quantile with upward rounding.
%
%   Q = QUANTILEHIGHER(X, P) returns the smallest order statistic of X whose
%   rank is at least P * numel(X).
%
%   Why not use QUANTILE
%   --------------------
%   Two reasons, both important.
%
%   First, MATLAB's QUANTILE lives in the Statistics and Machine Learning
%   Toolbox. Implementing this one function keeps the entire project runnable
%   on base MATLAB.
%
%   Second, and more fundamentally, QUANTILE interpolates between order
%   statistics. The finite-sample coverage guarantee of split conformal
%   prediction requires the ceil((n+1)(1-alpha))-th ORDER STATISTIC, not an
%   interpolated value. An interpolated quantile is slightly smaller, and the
%   guarantee is lost: coverage falls below nominal by O(1/n). That failure is
%   easy to miss because the shortfall is small at large n, so the correct
%   estimator is written out explicitly here and pinned by a unit test.

    x = sort(x(:));
    n = numel(x);
    if n == 0
        error('sfuq:quantileHigher:empty', 'Cannot take a quantile of an empty set.');
    end

    idx = ceil(p * n);
    idx = min(max(idx, 1), n);
    q = x(idx);
end
