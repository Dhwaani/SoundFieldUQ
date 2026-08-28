function [nominal, empirical] = reliabilityCurve(residCal, sigmaCal, ...
                                                 residTest, sigmaTest, ...
                                                 method, levels)
%RELIABILITYCURVE  Empirical coverage across a sweep of nominal levels.
%
%   [NOM, EMP] = RELIABILITYCURVE(RESIDCAL, SIGMACAL, RESIDTEST, SIGMATEST,
%   METHOD, LEVELS) recomputes the prediction radii at every nominal coverage
%   level in LEVELS and measures the coverage actually achieved.
%
%   METHOD is one of 'gp', 'conformal' or 'normalized'.
%
%   The resulting curve is the diagnostic that makes calibration visible: a
%   perfectly calibrated method lies on the diagonal. Reporting a single
%   coverage number at 90% hides whether a method is calibrated everywhere or
%   merely happens to cross the diagonal at one point, and that distinction
%   changes the conclusion.

    if nargin < 6 || isempty(levels)
        levels = 0.50:0.05:0.99;
    end

    nominal   = levels(:);
    empirical = zeros(size(nominal));

    for i = 1:numel(nominal)
        a = 1 - nominal(i);
        switch lower(method)
            case 'gp'
                R = sfuq.uq.gpRadius(sigmaTest .^ 2, a);
            case 'conformal'
                R = sfuq.uq.splitConformal(residCal, a);
            case 'normalized'
                R = sfuq.uq.normalizedConformal(residCal, sigmaCal, sigmaTest, a);
            otherwise
                error('sfuq:reliabilityCurve:badMethod', ...
                      'Unknown method "%s".', method);
        end
        empirical(i) = sfuq.eval.coverage(residTest, R);
    end
end
