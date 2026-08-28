function [ess, frac] = effectiveSampleSize(w)
%EFFECTIVESAMPLESIZE  Kish effective sample size of a weight vector.
%
%   [ESS, FRAC] = EFFECTIVESAMPLESIZE(W) returns
%
%       ESS = (sum w)^2 / sum(w^2),      FRAC = ESS / numel(w).
%
%   ESS equals numel(W) exactly when all weights are equal and falls toward 1
%   as the weights concentrate on a few points.
%
%   Why it is the right diagnostic here
%   -----------------------------------
%   Weighted conformal buys back validity under covariate shift, but it pays in
%   statistical efficiency: the (1-alpha) quantile is effectively estimated from
%   ESS points, not n. When the calibration microphones cluster far from the
%   listener positions being certified, ESS collapses and the certified
%   intervals become wide and unstable, even though coverage remains correct.
%
%   That gives the quantity a second, more useful life: since the weights are
%   w = q/p and ESS is maximised when p = q, MAXIMISING ESS OVER MICROPHONE
%   POSITIONS IS EQUIVALENT TO MATCHING THE CALIBRATION LAYOUT TO THE QUERY
%   DISTRIBUTION. It is therefore not merely a diagnostic but a design
%   objective, and it is the criterion SFUQ.DESIGN.GREEDYESS optimises.
%
%   Infinite weights give ESS = 0: a design that cannot certify part of the
%   query region has no usable effective sample there.

    w = w(:);
    if isempty(w)
        ess = 0; frac = 0; return
    end
    if any(~isfinite(w))
        ess = 0; frac = 0; return
    end

    s = sum(w);
    if s <= 0
        ess = 0; frac = 0; return
    end

    ess  = s^2 / sum(w.^2);
    frac = ess / numel(w);
end
