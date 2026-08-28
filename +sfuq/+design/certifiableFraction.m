function [frac, cert, X] = certifiableFraction(Pcal, Qquery, nSamples, stream)
%CERTIFIABLEFRACTION  What proportion of the listening area can be certified.
%
%   [FRAC, CERT, X] = CERTIFIABLEFRACTION(PCAL, QQUERY, N) draws N points from
%   the query density and reports the fraction at which the likelihood ratio is
%   finite -- that is, the fraction of listener positions for which weighted
%   conformal prediction can issue a finite, guaranteed interval.
%
%   This is the deployable output of the whole project.
%
%   A practitioner planning a measurement does not primarily want a coverage
%   number; they want to know which seats in the room their microphone layout
%   will be able to certify, BEFORE they measure anything. Because the
%   likelihood ratio comes from geometry alone, this can be answered with no
%   measurements, no field data and no hardware -- only the intended microphone
%   layout and the intended listening area.
%
%   FRAC < 1 is not a defect of the method. It is the honest statement that
%   part of the listening area lies outside what the microphones can support,
%   and it is exactly the quantity that SFUQ.DESIGN.GREEDYESS is designed to
%   drive to one.
%
%   Outputs
%     FRAC  scalar in [0, 1]
%     CERT  N-by-1 logical, per sampled query point
%     X     N-by-3 the sampled query points, for plotting the certifiable map

    if nargin < 3 || isempty(nSamples), nSamples = 5000; end
    if nargin < 4, stream = []; end

    X = Qquery.sample(nSamples, stream);
    [~, cert] = sfuq.geom.likelihoodRatio(Pcal, Qquery, X);
    frac = mean(cert);
end
