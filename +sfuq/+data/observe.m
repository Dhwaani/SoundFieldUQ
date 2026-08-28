function [y, pTrue] = observe(X, S, rngStream)
%OBSERVE  Simulated microphone measurements with additive complex noise.
%
%   [Y, PTRUE] = OBSERVE(X, S, STREAM) returns noisy observations Y at the
%   positions X for scenario S, together with the noiseless ground truth
%   PTRUE. Noise is circularly-symmetric complex Gaussian with standard
%   deviation S.NoiseRel times the RMS field magnitude over X.
%
%   Scaling the noise to the field level rather than fixing it absolutely
%   keeps the effective SNR constant as the source, frequency or region move,
%   so a coverage change across conditions reflects the estimator rather than
%   an accidental change in SNR.

    pTrue = sfuq.data.fieldAt(X, S.imgPos, S.imgAmp, S.k);
    lvl   = sqrt(mean(abs(pTrue) .^ 2));

    n = numel(pTrue);
    if nargin < 3 || isempty(rngStream)
        w = randn(n, 1) + 1i * randn(n, 1);
    else
        w = randn(rngStream, n, 1) + 1i * randn(rngStream, n, 1);
    end

    y = pTrue + S.NoiseRel * lvl * w / sqrt(2);
end
