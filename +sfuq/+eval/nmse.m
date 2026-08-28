function v = nmse(est, truth)
%NMSE  Normalised mean square error in decibels.
%
%   V = NMSE(EST, TRUTH) returns 10*log10( sum|est-truth|^2 / sum|truth|^2 ).
%
%   This is the standard reporting metric in the sound-field reconstruction
%   literature, which is why it is used here unchanged: comparability with
%   published numbers is worth more than a metric of one's own devising.
%
%   Note that it is computed against the NOISELESS ground truth, which is
%   available because the field is generated analytically. Evaluating against
%   noisy observations would confound estimator error with measurement noise
%   and flatter every method equally.

    est   = est(:);
    truth = truth(:);
    if numel(est) ~= numel(truth)
        error('sfuq:nmse:sizeMismatch', 'Inputs must have equal length.');
    end

    num = sum(abs(est - truth) .^ 2);
    den = sum(abs(truth) .^ 2);
    if den == 0
        error('sfuq:nmse:zeroReference', 'Reference field has zero energy.');
    end

    v = 10 * log10(num / den);
end
