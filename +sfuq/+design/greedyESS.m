function [sel, history] = greedyESS(candidates, Qquery, nSelect, opts)
%GREEDYESS  Choose calibration microphone positions by a distribution-free criterion.
%
%   [SEL, HISTORY] = GREEDYESS(CANDIDATES, QQUERY, NSELECT) greedily selects
%   NSELECT rows of CANDIDATES to serve as calibration microphone positions,
%   maximising the effective sample size of the resulting conformal weights
%   over the query density QQUERY.
%
%   How this differs from existing microphone placement
%   ---------------------------------------------------
%   Established sensor placement for sound-field estimation optimises a
%   model-internal quantity: Gaussian-process posterior variance, mutual
%   information, or conditional entropy. Those criteria answer "where do I
%   learn the most about the field?" and their guarantees hold only if the
%   model is correct.
%
%   This criterion answers a different and more operational question: "where do
%   I put microphones so that the reconstruction can be CERTIFIED, distribution
%   free, everywhere the listeners actually are?" It optimises the statistical
%   efficiency of the certificate, not the fit of the model, and it inherits no
%   assumption about the model being right.
%
%   The objective, and why it is the right one
%   ------------------------------------------
%   The conformal weights are w = q/p, and Kish effective sample size is
%   maximised when p = q. So maximising ESS drives the empirical calibration
%   layout toward the query distribution. The greedy search is therefore a
%   discrete approximation to a result that is exactly true in the continuum:
%   THE OPTIMAL CALIBRATION DENSITY IS THE QUERY DENSITY ITSELF. Placement
%   should follow where the listeners are, not where the array is convenient.
%
%   That contradicts common practice, in which calibration measurements are
%   taken wherever the array already sits, and it is a testable, practically
%   consequential claim -- see EXP08_PLACEMENT.
%
%   Options (struct)
%     .lambda       weight on the certifiable-fraction term (default 1.0)
%     .nQuerySamp   query points used to score a design      (default 2000)
%     .sigma        kernel width for the empirical density   (default auto)
%     .stream       RandStream for reproducibility
%
%   Outputs
%     SEL      NSELECT-by-1 indices into CANDIDATES
%     HISTORY  struct array, per step: .ess, .essFraction, .certFraction, .score

    if nargin < 4, opts = struct(); end
    if ~isfield(opts, 'lambda'),     opts.lambda = 1.0;      end
    if ~isfield(opts, 'nQuerySamp'), opts.nQuerySamp = 2000; end
    if ~isfield(opts, 'stream'),     opts.stream = [];       end

    Xq = Qquery.sample(opts.nQuerySamp, opts.stream);
    qAtQ = Qquery.pdf(Xq);

    nc = size(candidates, 1);
    if ~isfield(opts, 'sigma') || isempty(opts.sigma)
        % Bandwidth for the empirical calibration density: a fraction of the
        % candidate cloud's extent. Wide enough that a handful of microphones
        % still define a smooth density, narrow enough to distinguish layouts.
        spread = max(candidates, [], 1) - min(candidates, [], 1);
        opts.sigma = 0.25 * mean(spread(spread > 0));
    end

    sel = zeros(nSelect, 1);
    chosen = false(nc, 1);
    history = repmat(struct('ess', 0, 'essFraction', 0, ...
                            'certFraction', 0, 'score', -Inf), nSelect, 1);

    for step = 1:nSelect
        bestScore = -Inf; bestIdx = 0; bestStats = [];

        for c = find(~chosen).'
            trial = candidates([sel(1:step-1); c], :);
            [score, st] = scoreDesign(trial, Xq, qAtQ, opts);
            if score > bestScore
                bestScore = score; bestIdx = c; bestStats = st;
            end
        end

        sel(step) = bestIdx;
        chosen(bestIdx) = true;
        history(step) = bestStats;
        history(step).score = bestScore;
    end
end

function [score, st] = scoreDesign(P, Xq, qAtQ, opts)
%   Empirical calibration density from the chosen positions (a Gaussian kernel
%   estimate), then the induced weights at the query points.

    D2 = sfuq.kernels.pairwiseDist(Xq, P) .^ 2;
    pAtQ = mean(exp(-0.5 * D2 / opts.sigma^2), 2);

    live = pAtQ > 1e-12;
    w = inf(size(pAtQ));
    w(live) = qAtQ(live) ./ pAtQ(live);

    certFrac = mean(live);

    if any(live)
        [ess, essFrac] = sfuq.uq.effectiveSampleSize(w(live));
    else
        ess = 0; essFrac = 0;
    end

    % Certifiability first, efficiency second: a design that cannot certify a
    % listener position is worse than one that certifies it inefficiently.
    score = certFrac + opts.lambda * essFrac;

    st = struct('ess', ess, 'essFraction', essFrac, ...
                'certFraction', certFrac, 'score', score);
end
