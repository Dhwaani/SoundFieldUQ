function R = exp06_weighted_conformal(nTrials)
%EXP06_WEIGHTED_CONFORMAL  The project's main result.
%
%   Regime A: the calibration microphones and the listener positions occupy the
%   SAME region but with DIFFERENT densities. Microphones cluster near the
%   centre of the array, as they do in practice; the positions to be certified
%   are spread uniformly over the listening area.
%
%   Because both densities are designed by the experimenter, the likelihood
%   ratio dQ/dP is available in closed form -- no estimation, no estimation
%   error. This is the property that sound-field reconstruction has and almost
%   every other application of weighted conformal prediction lacks.
%
%   The claim under test
%   --------------------
%   Plain split conformal, whose guarantee assumes exchangeability, should
%   UNDER-COVER here: its calibration residuals are drawn disproportionately
%   from the easy centre, so its single constant radius is too small for the
%   harder edges. Weighted conformal with the exact geometric ratio should
%   restore coverage to nominal.
%
%   The sweep over cluster tightness turns this into a falsifiable curve rather
%   than a single number: as the calibration density diverges from the query
%   density, split conformal should degrade monotonically while weighted
%   conformal holds at nominal, and the price should show up in the effective
%   sample size and the interval width.
%
%   Produces figures 7 and 8.

    if nargin < 1 || isempty(nTrials), nTrials = 10; end

    S     = sfuq.data.makeScenario();
    nMic  = 120;
    alpha = 0.10;

    % Cluster tightness. Small sigma means the calibration microphones huddle
    % at the array centre; large sigma approaches the uniform query density.
    sigmas = [0.18 0.25 0.35 0.50 0.80];

    Q = sfuq.geom.uniformBox(S.Centre, S.HalfWidth);

    fprintf('\n=== Experiment 06: weighted conformal under designed density shift ===\n');
    fprintf('  query density: uniform over the %.1f m listening cube\n', 2*S.HalfWidth);
    fprintf('  %8s %8s %10s %10s %10s %11s %11s\n', ...
            'sigma_p', 'ESS/n', 'CP cov', 'nCP cov', 'wCP cov', 'CP width', 'wCP width');

    n = numel(sigmas);
    [covCP, covNC, covWC] = deal(zeros(n, nTrials));
    [wCP,   wWC,   essF ] = deal(zeros(n, nTrials));
    [isCP,  isWC        ] = deal(zeros(n, nTrials));

    for i = 1:n
        P = sfuq.geom.gaussianBox(S.Centre, S.HalfWidth, sigmas(i));

        for t = 1:nTrials
            T = sfuq.shiftTrial(S, nMic, P, Q, ...
                    struct('seed', 6000 + 100*i + t));

            rCP = sfuq.uq.splitConformal(T.resCal, alpha);
            rNC = sfuq.uq.normalizedConformal(T.resCal, T.sigCal, T.sigQuery, alpha);
            [rWC, info] = sfuq.uq.weightedConformal(T.resCal, T.wCal, T.wQuery, alpha);

            covCP(i, t) = sfuq.eval.coverage(T.resQuery, rCP);
            covNC(i, t) = sfuq.eval.coverage(T.resQuery, rNC);
            covWC(i, t) = sfuq.eval.coverage(T.resQuery, rWC);

            wCP(i, t) = rCP;
            wWC(i, t) = mean(rWC(isfinite(rWC)));
            essF(i, t) = info.essFraction;

            isCP(i, t) = sfuq.eval.intervalScore(T.resQuery, rCP, alpha);
            isWC(i, t) = sfuq.eval.intervalScore(T.resQuery, ...
                            rWC(isfinite(rWC)), alpha);
        end

        fprintf('  %8.2f %8.2f %10.3f %10.3f %10.3f %11.5f %11.5f\n', ...
                sigmas(i), mean(essF(i,:)), mean(covCP(i,:)), ...
                mean(covNC(i,:)), mean(covWC(i,:)), ...
                mean(wCP(i,:)), mean(wWC(i,:)));
    end

    R = struct();
    R.sigmas  = sigmas(:);
    R.alpha   = alpha;
    R.nMic    = nMic;
    R.covCP   = mean(covCP, 2);
    R.covNCP  = mean(covNC, 2);
    R.covWCP  = mean(covWC, 2);
    R.sdCP    = std(covCP, 0, 2);
    R.sdWCP   = std(covWC, 0, 2);
    R.widthCP = mean(wCP, 2);
    R.widthWCP= mean(wWC, 2);
    R.essFrac = mean(essF, 2);
    R.scoreCP = mean(isCP, 2);
    R.scoreWCP= mean(isWC, 2);
    R.scenario = S;

    fprintf('\n  Worst-case split-conformal coverage : %.3f (nominal %.2f)\n', ...
            min(R.covCP), 1 - alpha);
    fprintf('  Worst-case weighted coverage        : %.3f\n', min(R.covWCP));
    fprintf('  Max |weighted coverage - nominal|   : %.3f\n', ...
            max(abs(R.covWCP - (1 - alpha))));

    sfuq.viz.figWeightedCoverage(R, 'fig07_weighted_coverage');
    sfuq.viz.figEfficiencyCost(R, 'fig08_efficiency_cost');
    save(fullfile(sfuq.viz.projectRoot(), 'results', 'exp06_weighted.mat'), '-struct', 'R');
end
