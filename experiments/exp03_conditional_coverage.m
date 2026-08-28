function R = exp03_conditional_coverage(nTrials)
%EXP03_CONDITIONAL_COVERAGE  Where in space is each method actually valid?
%
%   This is the project's central experiment.
%
%   Marginal coverage is an average over test positions, and an average can be
%   exactly right while the underlying distribution is badly wrong. A method
%   that covers 100% of points near the microphones and 75% of points at the
%   edge of the aperture reports a flawless 90% overall. For anyone deciding
%   whether to trust a reconstruction at a PARTICULAR location, the marginal
%   number is close to meaningless.
%
%   Here coverage is binned by distance from the centre of the measurement
%   region, which is a proxy for how far a query point sits from the
%   supporting microphones.
%
%   The expected finding
%   --------------------
%   Split conformal, which assigns the same radius everywhere, is the WORST of
%   the three conditionally, despite being the only one with a formal
%   guarantee. It over-covers near the centre and under-covers at the edges.
%   The GP posterior, whose variance grows away from the data, is far more
%   uniform. Normalized conformal keeps the guarantee and inherits the shape.
%
%   That inversion, a guaranteed method behaving worse than an unguaranteed
%   one on the metric a practitioner cares about, is the result worth
%   publishing. It is also the standard theoretical position on conditional
%   validity, made concrete and quantitative for a real acoustics problem.
%
%   Produces figures 3 and 4.

    if nargin < 1 || isempty(nTrials)
        nTrials = 12;
    end

    S = sfuq.data.makeScenario();
    nMic  = 120;
    alpha = 0.10;
    edges = [0, 0.25, 0.40, 0.55, 0.70, 0.90];

    fprintf('\n=== Experiment 03: conditional (spatial) coverage ===\n');
    fprintf('  target coverage %.0f%%, %d microphones, %d trials\n', ...
            100 * (1 - alpha), nMic, nTrials);

    nb = numel(edges) - 1;
    covGP = nan(nb, nTrials);
    covCP = nan(nb, nTrials);
    covNC = nan(nb, nTrials);
    nBin  = zeros(nb, nTrials);

    margGP = zeros(nTrials, 1);
    margCP = zeros(nTrials, 1);
    margNC = zeros(nTrials, 1);
    isGP = zeros(nTrials, 1); isCP = zeros(nTrials, 1); isNC = zeros(nTrials, 1);

    for t = 1:nTrials
        T = sfuq.trial(S, nMic, struct('seed', 3000 + t, 'nTest', 1500));

        rGP = sfuq.uq.gpRadius(T.sigTest .^ 2, alpha);
        rCP = repmat(sfuq.uq.splitConformal(T.resCal, alpha), numel(T.resTest), 1);
        rNC = sfuq.uq.normalizedConformal(T.resCal, T.sigCal, T.sigTest, alpha);

        margGP(t) = sfuq.eval.coverage(T.resTest, rGP);
        margCP(t) = sfuq.eval.coverage(T.resTest, rCP);
        margNC(t) = sfuq.eval.coverage(T.resTest, rNC);

        isGP(t) = sfuq.eval.intervalScore(T.resTest, rGP, alpha);
        isCP(t) = sfuq.eval.intervalScore(T.resTest, rCP, alpha);
        isNC(t) = sfuq.eval.intervalScore(T.resTest, rNC, alpha);

        cG = sfuq.eval.conditionalCoverage(T.resTest, rGP, T.dTest, edges);
        cC = sfuq.eval.conditionalCoverage(T.resTest, rCP, T.dTest, edges);
        cN = sfuq.eval.conditionalCoverage(T.resTest, rNC, T.dTest, edges);

        covGP(:, t) = [cG.coverage];
        covCP(:, t) = [cC.coverage];
        covNC(:, t) = [cN.coverage];
        nBin(:, t)  = [cG.n];
    end

    R = struct();
    R.edges = edges;
    R.alpha = alpha;
    R.nPerBin = mean(nBin(:));
    R.apertureHalfWidth = S.HalfWidth;
    R.covGP  = mean(covGP,  2, 'omitnan');
    R.covCP  = mean(covCP,  2, 'omitnan');
    R.covNCP = mean(covNC,  2, 'omitnan');
    R.marginal = struct('GP', mean(margGP), 'CP', mean(margCP), 'NCP', mean(margNC));
    R.intervalScore = struct('GP', mean(isGP), 'CP', mean(isCP), 'NCP', mean(isNC));
    R.scenario = S;

    fprintf('  marginal coverage : GP %.3f | CP %.3f | nCP %.3f\n', ...
            R.marginal.GP, R.marginal.CP, R.marginal.NCP);
    fprintf('  interval score    : GP %.4f | CP %.4f | nCP %.4f  (lower better)\n', ...
            R.intervalScore.GP, R.intervalScore.CP, R.intervalScore.NCP);
    fprintf('  worst-bin coverage: GP %.3f | CP %.3f | nCP %.3f\n', ...
            min(R.covGP), min(R.covCP), min(R.covNCP));
    % max - min written out rather than RANGE(), which needs the Statistics
    % Toolbox. This spread is the headline number of the experiment: it is how
    % far coverage varies across space, and split conformal should be worst.
    spread = @(v) max(v) - min(v);
    fprintf('  spread over bins  : GP %.3f | CP %.3f | nCP %.3f  (max - min)\n', ...
            spread(R.covGP), spread(R.covCP), spread(R.covNCP));

    % --- spatial radius maps for figure 4 -------------------------------
    T = sfuq.trial(S, nMic, struct('seed', 3999, 'nTest', 10));
    g = linspace(-1.2 * S.HalfWidth, 1.2 * S.HalfWidth, 90);
    [GX, GY] = meshgrid(S.Centre(1) + g, S.Centre(2) + g);
    grid3 = [GX(:), GY(:), repmat(S.Centre(3), numel(GX), 1)];

    [~, vGrid] = sfuq.models.predictGP(T.model, grid3);
    rCPscalar = sfuq.uq.splitConformal(T.resCal, alpha);
    rNCgrid = sfuq.uq.normalizedConformal(T.resCal, T.sigCal, sqrt(vGrid), alpha);

    R.mapX = S.Centre(1) + g;
    R.mapY = S.Centre(2) + g;
    R.mapCP  = repmat(rCPscalar, size(GX));
    R.mapNCP = reshape(rNCgrid, size(GX));
    R.micXY  = T.Xmic(:, 1:2);
    R.apertureBox = [S.Centre(1) - S.HalfWidth, S.Centre(1) + S.HalfWidth, ...
                     S.Centre(2) - S.HalfWidth, S.Centre(2) + S.HalfWidth];

    sfuq.viz.figConditionalCoverage(R, 'fig03_conditional_coverage');
    sfuq.viz.figSpatialMap(R, 'fig04_spatial_radius');
    save(fullfile(sfuq.viz.projectRoot(), 'results', 'exp03_conditional.mat'), '-struct', 'R');
end
