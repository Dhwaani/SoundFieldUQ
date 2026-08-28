function R = exp08_placement(nTrials)
%EXP08_PLACEMENT  Microphone placement driven by a distribution-free criterion.
%
%   The corollary of experiments 06 and 07, and the most directly useful output
%   of the project.
%
%   The conformal weights are w = q/p and the Kish effective sample size is
%   maximised when p = q. Placement that maximises ESS therefore drives the
%   calibration layout toward the query distribution. In the continuum the
%   optimum is exact: THE BEST CALIBRATION LAYOUT IS A SAMPLE FROM THE QUERY
%   DENSITY ITSELF.
%
%   Why this is worth stating
%   -------------------------
%   Established placement criteria for sound-field estimation -- posterior
%   variance, mutual information, conditional entropy -- optimise how much the
%   model LEARNS. They are the right objective for fitting the field. They are
%   the wrong objective for certifying it, and the two optima differ: a compact
%   array can fit a field perfectly well while being nearly useless for
%   certifying positions away from it.
%
%   This experiment measures that gap directly. Three layouts, identical
%   microphone counts, identical reconstruction model:
%
%     compact   calibration clustered at the array centre (common practice)
%     uniform   calibration spread evenly over the region
%     matched   calibration drawn from the query density (the theoretical optimum)
%
%   The prediction is that all three give correctly calibrated weighted
%   intervals -- validity is not the differentiator -- but that they differ
%   sharply in effective sample size, interval width and certifiable fraction.
%   Compact layouts pay for their convenience in width, not in coverage, which
%   is precisely why the cost is easy to miss without measuring it.
%
%   Produces figure 10.

    if nargin < 1 || isempty(nTrials), nTrials = 10; end

    S     = sfuq.data.makeScenario();
    nMic  = 120;
    nCal  = 400;
    alpha = 0.10;

    Q = sfuq.geom.uniformBox(S.Centre, S.HalfWidth);

    designs = struct( ...
        'name',    {'compact', 'uniform', 'matched'}, ...
        'density', {sfuq.geom.gaussianBox(S.Centre, S.HalfWidth, 0.18), ...
                    sfuq.geom.uniformBox(S.Centre, S.HalfWidth), ...
                    sfuq.geom.uniformBox(S.Centre, S.HalfWidth)});

    fprintf('\n=== Experiment 08: placement by conformal efficiency ===\n');
    fprintf('  %10s %8s %10s %10s %11s %11s\n', ...
            'design', 'ESS/n', 'wCP cov', 'CP cov', 'wCP width', 'certifiable');

    nd = numel(designs);
    [essF, covWC, covCP, widWC, certF] = deal(zeros(nd, nTrials));

    for d = 1:nd
        for t = 1:nTrials
            T = sfuq.shiftTrial(S, nMic, designs(d).density, Q, ...
                    struct('seed', 8000 + 100*d + t, 'nCal', nCal));

            rCP = sfuq.uq.splitConformal(T.resCal, alpha);
            [rWC, info] = sfuq.uq.weightedConformal(T.resCal, T.wCal, T.wQuery, alpha);

            essF(d, t)  = info.essFraction;
            covCP(d, t) = sfuq.eval.coverage(T.resQuery, rCP);
            covWC(d, t) = sfuq.eval.coverage(T.resQuery, rWC);
            widWC(d, t) = mean(rWC(isfinite(rWC)));
            certF(d, t) = mean(T.certifiable);
        end

        fprintf('  %10s %8.2f %10.3f %10.3f %11.5f %11.3f\n', ...
                designs(d).name, mean(essF(d,:)), mean(covWC(d,:)), ...
                mean(covCP(d,:)), mean(widWC(d,:)), mean(certF(d,:)));
    end

    % --- greedy design over a candidate grid ---------------------------
    fprintf('\n  Greedy ESS placement over a candidate grid...\n');
    g = linspace(-S.HalfWidth, S.HalfWidth, 7);
    [GX, GY, GZ] = ndgrid(g, g, g);
    candidates = S.Centre + [GX(:), GY(:), GZ(:)];

    sel = sfuq.design.greedyESS(candidates, Q, 24, ...
              struct('nQuerySamp', 800, 'stream', ...
                     RandStream('mt19937ar', 'Seed', 8999)));
    greedyPts = candidates(sel, :);

    spread = @(X) mean(sqrt(sum((X - mean(X, 1)).^2, 2)));
    fprintf('  greedy layout mean radius %.3f m vs region half-width %.3f m\n', ...
            spread(greedyPts), S.HalfWidth);

    R = struct();
    R.designNames = {designs.name};
    R.essFrac     = mean(essF, 2);
    R.covWCP      = mean(covWC, 2);
    R.covCP       = mean(covCP, 2);
    R.widthWCP    = mean(widWC, 2);
    R.certFrac    = mean(certF, 2);
    R.greedyPts   = greedyPts;
    R.candidates  = candidates;
    R.alpha       = alpha;
    R.scenario    = S;

    fprintf('\n  Compact vs matched: coverage %.3f vs %.3f (both valid),\n', ...
            R.covWCP(1), R.covWCP(3));
    fprintf('  but effective sample size %.2f vs %.2f and width %.5f vs %.5f.\n', ...
            R.essFrac(1), R.essFrac(3), R.widthWCP(1), R.widthWCP(3));
    fprintf('  The cost of a convenient layout is paid in width, not coverage.\n');

    sfuq.viz.figPlacement(R, 'fig10_placement');
    save(fullfile(sfuq.viz.projectRoot(), 'results', 'exp08_placement.mat'), '-struct', 'R');
end
