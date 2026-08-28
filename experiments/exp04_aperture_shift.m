function R = exp04_aperture_shift(nTrials)
%EXP04_APERTURE_SHIFT  What happens when exchangeability is violated.
%
%   Split conformal's guarantee is conditional on the calibration and test
%   points being exchangeable. In practice they very often are not: an array
%   is calibrated with the microphones it has, and then the reconstruction is
%   queried at listener positions outside that aperture. This experiment makes
%   that violation explicit and measures the damage.
%
%   Calibration points stay INSIDE the microphone aperture. Test points are
%   drawn from a region progressively larger than it. The x-axis is therefore
%   a controlled dial on the severity of the covariate shift.
%
%   The expected finding
%   --------------------
%   Split conformal degrades sharply and without warning: its radius is a
%   single number fixed by the calibration set, and it has no mechanism to
%   notice that the test points are harder. Coverage can fall to a small
%   fraction of the nominal level while the method continues to report the
%   same confident radius. The GP variance and the normalized conformal radius
%   both grow with distance from the data and degrade far more gracefully.
%
%   Honesty note for the paper: normalization does NOT restore the guarantee.
%   Nothing does, under covariate shift, short of the weighted-conformal
%   machinery that needs a likelihood ratio nobody has here. The correct claim
%   is that it fails gracefully rather than catastrophically, and the paper
%   should say exactly that.
%
%   Produces figure 5.

    if nargin < 1 || isempty(nTrials)
        nTrials = 10;
    end

    S = sfuq.data.makeScenario();
    nMic  = 120;
    alpha = 0.10;
    ratios = [1.0 1.25 1.5 1.75 2.0 2.5];

    fprintf('\n=== Experiment 04: aperture extrapolation ===\n');
    fprintf('  calibration inside the aperture, test region scaled outward\n');
    fprintf('  %8s %10s %10s %10s\n', 'ratio', 'GP cov', 'CP cov', 'nCP cov');

    covGP = zeros(numel(ratios), nTrials);
    covCP = zeros(numel(ratios), nTrials);
    covNC = zeros(numel(ratios), nTrials);
    nmse  = zeros(numel(ratios), nTrials);

    for i = 1:numel(ratios)
        for t = 1:nTrials
            T = sfuq.trial(S, nMic, struct( ...
                'seed', 4000 + 100 * i + t, ...
                'nTest', 1200, ...
                'testHalf', ratios(i) * S.HalfWidth));

            rGP = sfuq.uq.gpRadius(T.sigTest .^ 2, alpha);
            rCP = sfuq.uq.splitConformal(T.resCal, alpha);
            rNC = sfuq.uq.normalizedConformal(T.resCal, T.sigCal, T.sigTest, alpha);

            covGP(i, t) = sfuq.eval.coverage(T.resTest, rGP);
            covCP(i, t) = sfuq.eval.coverage(T.resTest, rCP);
            covNC(i, t) = sfuq.eval.coverage(T.resTest, rNC);
            nmse(i, t)  = T.nmse;
        end
        fprintf('  %8.2f %10.3f %10.3f %10.3f\n', ratios(i), ...
                mean(covGP(i, :)), mean(covCP(i, :)), mean(covNC(i, :)));
    end

    R = struct('ratios', ratios(:), 'alpha', alpha, ...
               'covGP', mean(covGP, 2), 'covCP', mean(covCP, 2), ...
               'covNCP', mean(covNC, 2), 'nmse', mean(nmse, 2), ...
               'scenario', S, 'nMic', nMic);

    c = sfuq.viz.colors();
    fig = figure('Units', 'inches', 'Position', [1 1 3.5 2.6], 'Color', 'w');
    ax = axes(fig); hold(ax, 'on');
    yline(ax, 1 - alpha, '--', 'Color', c.ideal, 'LineWidth', 1.0, ...
          'DisplayName', 'nominal');
    plot(ax, R.ratios, R.covGP,  '-o', 'Color', c.gp, 'LineWidth', 1.2, ...
         'MarkerFaceColor', c.gp, 'MarkerSize', 4, 'DisplayName', 'GP posterior');
    plot(ax, R.ratios, R.covCP,  '-s', 'Color', c.conformal, 'LineWidth', 1.2, ...
         'MarkerFaceColor', c.conformal, 'MarkerSize', 4, ...
         'DisplayName', 'split conformal');
    plot(ax, R.ratios, R.covNCP, '-^', 'Color', c.normalized, 'LineWidth', 1.2, ...
         'MarkerFaceColor', c.normalized, 'MarkerSize', 4, ...
         'DisplayName', 'normalized conformal');
    xlabel(ax, 'Test region size / aperture size');
    ylabel(ax, 'Empirical coverage');
    ylim(ax, [0 1.02]);
    legend(ax, 'Location', 'southwest');
    sfuq.viz.style(ax);
    sfuq.viz.saveFig(fig, 'fig05_aperture_shift');

    save(fullfile(sfuq.viz.projectRoot(), 'results', 'exp04_shift.mat'), '-struct', 'R');
end
