function fig = figWeightedCoverage(R, outName)
%FIGWEIGHTEDCOVERAGE  The project's headline figure.
%
%   Coverage against calibration-cluster tightness. As the calibration density
%   diverges from the query density, split conformal degrades monotonically
%   while weighted conformal -- using the exact geometric likelihood ratio --
%   holds at nominal. The x-axis is a dial on how badly exchangeability is
%   violated, and the separation between the two curves is the contribution.

    if nargin < 2, outName = ''; end
    c = sfuq.viz.colors();

    fig = figure('Units', 'inches', 'Position', [1 1 3.5 2.7], 'Color', 'w');
    ax = axes(fig); hold(ax, 'on');

    target = 1 - R.alpha;
    yline(ax, target, '--', 'Color', c.ideal, 'LineWidth', 1.0, ...
          'DisplayName', 'nominal');

    errorbar(ax, R.sigmas, R.covCP, R.sdCP, '-s', 'Color', c.conformal, ...
             'LineWidth', 1.2, 'MarkerFaceColor', c.conformal, ...
             'MarkerSize', 4, 'CapSize', 3, 'DisplayName', 'split conformal');
    plot(ax, R.sigmas, R.covNCP, '-o', 'Color', c.normalized, ...
         'LineWidth', 1.2, 'MarkerFaceColor', c.normalized, ...
         'MarkerSize', 4, 'DisplayName', 'normalized conformal');
    errorbar(ax, R.sigmas, R.covWCP, R.sdWCP, '-^', 'Color', c.weighted, ...
             'LineWidth', 1.4, 'MarkerFaceColor', c.weighted, ...
             'MarkerSize', 5, 'CapSize', 3, ...
             'DisplayName', 'weighted (exact dQ/dP)');

    xlabel(ax, 'Calibration cluster width \sigma_p (m)');
    ylabel(ax, 'Empirical coverage');
    ylim(ax, [0.4 1.02]);
    title(ax, 'Designed density shift: mics cluster, listeners do not');
    legend(ax, 'Location', 'southeast');
    sfuq.viz.style(ax);

    if ~isempty(outName), sfuq.viz.saveFig(fig, outName); end
end
