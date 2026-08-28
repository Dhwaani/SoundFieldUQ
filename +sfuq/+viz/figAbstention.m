function fig = figAbstention(R, outName)
%FIGABSTENTION  Silent failure versus honest abstention.
%
%   As the region to be certified grows beyond the calibration support, split
%   conformal's coverage collapses while it continues to issue confident finite
%   intervals. Weighted conformal instead declines to certify an increasing
%   fraction of points, and remains correct on the ones it does certify.
%
%   The two curves answer different questions and are plotted together on
%   purpose: the left axis is "how often is the method right", the right axis
%   is "how often does it decline to answer". A method that never declines and
%   is wrong 70% of the time is strictly worse than one that declines 87% of
%   the time and is right on the rest.

    if nargin < 2, outName = ''; end
    c = sfuq.viz.colors();

    fig = figure('Units', 'inches', 'Position', [1 1 3.6 2.7], 'Color', 'w');
    ax = axes(fig); hold(ax, 'on');

    yyaxis(ax, 'left');
    ax.YColor = [0 0 0];
    yline(ax, 1 - R.alpha, '--', 'Color', c.ideal, 'LineWidth', 1.0, ...
          'HandleVisibility', 'off');
    plot(ax, R.ratios, R.covCP, '-s', 'Color', c.conformal, 'LineWidth', 1.2, ...
         'MarkerFaceColor', c.conformal, 'MarkerSize', 4, ...
         'DisplayName', 'split conformal coverage');
    plot(ax, R.ratios, R.covWCP, '-^', 'Color', c.weighted, 'LineWidth', 1.4, ...
         'MarkerFaceColor', c.weighted, 'MarkerSize', 5, ...
         'DisplayName', 'weighted coverage (certified pts)');
    ylabel(ax, 'Empirical coverage');
    ylim(ax, [0 1.05]);

    yyaxis(ax, 'right');
    ax.YColor = c.warn;
    plot(ax, R.ratios, R.abstention, ':d', 'Color', c.warn, 'LineWidth', 1.2, ...
         'MarkerFaceColor', c.warn, 'MarkerSize', 4, ...
         'DisplayName', 'weighted abstention rate');
    ylabel(ax, 'Abstention rate');
    ylim(ax, [0 1.05]);

    xlabel(ax, 'Query region size / calibration support');
    title(ax, 'Beyond the support, abstention beats a confident answer');
    legend(ax, 'Location', 'east');
    sfuq.viz.style(ax);

    if ~isempty(outName), sfuq.viz.saveFig(fig, outName); end
end
