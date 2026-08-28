function fig = figEfficiencyCost(R, outName)
%FIGEFFICIENCYCOST  What validity under shift actually costs.
%
%   Two panels sharing an x-axis. Top: effective sample size as a fraction of
%   the calibration set. Bottom: mean interval width for split and weighted
%   conformal.
%
%   Included so the method cannot be read as a free lunch. Weighted conformal
%   buys back coverage that split conformal loses, and it pays for it in
%   effective sample size and in width. Split conformal's intervals are
%   narrower precisely because they are wrong; comparing widths without
%   comparing coverage would invert the conclusion, which is why the two panels
%   share a figure.

    if nargin < 2, outName = ''; end
    c = sfuq.viz.colors();

    fig = figure('Units', 'inches', 'Position', [1 1 3.5 3.4], 'Color', 'w');
    tl = tiledlayout(fig, 2, 1, 'TileSpacing', 'compact', 'Padding', 'compact');

    ax1 = nexttile(tl); hold(ax1, 'on');
    plot(ax1, R.sigmas, R.essFrac, '-o', 'Color', c.weighted, ...
         'LineWidth', 1.3, 'MarkerFaceColor', c.weighted, 'MarkerSize', 4);
    yline(ax1, 1, ':', 'Color', c.ideal);
    ylabel(ax1, 'ESS / n');
    ylim(ax1, [0 1.05]);
    set(ax1, 'XTickLabel', []);
    title(ax1, 'Price of validity under covariate shift');
    sfuq.viz.style(ax1);

    ax2 = nexttile(tl); hold(ax2, 'on');
    plot(ax2, R.sigmas, R.widthCP, '-s', 'Color', c.conformal, ...
         'LineWidth', 1.2, 'MarkerFaceColor', c.conformal, ...
         'MarkerSize', 4, 'DisplayName', 'split conformal (under-covers)');
    plot(ax2, R.sigmas, R.widthWCP, '-^', 'Color', c.weighted, ...
         'LineWidth', 1.3, 'MarkerFaceColor', c.weighted, ...
         'MarkerSize', 5, 'DisplayName', 'weighted (valid)');
    xlabel(ax2, 'Calibration cluster width \sigma_p (m)');
    ylabel(ax2, 'Mean interval radius');
    legend(ax2, 'Location', 'northeast');
    sfuq.viz.style(ax2);

    if ~isempty(outName), sfuq.viz.saveFig(fig, outName); end
end
