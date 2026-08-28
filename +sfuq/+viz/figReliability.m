function fig = figReliability(R, outName)
%FIGRELIABILITY  Reliability diagram: nominal versus empirical coverage.
%
%   The headline calibration figure. A perfectly calibrated method lies on the
%   diagonal; a method below the diagonal is over-confident (its intervals are
%   too narrow and it lies about how often it is right), and one above is
%   over-cautious.
%
%   Input R is the struct produced by EXP02_MARGINAL_CALIBRATION.

    c = sfuq.viz.colors();
    fig = figure('Units', 'inches', 'Position', [1 1 3.5 2.9], 'Color', 'w');
    ax = axes(fig); hold(ax, 'on');

    plot(ax, [0.5 1], [0.5 1], '--', 'Color', c.ideal, 'LineWidth', 1.0, ...
         'DisplayName', 'perfect calibration');

    plot(ax, R.nominal, R.empGP, '-o', 'Color', c.gp, 'LineWidth', 1.2, ...
         'MarkerSize', 4, 'MarkerFaceColor', c.gp, ...
         'DisplayName', 'GP posterior (Helmholtz prior)');
    plot(ax, R.nominal, R.empCP, '-s', 'Color', c.conformal, 'LineWidth', 1.2, ...
         'MarkerSize', 4, 'MarkerFaceColor', c.conformal, ...
         'DisplayName', 'split conformal');
    plot(ax, R.nominal, R.empNCP, '-^', 'Color', c.normalized, 'LineWidth', 1.2, ...
         'MarkerSize', 4, 'MarkerFaceColor', c.normalized, ...
         'DisplayName', 'normalized conformal');

    xlabel(ax, 'Nominal coverage');
    ylabel(ax, 'Empirical coverage');
    xlim(ax, [0.5 1]); ylim(ax, [0.5 1]);
    axis(ax, 'square');
    legend(ax, 'Location', 'southeast');
    sfuq.viz.style(ax);

    if nargin > 1 && ~isempty(outName)
        sfuq.viz.saveFig(fig, outName);
    end
end
