function fig = figPlacement(R, outName)
%FIGPLACEMENT  Layouts differ in efficiency, not in validity.
%
%   Grouped bars across three calibration layouts. Coverage is essentially
%   identical -- weighted conformal is valid for all of them, which is the
%   point of a distribution-free guarantee. What differs is effective sample
%   size and interval width.
%
%   This is the figure that carries the practical recommendation: put the
%   calibration microphones where the listeners are. A compact array is not
%   invalid, it is expensive, and the expense is invisible unless the effective
%   sample size is reported.

    if nargin < 2, outName = ''; end
    c = sfuq.viz.colors();

    fig = figure('Units', 'inches', 'Position', [1 1 3.5 2.7], 'Color', 'w');
    ax = axes(fig); hold(ax, 'on');

    % Width normalised to the matched design so the two series share an axis.
    relWidth = R.widthWCP / R.widthWCP(end);
    data = [R.covWCP(:), R.essFrac(:), relWidth(:)];

    b = bar(ax, data, 'grouped', 'EdgeColor', 'k', 'LineWidth', 0.4);
    b(1).FaceColor = c.weighted;
    b(2).FaceColor = c.gp;
    b(3).FaceColor = c.krr;

    yline(ax, 1 - R.alpha, '--', 'nominal coverage', 'Color', c.ideal, ...
          'FontSize', 7, 'LabelHorizontalAlignment', 'left', ...
          'HandleVisibility', 'off');

    set(ax, 'XTick', 1:numel(R.designNames), 'XTickLabel', R.designNames);
    ylabel(ax, 'value');
    legend(ax, {'coverage', 'ESS / n', 'width (rel. to matched)'}, ...
           'Location', 'northwest');
    title(ax, 'All layouts valid; they differ in cost');
    sfuq.viz.style(ax);

    if ~isempty(outName), sfuq.viz.saveFig(fig, outName); end
end
