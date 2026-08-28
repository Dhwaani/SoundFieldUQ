function fig = figSpatialMap(R, outName)
%FIGSPATIALMAP  Slice through the region showing prediction radius in space.
%
%   Two panels on a shared colour scale: the constant radius that split
%   conformal assigns everywhere, and the spatially varying radius that
%   normalized conformal assigns. Microphone positions are overlaid. This is
%   the figure that makes the abstract argument about conditional validity
%   immediately visible.
%
%   Input R is the struct produced by EXP03_CONDITIONAL_COVERAGE.

    c = sfuq.viz.colors();
    fig = figure('Units', 'inches', 'Position', [1 1 6.5 2.6], 'Color', 'w');
    tl = tiledlayout(fig, 1, 2, 'TileSpacing', 'compact', 'Padding', 'compact');

    lims = [0, max([R.mapCP(:); R.mapNCP(:)])];

    panels = {R.mapCP, R.mapNCP};
    names  = {'split conformal (constant)', 'normalized conformal (adaptive)'};

    for i = 1:2
        ax = nexttile(tl);
        imagesc(ax, R.mapX, R.mapY, panels{i});
        set(ax, 'YDir', 'normal', 'CLim', lims);
        hold(ax, 'on');
        plot(ax, R.micXY(:, 1), R.micXY(:, 2), 'w.', 'MarkerSize', 7);
        plot(ax, R.micXY(:, 1), R.micXY(:, 2), 'k.', 'MarkerSize', 4);
        rectangle(ax, 'Position', [R.apertureBox(1), R.apertureBox(3), ...
                                   diff(R.apertureBox(1:2)), ...
                                   diff(R.apertureBox(3:4))], ...
                  'EdgeColor', c.warn, 'LineStyle', ':', 'LineWidth', 1.0);
        title(ax, names{i});
        xlabel(ax, 'x (m)');
        if i == 1
            ylabel(ax, 'y (m)');
        end
        axis(ax, 'image');
        sfuq.viz.style(ax);
    end

    cb = colorbar(ax);
    cb.Label.String = 'prediction radius';
    cb.Label.FontName = 'Times New Roman';

    if nargin > 1 && ~isempty(outName)
        sfuq.viz.saveFig(fig, outName);
    end
end
