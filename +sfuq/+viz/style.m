function style(ax)
%STYLE  Consistent publication styling for every figure in the project.
%
%   Sized and weighted for an IEEE two-column layout reduced to 3.5 inches.
%   Call at the end of each figure function so the whole paper speaks with one
%   visual voice.

    if nargin < 1 || isempty(ax)
        ax = gca;
    end

    set(ax, 'FontName', 'Times New Roman', 'FontSize', 9, ...
            'Box', 'off', 'TickDir', 'out', 'LineWidth', 0.75, ...
            'XGrid', 'on', 'YGrid', 'on', 'GridAlpha', 0.12, 'Layer', 'top');

    if ~isempty(ax.XLabel), ax.XLabel.FontSize = 10; end
    if ~isempty(ax.YLabel), ax.YLabel.FontSize = 10; end
    if ~isempty(ax.Title)
        ax.Title.FontSize = 10;
        ax.Title.FontWeight = 'normal';
    end

    lg = findobj(ancestor(ax, 'figure'), 'Type', 'Legend');
    if ~isempty(lg)
        set(lg, 'FontSize', 8, 'Box', 'off');
    end
end
