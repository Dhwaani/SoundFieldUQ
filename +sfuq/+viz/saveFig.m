function saveFig(fig, name, outDir)
%SAVEFIG  Write a figure as vector PDF and 600-dpi PNG side by side.
%
%   Vector for the paper, raster for the README and slides, written together
%   so the two can never drift apart.
%
%   Requires R2020a or newer for EXPORTGRAPHICS; falls back to PRINT on older
%   releases rather than failing outright.

    if nargin < 3 || isempty(outDir)
        outDir = fullfile(sfuq.viz.projectRoot(), 'figures');
    end
    if ~isfolder(outDir)
        mkdir(outDir);
    end

    base = fullfile(outDir, name);

    if exist('exportgraphics', 'file') == 2
        exportgraphics(fig, [base '.pdf'], 'ContentType', 'vector', ...
                       'BackgroundColor', 'white');
        exportgraphics(fig, [base '.png'], 'Resolution', 600, ...
                       'BackgroundColor', 'white');
    else
        set(fig, 'PaperPositionMode', 'auto', 'InvertHardcopy', 'off', ...
                 'Color', 'w');
        print(fig, base, '-dpdf', '-painters');
        print(fig, base, '-dpng', '-r600');
    end

    fprintf('  wrote %s.{pdf,png}\n', base);
end
