function fig = figConditionalCoverage(R, outName)
%FIGCONDITIONALCOVERAGE  Coverage as a function of distance from the aperture.
%
%   The figure that carries the project's main finding. Marginal coverage is
%   an average and can look perfect while the spatial breakdown is badly
%   uneven; this plot exposes that. The shaded band marks the target level.
%
%   Input R is the struct produced by EXP03_CONDITIONAL_COVERAGE.

    c = sfuq.viz.colors();
    fig = figure('Units', 'inches', 'Position', [1 1 3.5 2.6], 'Color', 'w');
    ax = axes(fig); hold(ax, 'on');

    % Force column orientation. The experiments build edges as a row and
    % coverage as a column; mixing the two silently transposes some plots.
    centres = ((R.edges(1:end-1) + R.edges(2:end)) / 2).';
    covGP  = R.covGP(:);
    covCP  = R.covCP(:);
    covNCP = R.covNCP(:);
    target = 1 - R.alpha;

    % A +/- 2 standard-error band around the target, so a reader can judge
    % whether a deviation is meaningful or just finite-sample noise.
    se = 2 * sqrt(target * (1 - target) / R.nPerBin);
    fill(ax, [centres; flipud(centres)], ...
             [repmat(target - se, numel(centres), 1); ...
              repmat(target + se, numel(centres), 1)], ...
         [0.85 0.85 0.85], 'EdgeColor', 'none', 'FaceAlpha', 0.6, ...
         'DisplayName', 'target \pm 2 s.e.');
    yline(ax, target, '--', 'Color', c.ideal, 'LineWidth', 1.0, ...
          'HandleVisibility', 'off');

    plot(ax, centres, covGP,  '-o', 'Color', c.gp, 'LineWidth', 1.2, ...
         'MarkerFaceColor', c.gp, 'MarkerSize', 4, 'DisplayName', 'GP posterior');
    plot(ax, centres, covCP,  '-s', 'Color', c.conformal, 'LineWidth', 1.2, ...
         'MarkerFaceColor', c.conformal, 'MarkerSize', 4, ...
         'DisplayName', 'split conformal');
    plot(ax, centres, covNCP, '-^', 'Color', c.normalized, 'LineWidth', 1.2, ...
         'MarkerFaceColor', c.normalized, 'MarkerSize', 4, ...
         'DisplayName', 'normalized conformal');

    if isfield(R, 'apertureHalfWidth')
        xline(ax, R.apertureHalfWidth, ':', 'aperture edge', 'FontSize', 7, ...
              'Color', c.warn, 'LabelOrientation', 'horizontal', ...
              'HandleVisibility', 'off');
    end

    xlabel(ax, 'Distance from aperture centre (m)');
    ylabel(ax, 'Empirical coverage');
    ylim(ax, [0 1.02]);
    legend(ax, 'Location', 'southwest');
    sfuq.viz.style(ax);

    if nargin > 1 && ~isempty(outName)
        sfuq.viz.saveFig(fig, outName);
    end
end
