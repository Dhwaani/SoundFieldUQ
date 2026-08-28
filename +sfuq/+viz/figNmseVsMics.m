function fig = figNmseVsMics(R, outName)
%FIGNMSEVSMICS  Reconstruction accuracy against microphone count.
%
%   Establishes that the estimators work before any uncertainty claim is made.
%   The Helmholtz kernel should beat the physics-free squared-exponential
%   control at every microphone count, and the size of that gap is the value
%   of encoding the wave equation in the prior.
%
%   Input R is the struct produced by EXP01_ACCURACY.

    c = sfuq.viz.colors();
    fig = figure('Units', 'inches', 'Position', [1 1 3.5 2.6], 'Color', 'w');
    ax = axes(fig); hold(ax, 'on');

    errorbar(ax, R.nMics, mean(R.nmseHelmholtz, 2), std(R.nmseHelmholtz, 0, 2), ...
             '-o', 'Color', c.gp, 'LineWidth', 1.2, 'MarkerFaceColor', c.gp, ...
             'MarkerSize', 4, 'CapSize', 3, 'DisplayName', 'Helmholtz kernel');
    errorbar(ax, R.nMics, mean(R.nmseGaussian, 2), std(R.nmseGaussian, 0, 2), ...
             '-s', 'Color', c.krr, 'LineWidth', 1.2, 'MarkerFaceColor', c.krr, ...
             'MarkerSize', 4, 'CapSize', 3, 'DisplayName', 'squared exponential');

    set(ax, 'XScale', 'log');
    xlabel(ax, 'Number of microphones');
    ylabel(ax, 'NMSE (dB)');
    title(ax, sprintf('%g Hz, T_{60} proxy \\beta = %g', R.frequency, R.beta));
    legend(ax, 'Location', 'northeast');
    sfuq.viz.style(ax);

    if nargin > 1 && ~isempty(outName)
        sfuq.viz.saveFig(fig, outName);
    end
end
