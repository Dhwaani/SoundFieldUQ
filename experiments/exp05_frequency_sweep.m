function R = exp05_frequency_sweep(nTrials)
%EXP05_FREQUENCY_SWEEP  Behaviour across the spatial sampling limit.
%
%   The Helmholtz kernel's effective correlation length is set by the
%   wavelength. As frequency rises, a fixed microphone layout samples the
%   field ever more sparsely, and above roughly half-wavelength spacing the
%   reconstruction problem becomes ill-posed. This experiment tracks both
%   accuracy and calibration through that transition.
%
%   It matters for the paper because it separates two very different kinds of
%   failure. Reconstruction NMSE degrading with frequency is expected physics
%   and not a defect. Uncertainty estimates staying confident while the
%   reconstruction falls apart WOULD be a defect, and a serious one: it is the
%   regime where a practitioner is most likely to be misled.
%
%   Produces figure 6.

    if nargin < 1 || isempty(nTrials)
        nTrials = 8;
    end

    freqs = [200 400 800 1200 1600 2000];
    nMic  = 120;
    alpha = 0.10;

    fprintf('\n=== Experiment 05: frequency sweep ===\n');
    fprintf('  %8s %10s %9s %9s %9s %9s\n', ...
            'f (Hz)', 'spacing', 'NMSE dB', 'GP cov', 'CP cov', 'nCP cov');

    nmse  = zeros(numel(freqs), nTrials);
    covGP = zeros(numel(freqs), nTrials);
    covCP = zeros(numel(freqs), nTrials);
    covNC = zeros(numel(freqs), nTrials);
    spacing = zeros(numel(freqs), 1);

    for i = 1:numel(freqs)
        S = sfuq.data.makeScenario('Frequency', freqs(i));

        % Mean nearest-neighbour spacing for a uniform layout in the cube,
        % compared against the half-wavelength criterion.
        vol = (2 * S.HalfWidth) ^ 3;
        spacing(i) = (vol / nMic) ^ (1 / 3) / S.nyquistSpacing;

        for t = 1:nTrials
            T = sfuq.trial(S, nMic, struct('seed', 5000 + 100 * i + t, ...
                                           'nTest', 1000));
            rGP = sfuq.uq.gpRadius(T.sigTest .^ 2, alpha);
            rCP = sfuq.uq.splitConformal(T.resCal, alpha);
            rNC = sfuq.uq.normalizedConformal(T.resCal, T.sigCal, T.sigTest, alpha);

            nmse(i, t)  = T.nmse;
            covGP(i, t) = sfuq.eval.coverage(T.resTest, rGP);
            covCP(i, t) = sfuq.eval.coverage(T.resTest, rCP);
            covNC(i, t) = sfuq.eval.coverage(T.resTest, rNC);
        end

        fprintf('  %8d %10.2f %9.2f %9.3f %9.3f %9.3f\n', freqs(i), ...
                spacing(i), mean(nmse(i, :)), mean(covGP(i, :)), ...
                mean(covCP(i, :)), mean(covNC(i, :)));
    end

    R = struct('freqs', freqs(:), 'spacingRatio', spacing, 'alpha', alpha, ...
               'nmse', mean(nmse, 2), 'covGP', mean(covGP, 2), ...
               'covCP', mean(covCP, 2), 'covNCP', mean(covNC, 2), 'nMic', nMic);

    c = sfuq.viz.colors();
    fig = figure('Units', 'inches', 'Position', [1 1 3.5 3.4], 'Color', 'w');
    tl = tiledlayout(fig, 2, 1, 'TileSpacing', 'compact', 'Padding', 'compact');

    ax1 = nexttile(tl); hold(ax1, 'on');
    plot(ax1, R.freqs, R.nmse, '-o', 'Color', c.krr, 'LineWidth', 1.2, ...
         'MarkerFaceColor', c.krr, 'MarkerSize', 4);
    ylabel(ax1, 'NMSE (dB)');
    set(ax1, 'XTickLabel', []);
    sfuq.viz.style(ax1);

    ax2 = nexttile(tl); hold(ax2, 'on');
    yline(ax2, 1 - alpha, '--', 'Color', c.ideal, 'DisplayName', 'nominal');
    plot(ax2, R.freqs, R.covGP,  '-o', 'Color', c.gp, 'LineWidth', 1.2, ...
         'MarkerFaceColor', c.gp, 'MarkerSize', 4, 'DisplayName', 'GP');
    plot(ax2, R.freqs, R.covCP,  '-s', 'Color', c.conformal, 'LineWidth', 1.2, ...
         'MarkerFaceColor', c.conformal, 'MarkerSize', 4, 'DisplayName', 'conformal');
    plot(ax2, R.freqs, R.covNCP, '-^', 'Color', c.normalized, 'LineWidth', 1.2, ...
         'MarkerFaceColor', c.normalized, 'MarkerSize', 4, 'DisplayName', 'normalized');
    xlabel(ax2, 'Frequency (Hz)');
    ylabel(ax2, 'Coverage');
    ylim(ax2, [0 1.02]);
    legend(ax2, 'Location', 'southwest');
    sfuq.viz.style(ax2);

    sfuq.viz.saveFig(fig, 'fig06_frequency_sweep');
    save(fullfile(sfuq.viz.projectRoot(), 'results', 'exp05_frequency.mat'), '-struct', 'R');
end
