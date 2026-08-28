function R = exp01_accuracy(nTrials)
%EXP01_ACCURACY  Reconstruction accuracy versus microphone count.
%
%   Establishes that the estimators work at all, before any uncertainty claim
%   is made. Two kernels are compared: the Helmholtz kernel, which encodes the
%   wave equation, and a squared-exponential kernel of comparable length
%   scale, which does not. The gap between them is the value of the physics.
%
%   This experiment is deliberately first. If the reconstruction NMSE does not
%   fall with microphone count, or the physics-informed kernel does not beat
%   the generic one, something is wrong with the data generation or the solver
%   and no amount of careful uncertainty analysis downstream would be
%   meaningful.
%
%   Produces figure 1.

    if nargin < 1 || isempty(nTrials)
        nTrials = 8;
    end

    S = sfuq.data.makeScenario();
    nMics = [20 40 80 160 320];

    fprintf('\n=== Experiment 01: reconstruction accuracy ===\n');
    fprintf('  %g Hz, wavelength %.3f m, half-wavelength spacing %.3f m\n', ...
            S.Frequency, S.wavelength, S.nyquistSpacing);
    fprintf('  aperture %.1f m cube, %d image sources\n', ...
            2 * S.HalfWidth, size(S.imgPos, 1));

    nmseH = zeros(numel(nMics), nTrials);
    nmseG = zeros(numel(nMics), nTrials);

    % Length scale for the control kernel, matched to the first zero of the
    % Helmholtz kernel so the comparison is about the SHAPE of the prior
    % rather than about an arbitrary bandwidth difference.
    ell = pi / S.k;

    for i = 1:numel(nMics)
        for t = 1:nTrials
            o = struct('seed', 1000 * i + t, 'nTest', 500, 'nCal', 50);

            Th = sfuq.trial(S, nMics(i), o);
            nmseH(i, t) = Th.nmse;

            o.kernelFcn = @(a, b) sfuq.kernels.gaussian(a, b, ell);
            Tg = sfuq.trial(S, nMics(i), o);
            nmseG(i, t) = Tg.nmse;
        end
        fprintf('  %4d mics: Helmholtz %7.2f dB   squared-exp %7.2f dB\n', ...
                nMics(i), mean(nmseH(i, :)), mean(nmseG(i, :)));
    end

    R = struct('nMics', nMics(:), 'nmseHelmholtz', nmseH, ...
               'nmseGaussian', nmseG, 'frequency', S.Frequency, ...
               'beta', S.Beta, 'scenario', S, 'ell', ell);

    sfuq.viz.figNmseVsMics(R, 'fig01_nmse_vs_mics');
    save(fullfile(sfuq.viz.projectRoot(), 'results', 'exp01_accuracy.mat'), '-struct', 'R');
end
