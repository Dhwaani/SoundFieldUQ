function R = exp02_marginal_calibration(nTrials)
%EXP02_MARGINAL_CALIBRATION  Are the uncertainty estimates honest on average?
%
%   Sweeps the nominal coverage level from 50% to 99% and measures the
%   coverage actually achieved by three methods:
%
%     GP posterior      Bayesian, exact if the model is correct, no guarantee
%                       if it is not
%     split conformal   distribution-free, guaranteed marginally, constant
%                       radius everywhere
%     normalized        distribution-free, guaranteed marginally, radius
%     conformal         shaped by the GP's own variance
%
%   Expected result, and why it is worth stating up front
%   ----------------------------------------------------
%   A properly fitted Helmholtz GP is close to calibrated for interpolation
%   inside the aperture. That is not the "GP uncertainty is broken" story one
%   might expect, and it is worth reporting honestly: when the prior genuinely
%   matches the physics and its scale is fitted rather than assumed, the
%   Bayesian intervals are good. The interesting failures appear in the
%   conditional analysis (experiment 03) and under aperture shift
%   (experiment 04), not here.
%
%   Reporting this negative-looking result strengthens rather than weakens the
%   paper: it shows the comparison was run fairly.
%
%   Produces figure 2.

    if nargin < 1 || isempty(nTrials)
        nTrials = 10;
    end

    S = sfuq.data.makeScenario();
    nMic = 120;
    levels = (0.50:0.05:0.95).';

    fprintf('\n=== Experiment 02: marginal calibration ===\n');

    empGP = zeros(numel(levels), nTrials);
    empCP = zeros(numel(levels), nTrials);
    empNC = zeros(numel(levels), nTrials);

    for t = 1:nTrials
        T = sfuq.trial(S, nMic, struct('seed', 2000 + t));
        [~, empGP(:, t)] = sfuq.eval.reliabilityCurve( ...
            T.resCal, T.sigCal, T.resTest, T.sigTest, 'gp', levels);
        [~, empCP(:, t)] = sfuq.eval.reliabilityCurve( ...
            T.resCal, T.sigCal, T.resTest, T.sigTest, 'conformal', levels);
        [~, empNC(:, t)] = sfuq.eval.reliabilityCurve( ...
            T.resCal, T.sigCal, T.resTest, T.sigTest, 'normalized', levels);
    end

    R = struct();
    R.nominal = levels;
    R.empGP   = mean(empGP, 2);
    R.empCP   = mean(empCP, 2);
    R.empNCP  = mean(empNC, 2);
    R.scenario = S;
    R.nMic = nMic;

    % Calibration error: mean absolute deviation from the diagonal, the single
    % number that summarises the whole curve.
    R.eceGP  = mean(abs(R.empGP  - levels));
    R.eceCP  = mean(abs(R.empCP  - levels));
    R.eceNCP = mean(abs(R.empNCP - levels));

    fprintf('  mean |empirical - nominal| over the sweep:\n');
    fprintf('    GP posterior        %.4f\n', R.eceGP);
    fprintf('    split conformal     %.4f\n', R.eceCP);
    fprintf('    normalized conformal %.4f\n', R.eceNCP);

    sfuq.viz.figReliability(R, 'fig02_reliability');
    save(fullfile(sfuq.viz.projectRoot(), 'results', 'exp02_marginal.mat'), '-struct', 'R');
end
