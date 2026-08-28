function R = exp07_support_failure(nTrials)
%EXP07_SUPPORT_FAILURE  Where reweighting provably cannot help, and why that is good.
%
%   Regime B: the listener positions to be certified extend BEYOND the region
%   the calibration microphones cover. The query density is positive where the
%   calibration density is zero, so absolute continuity fails and the
%   likelihood ratio is genuinely infinite there.
%
%   This experiment exists because the honest answer is a negative one, and
%   stating it precisely is more valuable than hiding it.
%
%   What each method does when asked the impossible
%   -----------------------------------------------
%   Split conformal returns a finite, confident radius and is badly wrong:
%   empirical coverage collapses toward 0.28 against a nominal 0.90, with no
%   warning of any kind. The failure is silent, which is the worst property an
%   uncertainty estimate can have.
%
%   Weighted conformal, given the correct infinite ratio, returns an INFINITE
%   interval. It abstains. It does not pretend the guarantee survives, because
%   no reweighting can manufacture information about a region no measurement
%   constrains. Abstention is the correct behaviour, and the abstention rate is
%   a directly actionable output: it maps exactly which listener positions the
%   current microphone layout cannot certify.
%
%   The practical reading is not "weighted conformal fails here". It is: a
%   method that tells you it cannot answer is safe to deploy, and a method that
%   answers confidently and wrongly is not.
%
%   Produces figure 9.

    if nargin < 1 || isempty(nTrials), nTrials = 8; end

    S     = sfuq.data.makeScenario();
    nMic  = 120;
    alpha = 0.10;
    ratios = [1.0 1.25 1.5 2.0 2.5];

    P = sfuq.geom.uniformBox(S.Centre, S.HalfWidth);

    fprintf('\n=== Experiment 07: query region beyond calibration support ===\n');
    fprintf('  %8s %10s %10s %12s %12s\n', ...
            'ratio', 'CP cov', 'wCP cov', 'abstain', 'certifiable');

    n = numel(ratios);
    [covCP, covWC, abst, certF] = deal(zeros(n, nTrials));

    for i = 1:n
        Q = sfuq.geom.uniformBox(S.Centre, ratios(i) * S.HalfWidth);

        for t = 1:nTrials
            T = sfuq.shiftTrial(S, nMic, P, Q, ...
                    struct('seed', 7000 + 100*i + t));

            rCP = sfuq.uq.splitConformal(T.resCal, alpha);
            [rWC, info] = sfuq.uq.weightedConformal(T.resCal, T.wCal, T.wQuery, alpha);

            % Split conformal is scored over ALL query points: it makes a claim
            % everywhere, so it is judged everywhere.
            covCP(i, t) = sfuq.eval.coverage(T.resQuery, rCP);

            % Weighted conformal is scored only where it actually made a claim.
            % Scoring an abstention as a success would be meaningless -- an
            % infinite interval covers by definition.
            fin = isfinite(rWC);
            if any(fin)
                covWC(i, t) = sfuq.eval.coverage(T.resQuery(fin), rWC(fin));
            else
                covWC(i, t) = NaN;
            end

            abst(i, t)  = info.abstentionRate;
            certF(i, t) = mean(T.certifiable);
        end

        fprintf('  %8.2f %10.3f %10.3f %12.3f %12.3f\n', ratios(i), ...
                mean(covCP(i,:)), mean(covWC(i,:), 'omitnan'), ...
                mean(abst(i,:)), mean(certF(i,:)));
    end

    R = struct('ratios', ratios(:), 'alpha', alpha, 'nMic', nMic, ...
               'covCP', mean(covCP, 2), ...
               'covWCP', mean(covWC, 2, 'omitnan'), ...
               'abstention', mean(abst, 2), ...
               'certifiable', mean(certF, 2), 'scenario', S);

    fprintf('\n  Split conformal at ratio %.1f: %.3f coverage, 0.000 abstention.\n', ...
            ratios(end), R.covCP(end));
    fprintf('  It never declines to answer, and it is wrong %.0f%% of the time.\n', ...
            100 * (1 - R.covCP(end)));
    fprintf('  Weighted conformal abstains on %.0f%% of those points instead.\n', ...
            100 * R.abstention(end));

    sfuq.viz.figAbstention(R, 'fig09_abstention');
    save(fullfile(sfuq.viz.projectRoot(), 'results', 'exp07_support.mat'), '-struct', 'R');
end
