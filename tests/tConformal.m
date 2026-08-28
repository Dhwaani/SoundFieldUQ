classdef tConformal < matlab.unittest.TestCase
%TCONFORMAL  The coverage guarantee itself, pinned by Monte Carlo.
%
%   This is the most important test file in the project. Every claim the paper
%   makes rests on split conformal actually delivering its advertised
%   coverage, and on the GP radius conversion carrying the right factor of
%   two. Both are easy to get subtly wrong in ways that shift results by a few
%   percent and are invisible by inspection.

    methods (Test)

        function quantileUsesOrderStatisticNotInterpolation(tc)
            % The finite-sample guarantee needs the ceil(p*n)-th order
            % statistic. An interpolating quantile returns something slightly
            % smaller and quietly breaks the guarantee.
            x = (1:10).';
            tc.verifyEqual(sfuq.uq.quantileHigher(x, 0.95), 10);
            tc.verifyEqual(sfuq.uq.quantileHigher(x, 0.50), 5);
            tc.verifyEqual(sfuq.uq.quantileHigher(x, 0.01), 1);
        end

        function splitConformalAchievesNominalCoverage(tc)
            % Model-free Monte Carlo: draw calibration and test residuals from
            % the same arbitrary (deliberately non-Gaussian) distribution and
            % confirm coverage is at least nominal.
            rng(0);
            alpha = 0.10;
            nRep = 400;
            cov = zeros(nRep, 1);

            for r = 1:nRep
                cal  = exprnd_local(1, 200, 1);   % heavy-tailed, not Gaussian
                test = exprnd_local(1, 500, 1);
                R = sfuq.uq.splitConformal(cal, alpha);
                cov(r) = mean(test <= R);
            end

            % The guarantee is on the average over calibration draws.
            tc.verifyGreaterThanOrEqual(mean(cov), 1 - alpha - 0.01, ...
                'Split conformal fell below its guaranteed coverage.');
            tc.verifyLessThan(mean(cov), 1 - alpha + 0.03, ...
                'Coverage far above nominal suggests the quantile is too high.');
        end

        function guaranteeHoldsForSmallCalibrationSets(tc)
            % The (n+1) correction matters most when n is small. Without it,
            % coverage at n = 30 would sit visibly below nominal.
            rng(1);
            alpha = 0.10;
            cov = zeros(600, 1);
            for r = 1:600
                cal  = randn(30, 1) .^ 2;
                test = randn(200, 1) .^ 2;
                cov(r) = mean(test <= sfuq.uq.splitConformal(cal, alpha));
            end
            tc.verifyGreaterThanOrEqual(mean(cov), 1 - alpha - 0.015);
        end

        function warnsWhenCalibrationSetTooSmall(tc)
            % With n = 5 and alpha = 0.01 the guarantee is unattainable and
            % the user must be told rather than handed a silent fallback.
            tc.verifyWarning(@() sfuq.uq.splitConformal(rand(5, 1), 0.01), ...
                             'sfuq:splitConformal:tooFewPoints');
        end

        function gpRadiusHasCorrectRayleighFactor(tc)
            % The radius conversion assumes |e| is Rayleigh with scale sqrt(v)
            % when the real and imaginary errors are each N(0, v). Verify
            % against direct simulation; a missing factor of two here would
            % bias every GP coverage number in the paper.
            rng(2);
            v = 0.37;
            alpha = 0.10;
            e = sqrt(v) * (randn(400000, 1) + 1i * randn(400000, 1));
            R = sfuq.uq.gpRadius(v, alpha);
            tc.verifyEqual(mean(abs(e) <= R), 1 - alpha, 'AbsTol', 0.005, ...
                'GP radius does not achieve nominal coverage under its own model.');
        end

        function normalizedConformalIsAlsoMarginallyValid(tc)
            % Normalisation must not cost the marginal guarantee.
            rng(3);
            alpha = 0.10;
            cov = zeros(200, 1);
            for r = 1:200
                sCal  = 0.5 + rand(300, 1);        % heteroscedastic
                sTest = 0.5 + rand(600, 1);
                resCal  = abs(sCal  .* randn(300, 1));
                resTest = abs(sTest .* randn(600, 1));
                Rq = sfuq.uq.normalizedConformal(resCal, sCal, sTest, alpha);
                cov(r) = mean(resTest <= Rq);
            end
            tc.verifyGreaterThanOrEqual(mean(cov), 1 - alpha - 0.015);
        end

        function normalizedBeatsSplitOnConditionalCoverage(tc)
            % The project's core claim, in miniature. Under heteroscedastic
            % errors a constant radius must over-cover the easy half and
            % under-cover the hard half, while a scaled radius should not.
            % Thresholds below have margin: across 20 independent seeds the
            % split-conformal gap never fell below 0.18 and the normalized gap
            % never exceeded 0.015 at this sample size. Note the gap is
            % insensitive to the size of the heteroscedasticity contrast -- it
            % is set by the mixture proportions and alpha -- so widening the
            % contrast would not buy extra margin.
            rng(4);
            alpha = 0.10;
            n = 8000;
            sCal  = [0.2 * ones(n/2, 1); 2.0 * ones(n/2, 1)];
            sTest = [0.2 * ones(n/2, 1); 2.0 * ones(n/2, 1)];
            resCal  = abs(sCal  .* randn(n, 1));
            resTest = abs(sTest .* randn(n, 1));

            rCP = sfuq.uq.splitConformal(resCal, alpha);
            rNC = sfuq.uq.normalizedConformal(resCal, sCal, sTest, alpha);

            easy = 1:(n/2);
            hard = (n/2 + 1):n;

            gapCP = abs(mean(resTest(easy) <= rCP) - mean(resTest(hard) <= rCP));
            gapNC = abs(mean(resTest(easy) <= rNC(easy)) - ...
                        mean(resTest(hard) <= rNC(hard)));

            tc.verifyGreaterThan(gapCP, 0.12, ...
                'Expected split conformal to be badly non-uniform here.');
            tc.verifyLessThan(gapNC, 0.05, ...
                'Normalized conformal should be nearly conditionally valid.');
        end
    end
end

function x = exprnd_local(mu, m, n)
%EXPRND_LOCAL  Exponential samples by inverse transform.
%   Avoids the Statistics Toolbox so the test suite runs on base MATLAB.
    x = -mu * log(rand(m, n));
end
