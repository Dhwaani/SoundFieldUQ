classdef tWeightedConformal < matlab.unittest.TestCase
%TWEIGHTEDCONFORMAL  The new guarantee, pinned by Monte Carlo.
%
%   Everything the project now claims rests on weighted conformal delivering
%   its advertised coverage under a KNOWN likelihood ratio, and on it abstaining
%   rather than lying when the ratio is infinite. Both are tested here against
%   constructed cases with known ground truth.

    methods (Test)

        function reducesToSplitConformalWhenWeightsAreEqual(tc)
            % With w identical everywhere there is no shift, and the weighted
            % construction must collapse onto ordinary split conformal. If this
            % fails, the normalisation or the atom at infinity is wrong.
            rng(0);
            res = abs(randn(300, 1));
            alpha = 0.10;
            w = ones(300, 1);

            rSplit = sfuq.uq.splitConformal(res, alpha);
            rW = sfuq.uq.weightedConformal(res, w, ones(50, 1), alpha);

            tc.verifyEqual(unique(rW), rSplit, 'RelTol', 1e-12, ...
                'Equal weights must reproduce split conformal exactly.');
        end

        function achievesNominalCoverageUnderKnownShift(tc)
            % The central guarantee. Calibration is drawn from a density that
            % over-samples the easy region; the test set is uniform. Unweighted
            % conformal must under-cover, weighted must not.
            rng(1);
            alpha = 0.10;
            nRep = 250;
            covPlain = zeros(nRep, 1);
            covW = zeros(nRep, 1);

            for r = 1:nRep
                % Covariate u in [0,1]; error scale grows with u, so large u is
                % "hard". Calibration over-samples small u (Beta-like via u^2).
                uCal  = rand(400, 1) .^ 2;      % concentrated near 0
                uTest = rand(800, 1);           % uniform

                sd = @(u) 0.2 + 2 * u;
                resCal  = abs(sd(uCal)  .* randn(400, 1));
                resTest = abs(sd(uTest) .* randn(800, 1));

                % p(u) = 1/(2 sqrt(u)) for u = v^2 with v uniform; q(u) = 1.
                wOf = @(u) 2 * sqrt(u);
                rPlain = sfuq.uq.splitConformal(resCal, alpha);
                rW = sfuq.uq.weightedConformal(resCal, wOf(uCal), wOf(uTest), alpha);

                covPlain(r) = mean(resTest <= rPlain);
                fin = isfinite(rW);
                covW(r) = mean(resTest(fin) <= rW(fin));
            end

            tc.verifyLessThan(mean(covPlain), 1 - alpha - 0.03, ...
                'Unweighted conformal should visibly under-cover under this shift.');
            tc.verifyGreaterThanOrEqual(mean(covW), 1 - alpha - 0.02, ...
                'Weighted conformal lost its coverage guarantee.');
            tc.verifyLessThan(mean(covW), 1 - alpha + 0.06, ...
                'Weighted conformal is far over-covering; check the weights.');
        end

        function abstainsRatherThanUnderCovering(tc)
            % An infinite query weight means the calibration set says nothing
            % about that point. The radius must be Inf, not a finite guess.
            rng(2);
            res = abs(randn(200, 1));
            wCal = ones(200, 1);
            wQuery = [ones(10, 1); inf(10, 1)];

            [R, info] = sfuq.uq.weightedConformal(res, wCal, wQuery, 0.10);

            tc.verifyTrue(all(isfinite(R(1:10))));
            tc.verifyTrue(all(isinf(R(11:20))));
            tc.verifyEqual(info.abstentionRate, 0.5, 'AbsTol', 1e-12);
        end

        function extremeWeightConcentrationForcesAbstention(tc)
            % When a single query point dominates the total weight, its own
            % atom at infinity exceeds alpha and the quantile is genuinely
            % infinite. The method must not silently truncate to the largest
            % calibration score.
            res = abs(randn(50, 1));
            wCal = ones(50, 1) * 1e-6;
            wQuery = 1e6;

            R = sfuq.uq.weightedConformal(res, wCal, wQuery, 0.10);
            tc.verifyTrue(isinf(R), ...
                'A dominating query weight must produce an infinite interval.');
        end

        function effectiveSampleSizeBehavesCorrectly(tc)
            tc.verifyEqual(sfuq.uq.effectiveSampleSize(ones(100, 1)), 100, ...
                'RelTol', 1e-12);

            [ess, frac] = sfuq.uq.effectiveSampleSize([1; zeros(99, 1)]);
            tc.verifyEqual(ess, 1, 'RelTol', 1e-12);
            tc.verifyEqual(frac, 0.01, 'RelTol', 1e-12);

            % Concentrated weights must give ESS strictly below n.
            w = exp(3 * rand(200, 1));
            tc.verifyLessThan(sfuq.uq.effectiveSampleSize(w), 200);

            tc.verifyEqual(sfuq.uq.effectiveSampleSize([1; Inf]), 0, ...
                'Infinite weights mean no usable effective sample.');
        end

        function rejectsInfiniteCalibrationWeights(tc)
            % An infinite CALIBRATION weight is contradictory: the point exists
            % in the calibration sample, so its density there cannot be zero.
            % Usually it means the two densities were passed the wrong way round.
            tc.verifyError( ...
                @() sfuq.uq.weightedConformal(rand(10,1), [ones(9,1); Inf], 1, 0.1), ...
                'sfuq:weightedConformal:infiniteCalWeight');
        end

        function widerWeightSpreadCostsEffectiveSampleSize(tc)
            % The efficiency claim reported in figure 8.
            rng(3);
            narrow = exp(0.2 * randn(500, 1));
            wide   = exp(1.5 * randn(500, 1));
            tc.verifyGreaterThan(sfuq.uq.effectiveSampleSize(narrow), ...
                                 sfuq.uq.effectiveSampleSize(wide));
        end
    end
end
