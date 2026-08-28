classdef tGeometry < matlab.unittest.TestCase
%TGEOMETRY  Densities and the exact likelihood ratio.
%
%   The ratio being EXACT rather than estimated is the whole claim to novelty,
%   so these tests check it against closed-form values rather than against
%   another implementation.

    methods (Test)

        function uniformBoxIntegratesToOne(tc)
            D = sfuq.geom.uniformBox([0 0 0], 0.5);
            tc.verifyEqual(D.pdf([0 0 0]), 1 / D.volume, 'RelTol', 1e-12);
            tc.verifyEqual(D.pdf([10 0 0]), 0);
            tc.verifyEqual(D.volume, 1, 'RelTol', 1e-12);
        end

        function uniformBoxSamplesStayInside(tc)
            D = sfuq.geom.uniformBox([1 2 3], 0.25);
            X = D.sample(2000, RandStream('mt19937ar', 'Seed', 1));
            tc.verifyTrue(all(D.support(X)));
        end

        function gaussianBoxSamplesStayInside(tc)
            D = sfuq.geom.gaussianBox([0 0 0], 0.5, 0.3);
            X = D.sample(2000, RandStream('mt19937ar', 'Seed', 2));
            tc.verifyTrue(all(D.support(X)));
        end

        function gaussianBoxConcentratesAtCentre(tc)
            D = sfuq.geom.gaussianBox([0 0 0], 0.5, 0.15);
            X = D.sample(4000, RandStream('mt19937ar', 'Seed', 3));
            U = sfuq.geom.uniformBox([0 0 0], 0.5);
            Y = U.sample(4000, RandStream('mt19937ar', 'Seed', 4));
            tc.verifyLessThan(mean(vecnorm(X, 2, 2)), mean(vecnorm(Y, 2, 2)), ...
                'Truncated Gaussian must sit closer to the centre than uniform.');
        end

        function likelihoodRatioIsExactForUniformOverUniform(tc)
            % Two uniform boxes, one inside the other. On the overlap the ratio
            % is the volume ratio exactly; outside the smaller box it is Inf.
            P = sfuq.geom.uniformBox([0 0 0], 0.5);   % calibration
            Q = sfuq.geom.uniformBox([0 0 0], 1.0);   % query, larger

            wIn = sfuq.geom.likelihoodRatio(P, Q, [0 0 0]);
            tc.verifyEqual(wIn, P.volume / Q.volume, 'RelTol', 1e-12);

            [wOut, cert] = sfuq.geom.likelihoodRatio(P, Q, [0.8 0 0]);
            tc.verifyTrue(isinf(wOut));
            tc.verifyFalse(cert);
        end

        function likelihoodRatioIsOneForIdenticalDensities(tc)
            P = sfuq.geom.uniformBox([1 1 1], 0.4);
            X = P.sample(200, RandStream('mt19937ar', 'Seed', 5));
            w = sfuq.geom.likelihoodRatio(P, P, X);
            tc.verifyEqual(w, ones(200, 1), 'RelTol', 1e-12, ...
                'No shift must give unit weights.');
        end

        function gaussianOverUniformRatioMatchesClosedForm(tc)
            % w = q/p with q uniform and p an unnormalised truncated Gaussian
            % must be proportional to exp(+d^2 / 2 sigma^2).
            sigma = 0.3;
            P = sfuq.geom.gaussianBox([0 0 0], 0.5, sigma);
            Q = sfuq.geom.uniformBox([0 0 0], 0.5);

            X = [0 0 0; 0.2 0 0; 0.4 0 0];
            w = sfuq.geom.likelihoodRatio(P, Q, X);

            d2 = sum(X.^2, 2);
            expected = exp(0.5 * d2 / sigma^2);
            tc.verifyEqual(w / w(1), expected / expected(1), 'RelTol', 1e-10);
        end

        function certifiableFractionIsOneWhenSupportsMatch(tc)
            P = sfuq.geom.uniformBox([0 0 0], 0.5);
            Q = sfuq.geom.uniformBox([0 0 0], 0.5);
            f = sfuq.design.certifiableFraction(P, Q, 1000, ...
                    RandStream('mt19937ar', 'Seed', 6));
            tc.verifyEqual(f, 1, 'AbsTol', 1e-12);
        end

        function certifiableFractionFallsWhenQueryRegionGrows(tc)
            P = sfuq.geom.uniformBox([0 0 0], 0.5);
            Q = sfuq.geom.uniformBox([0 0 0], 1.0);
            f = sfuq.design.certifiableFraction(P, Q, 4000, ...
                    RandStream('mt19937ar', 'Seed', 7));
            % Volume ratio of the boxes is (0.5/1.0)^3 = 1/8.
            tc.verifyEqual(f, 0.125, 'AbsTol', 0.02);
        end

        function matchedDesignMaximisesEffectiveSampleSize(tc)
            % The theorem behind exp08: p = q gives unit weights and ESS = n.
            Q = sfuq.geom.uniformBox([0 0 0], 0.5);
            X = sfuq.design.matchedDesign(Q, 500, ...
                    RandStream('mt19937ar', 'Seed', 8));
            w = sfuq.geom.likelihoodRatio(Q, Q, X);
            [~, frac] = sfuq.uq.effectiveSampleSize(w);
            tc.verifyEqual(frac, 1, 'AbsTol', 1e-10);

            % A clustered design must do strictly worse.
            P = sfuq.geom.gaussianBox([0 0 0], 0.5, 0.18);
            Xc = P.sample(500, RandStream('mt19937ar', 'Seed', 9));
            wc = sfuq.geom.likelihoodRatio(P, Q, Xc);
            [~, fracC] = sfuq.uq.effectiveSampleSize(wc);
            tc.verifyLessThan(fracC, 0.9);
        end
    end
end
