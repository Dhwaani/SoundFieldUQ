classdef tModels < matlab.unittest.TestCase
%TMODELS  The estimators must be correct before their uncertainty means anything.

    methods (Test)

        function gpInterpolatesNoiselessDataExactly(tc)
            % With vanishing assumed noise the posterior mean must pass
            % through the training observations.
            rng(0);
            S = sfuq.data.makeScenario();
            X = sfuq.data.sampleRegion(40, S.Centre, S.HalfWidth, ...
                                       RandStream('mt19937ar', 'Seed', 5));
            y = sfuq.data.fieldAt(X, S.imgPos, S.imgAmp, S.k);

            M = sfuq.models.fitGP(X, y, S.k, 1e-8);
            mu = sfuq.models.predictGP(M, X);

            tc.verifyEqual(mu, y, 'RelTol', 1e-4, ...
                'GP does not interpolate its own training data.');
        end

        function posteriorVarianceIsNonNegative(tc)
            rng(1);
            S = sfuq.data.makeScenario();
            X = sfuq.data.sampleRegion(60, S.Centre, S.HalfWidth);
            y = sfuq.data.observe(X, S);
            M = sfuq.models.fitGP(X, y, S.k, S.NoiseRel);

            Xq = sfuq.data.sampleRegion(200, S.Centre, 2 * S.HalfWidth);
            [~, v] = sfuq.models.predictGP(M, Xq);
            tc.verifyGreaterThanOrEqual(min(v), 0);
            tc.verifyTrue(all(isfinite(v)));
        end

        function varianceGrowsAwayFromTheData(tc)
            % The property the whole normalized-conformal idea depends on.
            rng(2);
            S = sfuq.data.makeScenario();
            X = sfuq.data.sampleRegion(80, S.Centre, S.HalfWidth, ...
                                       RandStream('mt19937ar', 'Seed', 9));
            y = sfuq.data.observe(X, S, RandStream('mt19937ar', 'Seed', 10));
            M = sfuq.models.fitGP(X, y, S.k, S.NoiseRel);

            near = sfuq.data.sampleRegion(300, S.Centre, 0.3 * S.HalfWidth);
            far  = S.Centre + [3 * S.HalfWidth, 0, 0] + ...
                   sfuq.data.sampleRegion(300, [0 0 0], 0.1);

            [~, vNear] = sfuq.models.predictGP(M, near);
            [~, vFar ] = sfuq.models.predictGP(M, far);

            tc.verifyGreaterThan(mean(vFar), mean(vNear), ...
                'Posterior variance must grow away from the microphones.');
        end

        function krrMeanMatchesGpMean(tc)
            % Kernel ridge with lambda = sigma^2 is the GP posterior mean.
            % Verifying the identity guards against a scaling error in the
            % fitted signal variance leaking into the mean, where it must
            % cancel exactly.
            rng(3);
            S = sfuq.data.makeScenario();
            X = sfuq.data.sampleRegion(50, S.Centre, S.HalfWidth);
            y = sfuq.data.observe(X, S);

            lambda = S.NoiseRel ^ 2;
            Mk = sfuq.models.fitKRR(X, y, S.k, lambda);
            Mg = sfuq.models.fitGP(X,  y, S.k, S.NoiseRel);

            Xq = sfuq.data.sampleRegion(100, S.Centre, S.HalfWidth);
            tc.verifyEqual(sfuq.models.predictKRR(Mk, Xq), ...
                           sfuq.models.predictGP(Mg, Xq), 'RelTol', 1e-8);
        end

        function accuracyImprovesWithMoreMicrophones(tc)
            % Sanity: the estimator must actually estimate.
            S = sfuq.data.makeScenario();
            few  = sfuq.trial(S, 25,  struct('seed', 11, 'nTest', 400, 'nCal', 30));
            many = sfuq.trial(S, 200, struct('seed', 11, 'nTest', 400, 'nCal', 30));
            tc.verifyLessThan(many.nmse, few.nmse - 3, ...
                'NMSE should improve substantially with 8x the microphones.');
        end

        function helmholtzKernelBeatsGenericKernel(tc)
            % The physics has to be worth something, otherwise the framing of
            % the whole project is wrong.
            S = sfuq.data.makeScenario();
            ell = pi / S.k;
            o = struct('seed', 12, 'nTest', 500, 'nCal', 30);
            th = sfuq.trial(S, 80, o);
            o.kernelFcn = @(a, b) sfuq.kernels.gaussian(a, b, ell);
            tg = sfuq.trial(S, 80, o);
            tc.verifyLessThan(th.nmse, tg.nmse, ...
                'Helmholtz prior should outperform a generic RBF kernel.');
        end
    end
end
