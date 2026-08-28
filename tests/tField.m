classdef tField < matlab.unittest.TestCase
%TFIELD  The synthetic field must be a genuine solution of the physics.
%
%   Everything in this project is measured against this field as ground
%   truth. If it does not satisfy the Helmholtz equation, then the kernel's
%   prior is mismatched to the data by construction and every calibration
%   conclusion is an artefact of that mismatch rather than a property of the
%   estimators.

    methods (Test)

        function satisfiesHelmholtzEquation(tc)
            % Central finite differences of the Laplacian, checked against
            % -k^2 p at an interior point well away from any image source.
            S = sfuq.data.makeScenario('Frequency', 500);
            c0 = S.Centre;
            h = 4e-3;

            pts = c0;
            for ax = 1:3
                for s = [1 -1]
                    d = zeros(1, 3);
                    d(ax) = s * h;
                    pts(end + 1, :) = c0 + d; %#ok<AGROW>
                end
            end

            v = sfuq.data.fieldAt(pts, S.imgPos, S.imgAmp, S.k);

            lap = 0;
            for ax = 1:3
                lap = lap + v(2 * ax) + v(2 * ax + 1) - 2 * v(1);
            end
            lap = lap / h^2;

            residual = abs(lap + S.k^2 * v(1)) / (S.k^2 * abs(v(1)));
            tc.verifyLessThan(residual, 1e-2, ...
                'Synthetic field does not satisfy the Helmholtz equation.');
        end

        function obeysInverseSquareLawInFreeField(tc)
            % A single source with no reflections must fall off as 1/r.
            pos = [0 0 0];
            amp = 1;
            k = 10;
            r1 = sfuq.data.fieldAt([1 0 0], pos, amp, k);
            r2 = sfuq.data.fieldAt([2 0 0], pos, amp, k);
            tc.verifyEqual(abs(r1) / abs(r2), 2, 'RelTol', 1e-9);
        end

        function imageSourceCountGrowsWithOrder(tc)
            [p0, ~] = sfuq.data.imageSources([1 1 1], [5 4 3], 0.5, 0);
            [p2, ~] = sfuq.data.imageSources([1 1 1], [5 4 3], 0.5, 2);
            tc.verifyGreaterThan(size(p2, 1), size(p0, 1));
        end

        function reflectionAmplitudesDecay(tc)
            [~, amp] = sfuq.data.imageSources([1 1 1], [5 4 3], 0.5, 3);
            tc.verifyLessThanOrEqual(max(amp), 1);
            tc.verifyGreaterThan(numel(unique(amp)), 1, ...
                'All amplitudes identical: reflection order is not being counted.');
        end

        function observationsAreReproducible(tc)
            S = sfuq.data.makeScenario();
            X = sfuq.data.sampleRegion(10, S.Centre, S.HalfWidth, ...
                                       RandStream('mt19937ar', 'Seed', 1));
            y1 = sfuq.data.observe(X, S, RandStream('mt19937ar', 'Seed', 7));
            y2 = sfuq.data.observe(X, S, RandStream('mt19937ar', 'Seed', 7));
            tc.verifyEqual(y1, y2, ...
                'Same seed must give identical observations.');
        end
    end
end
