classdef tKernel < matlab.unittest.TestCase
%TKERNEL  Properties the Helmholtz kernel must satisfy.
%
%   The kernel is the physics in this project. If it is not positive
%   semi-definite the Cholesky factorisation fails; if its value at zero
%   distance is wrong every variance is wrong. These are cheap tests that
%   catch expensive mistakes.

    methods (Test)

        function diagonalIsUnity(tc)
            X = rand(20, 3);
            K = sfuq.kernels.helmholtz(X, X, 10);
            tc.verifyEqual(diag(K), ones(20, 1), 'AbsTol', 1e-12, ...
                'k(r,r) must equal 1 for the normalised Helmholtz kernel.');
        end

        function isSymmetric(tc)
            X = rand(30, 3);
            K = sfuq.kernels.helmholtz(X, X, 7.5);
            tc.verifyEqual(K, K.', 'AbsTol', 1e-12);
        end

        function isPositiveSemiDefinite(tc)
            % Checked across a range of wavenumbers: the kernel becomes
            % increasingly oscillatory with k, and a sign error in the series
            % expansion would show up as a negative eigenvalue at high k.
            for k = [1 5 15 40]
                X = rand(120, 3) - 0.5;
                K = sfuq.kernels.helmholtz(X, X, k);
                ev = eig((K + K.') / 2);
                tc.verifyGreaterThan(min(ev), -1e-8, ...
                    sprintf('Kernel not PSD at k = %g.', k));
            end
        end

        function matchesClosedFormSinc(tc)
            a = [0 0 0];
            b = [0.3 0 0];
            k = 12;
            x = k * 0.3;
            tc.verifyEqual(sfuq.kernels.helmholtz(a, b, k), sin(x) / x, ...
                'RelTol', 1e-12);
        end

        function handlesTinyDistanceWithoutNaN(tc)
            % sin(x)/x is 0/0 at the origin. A naive implementation returns
            % NaN, which then propagates through every downstream matrix.
            a = [1 1 1];
            b = a + 1e-15;
            K = sfuq.kernels.helmholtz(a, b, 20);
            tc.verifyFalse(isnan(K));
            tc.verifyEqual(K, 1, 'AbsTol', 1e-9);
        end

        function distanceMatrixNeverComplex(tc)
            % The |a|^2 - 2ab + |b|^2 expansion can go microscopically
            % negative for coincident points; sqrt of that is complex and
            % would poison everything silently.
            X = repmat([0.5 0.5 0.5], 50, 1);
            D = sfuq.kernels.pairwiseDist(X, X);
            tc.verifyTrue(isreal(D));
            tc.verifyEqual(D, zeros(50), 'AbsTol', 1e-12);
        end
    end
end
