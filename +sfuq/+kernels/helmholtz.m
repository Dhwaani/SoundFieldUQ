function K = helmholtz(A, B, k)
%HELMHOLTZ  Reproducing kernel for interior sound fields.
%
%   K = HELMHOLTZ(A, B, K_WAVENUMBER) returns the kernel matrix with
%   K(i,j) = j0(k * ||A(i,:) - B(j,:)||), where j0 is the zeroth-order
%   spherical Bessel function, j0(x) = sin(x)/x.
%
%   Why this kernel
%   ---------------
%   Any sound field that satisfies the homogeneous Helmholtz equation in a
%   source-free region can be written as a superposition of plane waves whose
%   wavevectors lie on the sphere of radius k. Placing a uniform measure on
%   that sphere and integrating exp(-i k'*(r - r')) over it gives exactly
%   sin(k||r-r'||) / (k||r-r'||). The kernel therefore encodes the wave
%   physics: functions in its RKHS satisfy the Helmholtz equation by
%   construction, which is what makes it a far stronger prior than a generic
%   squared-exponential kernel for the same number of measurements.
%
%   The kernel is real and positive semi-definite even though the sound field
%   it models is complex; the complex field is handled by solving the real
%   linear system with complex right-hand sides.
%
%   Note K(r,r) = 1 for all r, so the kernel is already normalised; the signal
%   power is carried by a separate scale factor fitted in sfuq.models.fitGP.
%
%   See also SFUQ.KERNELS.PAIRWISEDIST, SFUQ.MODELS.FITGP.

    if ~isscalar(k) || ~isreal(k) || k <= 0
        error('sfuq:helmholtz:badWavenumber', ...
              'Wavenumber k must be a positive real scalar.');
    end

    x = k * sfuq.kernels.pairwiseDist(A, B);

    % sin(x)/x has a removable singularity at x = 0 where it equals 1.
    % Use a series-safe threshold rather than testing x == 0 exactly: for very
    % small but non-zero x the naive quotient loses precision.
    K = ones(size(x));
    big = x > 1e-8;
    K(big) = sin(x(big)) ./ x(big);
end
