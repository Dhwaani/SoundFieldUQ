function p = fieldAt(pts, pos, amp, k)
%FIELDAT  Complex sound pressure at arbitrary points, in the frequency domain.
%
%   P = FIELDAT(PTS, POS, AMP, K) evaluates
%
%       p(r) = sum_i  amp_i * exp(-1i*k*||r - pos_i||) / (4*pi*||r - pos_i||)
%
%   the superposition of monopole (free-field Green's function) contributions
%   from every image source.
%
%   Inputs
%     pts  N-by-3 evaluation points
%     pos  P-by-3 image source positions (from SFUQ.DATA.IMAGESOURCES)
%     amp  P-by-1 image source amplitudes
%     k    wavenumber, 2*pi*f/c
%
%   Output
%     p    N-by-1 complex pressure
%
%   Working in the frequency domain at a single wavenumber, rather than
%   synthesising a time-domain impulse response, avoids fractional-delay
%   interpolation entirely: this expression is an exact solution of the
%   Helmholtz equation everywhere except at the source points themselves.
%   The unit test verifies that numerically to about 1e-5 relative residual.
%
%   A time-domain room impulse response, when needed, is obtained by
%   evaluating this over a grid of wavenumbers and inverse-transforming.

    D = sfuq.kernels.pairwiseDist(pts, pos);

    % Guard the singularity at a source point. A microphone placed exactly on
    % an image source is unphysical, and without this the field is Inf.
    D = max(D, 1e-6);

    G = exp(-1i * k * D) ./ (4 * pi * D);
    p = G * amp(:);
end
