function X = sampleRegion(n, centre, halfWidth, rngStream)
%SAMPLEREGION  Uniform random points in an axis-aligned cuboid.
%
%   X = SAMPLEREGION(N, CENTRE, HALFWIDTH) returns N-by-3 points drawn
%   uniformly from the cuboid centred at CENTRE with the given half-width
%   (scalar, or 1-by-3 for an anisotropic region).
%
%   X = SAMPLEREGION(N, CENTRE, HALFWIDTH, S) draws from the RandStream S,
%   which is how the experiments keep the microphone layout, the calibration
%   set and the noise realisation independently reproducible.

    if nargin < 4 || isempty(rngStream)
        u = rand(n, 3);
    else
        u = rand(rngStream, n, 3);
    end

    centre = centre(:).';
    if isscalar(halfWidth)
        halfWidth = repmat(halfWidth, 1, 3);
    end
    halfWidth = halfWidth(:).';

    X = centre + (2 * u - 1) .* halfWidth;
end
