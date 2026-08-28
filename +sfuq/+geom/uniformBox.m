function D = uniformBox(centre, halfWidth)
%UNIFORMBOX  Uniform spatial density over an axis-aligned box.
%
%   D = UNIFORMBOX(CENTRE, HALFWIDTH) returns a spatial-density struct with
%   fields:
%     .sample(n, stream)  draw n points, N-by-3
%     .pdf(X)             density at each row of X (0 outside the box)
%     .support(X)         logical, true inside the box
%     .name, .centre, .halfWidth, .volume
%
%   Why densities are first-class objects here
%   ------------------------------------------
%   The contribution of this project rests on the likelihood ratio between the
%   calibration and query distributions being EXACTLY known rather than
%   estimated. That is only true because both distributions are *designed* --
%   the experimenter chooses where to put calibration microphones and which
%   listener positions to certify. Representing each as an explicit density
%   object, with an analytic pdf, is what makes the ratio exact and auditable
%   instead of an approximation buried in a fitting step.
%
%   See also SFUQ.GEOM.GAUSSIANBOX, SFUQ.GEOM.LIKELIHOODRATIO.

    centre = centre(:).';
    if isscalar(halfWidth)
        halfWidth = repmat(halfWidth, 1, 3);
    end
    halfWidth = halfWidth(:).';

    vol = prod(2 * halfWidth);

    D = struct();
    D.name      = 'uniformBox';
    D.centre    = centre;
    D.halfWidth = halfWidth;
    D.volume    = vol;

    D.support = @(X) all(abs(X - centre) <= halfWidth + 1e-12, 2);
    D.pdf     = @(X) double(D.support(X)) / vol;
    D.sample  = @(n, st) sampleBox(n, centre, halfWidth, st);
end

function X = sampleBox(n, centre, halfWidth, st)
    if nargin < 4 || isempty(st)
        u = rand(n, 3);
    else
        u = rand(st, n, 3);
    end
    X = centre + (2 * u - 1) .* halfWidth;
end
