function D = gaussianBox(centre, halfWidth, sigma)
%GAUSSIANBOX  Gaussian density truncated to an axis-aligned box.
%
%   D = GAUSSIANBOX(CENTRE, HALFWIDTH, SIGMA) models the realistic case where
%   calibration microphones CLUSTER near the centre of an array rather than
%   covering the listening area uniformly.
%
%   The truncation is the whole point
%   ---------------------------------
%   Truncating to the box keeps the SUPPORT identical to a uniform density
%   over the same box. That matters more than it may appear: weighted conformal
%   prediction requires the query distribution to be absolutely continuous with
%   respect to the calibration distribution. Two densities on the same support
%   satisfy this and the likelihood ratio stays finite everywhere, so coverage
%   can be restored. Two densities on DIFFERENT supports do not, and no
%   reweighting can repair it -- see SFUQ.GEOM.LIKELIHOODRATIO and
%   EXP07_SUPPORT_FAILURE.
%
%   The normalising constant of the truncated Gaussian cancels in the ratio
%   against another density on the same box, so it is never needed and is not
%   computed; PDF therefore returns an UNNORMALISED density. This is flagged in
%   the struct as .normalised = false so nothing downstream mistakes it for a
%   probability density.
%
%   SIGMA may be a scalar or a 1-by-3 vector for an anisotropic cluster.

    centre = centre(:).';
    if isscalar(halfWidth), halfWidth = repmat(halfWidth, 1, 3); end
    if isscalar(sigma),     sigma     = repmat(sigma, 1, 3);     end
    halfWidth = halfWidth(:).';
    sigma     = sigma(:).';

    D = struct();
    D.name       = 'gaussianBox';
    D.centre     = centre;
    D.halfWidth  = halfWidth;
    D.sigma      = sigma;
    D.normalised = false;

    D.support = @(X) all(abs(X - centre) <= halfWidth + 1e-12, 2);
    D.pdf     = @(X) double(D.support(X)) .* ...
                     exp(-0.5 * sum(((X - centre) ./ sigma).^2, 2));
    D.sample  = @(n, st) rejectionSample(n, centre, halfWidth, sigma, st);
end

function X = rejectionSample(n, centre, halfWidth, sigma, st)
%   Rejection sampling: draw from the untruncated Gaussian, keep what lands
%   inside the box. Simple and exact. The acceptance rate is checked so a
%   pathological sigma/halfWidth combination fails loudly rather than spinning.

    X = zeros(0, 3);
    attempts = 0;
    while size(X, 1) < n
        m = max(2 * (n - size(X, 1)), 256);
        if nargin < 5 || isempty(st)
            C = centre + randn(m, 3) .* sigma;
        else
            C = centre + randn(st, m, 3) .* sigma;
        end
        ok = all(abs(C - centre) <= halfWidth, 2);
        X = [X; C(ok, :)]; %#ok<AGROW>

        attempts = attempts + 1;
        if attempts > 200
            error('sfuq:gaussianBox:lowAcceptance', ...
                  ['Rejection sampling is not converging. sigma is probably ' ...
                   'far larger than halfWidth; use uniformBox instead.']);
        end
    end
    X = X(1:n, :);
end
