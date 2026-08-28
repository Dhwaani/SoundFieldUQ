function [pos, amp] = imageSources(src, room, beta, order)
%IMAGESOURCES  Mirror-image source positions and amplitudes for a shoebox room.
%
%   [POS, AMP] = IMAGESOURCES(SRC, ROOM, BETA, ORDER) returns the image-source
%   expansion of a point source at SRC in a rigid-walled rectangular room.
%
%   Inputs
%     src    1-by-3 source position [x y z] in metres
%     room   1-by-3 room dimensions [Lx Ly Lz] in metres
%     beta   scalar wall reflection coefficient in [0, 1)
%     order  maximum reflection order
%
%   Outputs
%     pos    P-by-3 image positions
%     amp    P-by-1 reflection amplitudes, beta^(number of reflections)
%
%   This is the Allen-Berkley construction. The eight sign/parity
%   combinations generate the mirrored positions within one "tile", and the
%   integer offsets tile the mirrored rooms outward.
%
%   Using image sources rather than a measured impulse response is what makes
%   this project hardware independent: the ground-truth field is available
%   analytically at ANY point, so held-out test positions are exact rather
%   than interpolated, and coverage can be measured against truth instead of
%   against another estimate.

    validateattributes(src,   {'numeric'}, {'numel', 3, 'finite', 'real'});
    validateattributes(room,  {'numeric'}, {'numel', 3, 'positive', 'finite'});
    validateattributes(beta,  {'numeric'}, {'scalar', '>=', 0, '<', 1});
    validateattributes(order, {'numeric'}, {'scalar', 'integer', 'nonnegative'});

    src  = src(:).';
    room = room(:).';

    n = -order:order;
    [NX, NY, NZ, PX, PY, PZ] = ndgrid(n, n, n, [0 1], [0 1], [0 1]);
    NX = NX(:); NY = NY(:); NZ = NZ(:);
    PX = PX(:); PY = PY(:); PZ = PZ(:);

    % Mirrored coordinates: parity p flips the source about the near wall,
    % the integer offset n translates by whole room widths.
    x = (1 - 2 * PX) * src(1) + 2 * NX * room(1);
    y = (1 - 2 * PY) * src(2) + 2 * NY * room(2);
    z = (1 - 2 * PZ) * src(3) + 2 * NZ * room(3);

    nRefl = abs(2 * NX - PX) + abs(PX) + ...
            abs(2 * NY - PY) + abs(PY) + ...
            abs(2 * NZ - PZ) + abs(PZ);

    keep = nRefl <= 3 * order;
    pos  = [x(keep), y(keep), z(keep)];
    amp  = beta .^ nRefl(keep);
end
