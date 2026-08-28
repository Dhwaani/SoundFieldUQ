function D = pairwiseDist(A, B)
%PAIRWISEDIST  Euclidean distance matrix between two sets of points.
%
%   D = PAIRWISEDIST(A, B) returns the N-by-M matrix with D(i,j) equal to the
%   Euclidean distance between A(i,:) and B(j,:).
%
%   Inputs
%     A  N-by-d matrix of points (one point per row)
%     B  M-by-d matrix of points
%
%   Written out by hand rather than calling PDIST2 so the whole project runs
%   on base MATLAB with no toolbox dependency. The expansion
%   |a-b|^2 = |a|^2 - 2 a.b + |b|^2 is used, and the result is clamped at zero
%   before the square root: catastrophic cancellation can make the expansion
%   very slightly negative for coincident points, and sqrt of -1e-17 is
%   complex, which would silently poison every downstream kernel.

    if size(A, 2) ~= size(B, 2)
        error('sfuq:pairwiseDist:dimMismatch', ...
              'A and B must have the same number of columns.');
    end

    aa = sum(A .^ 2, 2);            % N-by-1
    bb = sum(B .^ 2, 2).';          % 1-by-M
    D2 = aa + bb - 2 * (A * B.');
    D  = sqrt(max(D2, 0));
end
