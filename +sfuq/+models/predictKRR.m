function mu = predictKRR(M, Xq)
%PREDICTKRR  Evaluate a fitted kernel ridge sound-field model.
%
%   MU = PREDICTKRR(M, XQ) returns the complex pressure estimate at XQ.

    mu = M.kernelFcn(Xq, M.X) * M.alpha;
end
