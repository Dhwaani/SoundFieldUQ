function K = gaussian(A, B, ell)
%GAUSSIAN  Squared-exponential kernel, used only as a physics-free control.
%
%   K = GAUSSIAN(A, B, ELL) with length scale ELL.
%
%   Included so the experiments can quantify how much the Helmholtz prior is
%   actually worth: the same estimator run with this kernel is the
%   "no physics" ablation. It should need substantially more microphones to
%   reach the same NMSE, and that gap is one of the paper's results.

    D = sfuq.kernels.pairwiseDist(A, B);
    K = exp(-0.5 * (D / ell) .^ 2);
end
