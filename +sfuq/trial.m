function T = trial(S, nMic, opts)
%TRIAL  One complete fit / calibrate / test cycle for a scenario.
%
%   T = TRIAL(S, NMIC) draws a microphone layout, a calibration set and a test
%   set for scenario S, fits the Helmholtz GP, and returns everything the
%   uncertainty analyses need.
%
%   T = TRIAL(S, NMIC, OPTS) accepts a struct with fields:
%     .nCal        calibration points                      (default 300)
%     .nTest       test points                             (default 800)
%     .testHalf    half-width of the TEST region           (default S.HalfWidth)
%     .seed        trial seed                              (default S.Seed)
%     .kernelFcn   custom kernel                           (default Helmholtz)
%     .sigmaScale  multiplies the assumed noise level      (default 1)
%
%   The three-way split is the essential discipline
%   -----------------------------------------------
%   Training, calibration and test sets are disjoint and independently drawn.
%   Split conformal's guarantee holds only if the calibration residuals were
%   NOT used to fit the model; reusing training points would make the
%   residuals optimistically small and the radii too narrow, silently voiding
%   the guarantee while still producing plausible-looking numbers.
%
%   Setting opts.testHalf larger than S.HalfWidth moves test points outside
%   the microphone aperture while the calibration set stays inside it. That
%   deliberately breaks the exchangeability assumption, and is how experiment
%   04 probes what happens when conformal prediction is used outside the
%   regime it is licensed for.
%
%   Returned struct T
%     .Xmic .Xcal .Xtest   positions
%     .resCal .resTest     absolute residuals against the NOISY observations
%     .resTestTrue         absolute residuals against noiseless ground truth
%     .sigCal .sigTest     GP predictive standard deviations
%     .muTest .pTest       predicted and true complex pressure at test points
%     .dTest               distance of each test point from the region centre
%     .nmse                reconstruction NMSE in dB

    if nargin < 3
        opts = struct();
    end
    if ~isfield(opts, 'nCal'),       opts.nCal = 300;              end
    if ~isfield(opts, 'nTest'),      opts.nTest = 800;             end
    if ~isfield(opts, 'testHalf'),   opts.testHalf = S.HalfWidth;  end
    if ~isfield(opts, 'seed'),       opts.seed = S.Seed;           end
    if ~isfield(opts, 'kernelFcn'),  opts.kernelFcn = [];          end
    if ~isfield(opts, 'sigmaScale'), opts.sigmaScale = 1;          end

    st = RandStream('mt19937ar', 'Seed', opts.seed);

    Xmic  = sfuq.data.sampleRegion(nMic,       S.Centre, S.HalfWidth,    st);
    Xcal  = sfuq.data.sampleRegion(opts.nCal,  S.Centre, S.HalfWidth,    st);
    Xtest = sfuq.data.sampleRegion(opts.nTest, S.Centre, opts.testHalf,  st);

    [yMic,  ~]      = sfuq.data.observe(Xmic,  S, st);
    [yCal,  ~]      = sfuq.data.observe(Xcal,  S, st);
    [yTest, pTest]  = sfuq.data.observe(Xtest, S, st);

    M = sfuq.models.fitGP(Xmic, yMic, S.k, ...
                          S.NoiseRel * opts.sigmaScale, opts.kernelFcn);

    [muCal,  varCal ] = sfuq.models.predictGP(M, Xcal);
    [muTest, varTest] = sfuq.models.predictGP(M, Xtest);

    T = struct();
    T.Xmic  = Xmic;
    T.Xcal  = Xcal;
    T.Xtest = Xtest;

    T.resCal      = abs(yCal  - muCal);
    T.resTest     = abs(yTest - muTest);
    T.resTestTrue = abs(pTest - muTest);

    T.sigCal  = sqrt(varCal);
    T.sigTest = sqrt(varTest);

    T.muTest = muTest;
    T.pTest  = pTest;
    T.model  = M;

    T.dTest = sqrt(sum((Xtest - S.Centre) .^ 2, 2));
    T.nmse  = sfuq.eval.nmse(muTest, pTest);
end
