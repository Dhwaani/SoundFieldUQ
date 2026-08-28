function T = shiftTrial(S, nMic, Pcal, Qquery, opts)
%SHIFTTRIAL  One fit / calibrate / certify cycle under a designed covariate shift.
%
%   T = SHIFTTRIAL(S, NMIC, PCAL, QQUERY) fits the Helmholtz GP to NMIC
%   microphones, draws calibration points from the density PCAL and query
%   points from the density QQUERY, and returns everything the conformal
%   methods need -- including the exact likelihood ratios.
%
%   The distinction from SFUQ.TRIAL
%   -------------------------------
%   SFUQ.TRIAL draws calibration and test points from the SAME distribution, so
%   exchangeability holds and split conformal is valid. This function
%   deliberately breaks that, in the way a real deployment breaks it: the
%   calibration measurements sit where the microphones are, and the positions
%   to be certified are where the listeners are. Those are different
%   distributions, both designed by the experimenter, so the likelihood ratio
%   is exact.
%
%   Options
%     .nCal    calibration points (default 400)
%     .nQuery  query points       (default 1500)
%     .seed    trial seed         (default S.Seed)
%
%   Returned fields
%     .resCal, .resQuery    absolute residuals
%     .sigCal, .sigQuery    GP predictive standard deviations
%     .wCal,  .wQuery       EXACT likelihood ratios dQ/dP
%     .certifiable          logical, query points with finite weight
%     .Xmic, .Xcal, .Xquery positions
%     .nmse                 reconstruction NMSE in dB (certifiable points only)

    if nargin < 5, opts = struct(); end
    if ~isfield(opts, 'nCal'),   opts.nCal = 400;    end
    if ~isfield(opts, 'nQuery'), opts.nQuery = 1500; end
    if ~isfield(opts, 'seed'),   opts.seed = S.Seed; end

    st = RandStream('mt19937ar', 'Seed', opts.seed);

    Xmic   = sfuq.data.sampleRegion(nMic, S.Centre, S.HalfWidth, st);
    Xcal   = Pcal.sample(opts.nCal, st);
    Xquery = Qquery.sample(opts.nQuery, st);

    yMic              = sfuq.data.observe(Xmic,   S, st);
    yCal              = sfuq.data.observe(Xcal,   S, st);
    [yQuery, pQuery]  = sfuq.data.observe(Xquery, S, st);

    M = sfuq.models.fitGP(Xmic, yMic, S.k, S.NoiseRel);

    [muCal,   varCal]   = sfuq.models.predictGP(M, Xcal);
    [muQuery, varQuery] = sfuq.models.predictGP(M, Xquery);

    T = struct();
    T.Xmic   = Xmic;
    T.Xcal   = Xcal;
    T.Xquery = Xquery;
    T.model  = M;

    T.resCal   = abs(yCal   - muCal);
    T.resQuery = abs(yQuery - muQuery);
    T.sigCal   = sqrt(varCal);
    T.sigQuery = sqrt(varQuery);

    % The exact weights. Calibration weights must be finite by construction --
    % a calibration point always lies in the support of its own density.
    T.wCal   = sfuq.geom.likelihoodRatio(Pcal, Qquery, Xcal);
    [T.wQuery, T.certifiable] = sfuq.geom.likelihoodRatio(Pcal, Qquery, Xquery);

    % Calibration points outside the query support carry zero weight: they are
    % informative about nowhere anyone is listening.
    T.wCal(~isfinite(T.wCal)) = 0;

    T.dQuery = sqrt(sum((Xquery - S.Centre).^2, 2));

    if any(T.certifiable)
        T.nmse = sfuq.eval.nmse(muQuery(T.certifiable), pQuery(T.certifiable));
    else
        T.nmse = NaN;
    end
end
