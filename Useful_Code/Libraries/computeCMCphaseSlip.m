function obsSolutions = computeCMCphaseSlip(obsSolutions, opts)
% COMPUTECMCPHASESLIP  Phase-slip flagging from the Code-minus-Carrier (CMC)
%                      time series.
%
% METHOD
%   For each observable channel:
%     1. Build the Code-minus-Carrier combination  cmc = P - L  [m], which is
%        geometry- and clock-free, so that it retains only the (constant)
%        phase ambiguity, twice the ionospheric delay, the code/phase hardware
%        biases, multipath and noise.
%     2. Differentiate cmc w.r.t. receiver time and flag outliers of
%        d(cmc)/dt through Grubbs' test. A cycle slip appears as an impulsive
%        discontinuity in this rate, hence as a statistical outlier.
%
% DIFFERENCING CONVENTION
%   A BACKWARD difference labelled at the current epoch is used,
%
%       d(k) = [ cmc(k) - cmc(k-1) ] / [ t(k) - t(k-1) ]              [m/s]
%
%   so that a slip whose first corrupted sample is k produces a single
%   impulse exactly at index k. Central differencing must NOT be used here:
%   a step of amplitude D at epoch k appears in the central-difference series
%   as TWO impulses of amplitude D/(2*dt) at indices k-1 and k, so the flag is
%   smeared over two samples and is not centred on the event. The per-sample
%   signal-to-noise ratio is the same in both conventions (0.707*D/sigma),
%   because central differencing halves the noise as well as the step, so
%   nothing is lost in sensitivity by concentrating the step on one sample.
%
%   The resulting <bool_phase_slip> is therefore a POINT flag: exactly one
%   true sample per detected event, sitting on the first corrupted epoch.
%   Downstream masks are expected to dilate it. A +/-1 dilation is the natural
%   choice for ROT exclusion, because <computeROT.m> itself uses
%   <finite_diff_central.m>, so a slip at k corrupts ROT(k-1) and ROT(k); the
%   dilated set {k-1, k, k+1} covers both with one sample of margin. Set
%   <Dilate> to have that mask returned alongside the point flag.
%
% INPUTS
%   obsSolutions  [1xNch] struct array of raw observables. Required fields:
%                   .rawP          [1xN] code pseudorange                  [m]
%                   .CarrierPhase  [1xN] carrier phase                     [m]
%                   .rxTime        [1xN] receiver time, strictly increasing [s]
%   opts          (optional) name-value options:
%                   .OutlierMethod    method passed to ISOUTLIER
%                                                        (default "grubbs")
%                   .CollapseClusters keep only the leading sample of each run
%                                     of consecutive flags    (default true)
%                   .Dilate           half-width of the additional dilated
%                                     mask, in samples; 0 disables it
%                                                               (default 0)
%                   .StoreIndices     also return the CMC rate series
%                                                            (default true)
%
% OUTPUT
%   obsSolutions  same struct array, augmented with:
%                   .bool_phase_slip [1xN] logical, point flag on the first
%                                    corrupted epoch of each detected slip
%                 and, when the corresponding option is active:
%                   .bool_phase_slip_mask [1xN] logical, <Dilate>-dilated flag
%                   .dcmc_dtime           [1xN] d(cmc)/dt                [m/s]
%
% NOTE
%   Grubbs' test assumes normality and detects one outlier at a time; ISOUTLIER
%   iterates it until no further detection occurs. On long series with closely
%   spaced slips, masking remains possible. A robust alternative is
%   OutlierMethod = "movmedian" (MAD-based, scale-adaptive).
%
%   Grubbs is also a GLOBAL test: it compares every sample against the mean and
%   standard deviation of the whole series. An occultation arc is far from
%   stationary, so the threshold ends up dominated by the disturbed portion and
%   sensitivity drops in the quiet portion. OutlierMethod = "movmedian" is the
%   cheapest way to make the test scale-adaptive without changing anything else.
%
%   Because the statistic is a RATE, a slip occurring across a long data hole
%   is divided by a large dt and may fall below the detection threshold. That
%   case is already covered downstream by the <t_gap> re-initialisation of
%   <STECleveling.m>.
%
%
% DEPENDENCIES
%   none
%
% VERSION
%   v1.0 - Grubbs' test on the central-difference CMC rate.
%   v1.1 - Backward difference so the flag lands exactly on the first
%          corrupted epoch; consecutive flags collapsed to a point flag;
%          optional dilated mask; output orientation made [1xN] as documented.
%          Detection statistic and outlier test unchanged.
%   v1.2 - Difference taken against the previous VALID sample, removing the
%          phantom rate that central differencing produced at missing epochs
%          and recovering the sample after an isolated gap.
%
%   Behaviour change w.r.t. v1.0 worth checking: v1.0 returned a first-order
%   forward difference at index 1 and a backward difference at index N, so
%   sample 1 carried a value and could be flagged. It is now NaN, since it has
%   no predecessor to be compared against. Sample 1 is re-initialised by
%   <STECleveing.m> at <firstOccurrence> regardless.
%
% -------------------------------------------------------------------------
% NavSAS Group - Dipartimento di Elettronica e Telecomunicazioni
% Politecnico di Torino
% Francesco Fiorina
% -------------------------------------------------------------------------

arguments
    obsSolutions                      struct
    opts.OutlierMethod    (1,1) string  = "grubbs"
    opts.CollapseClusters (1,1) logical = true
    opts.Dilate           (1,1) double  = 0
    opts.StoreIndices     (1,1) logical = true
end

% Guard against being fed a RINEX-style dataset, where the carrier phase is
% stored in <rawCarrP> in cycles and would need a wavelength scaling first.
if ~isfield(obsSolutions,'CarrierPhase')
    error(['Field "CarrierPhase" [m] not found in <obsSolutions>. ' ...
           'If the carrier phase is stored as "rawCarrP" [cycles], scale it ' ...
           'by the channel wavelength before calling this function.']);
end

% --- Outliers detection on rate of change of code-minus-carrier (CMC) time-series

for ii = 1:numel(obsSolutions)

    % observation time series - both observables already expressed in [m]
    t             = double(obsSolutions(ii).rxTime(:));
    code_vs_time  = double(obsSolutions(ii).rawP(:));
    phase_vs_time = double(obsSolutions(ii).CarrierPhase(:));
    N             = numel(t);

    cmc = code_vs_time - phase_vs_time;

    % Backward difference labelled at the current epoch: the whole step falls
    % on the first corrupted sample instead of being split between k-1 and k.
    %
    % The difference is taken against the previous VALID sample rather than
    % the previous index, with the true elapsed time as denominator. An
    % isolated missing epoch then costs a single rate sample - the missing one
    % itself, which is untestable anyway - instead of two, and a gap is never
    % bridged without its dt being accounted for. Non-monotonic time steps are
    % invalidated rather than producing a sign flip or a division by zero.
    v          = find(isfinite(cmc) & isfinite(t));
    dcmc_dtime = nan(N,1);                                            % [m/s]
    if numel(v) >= 2
        dtv           = diff(t(v));
        dtv(dtv <= 0) = NaN;
        dcmc_dtime(v(2:end)) = diff(cmc(v)) ./ dtv;
    end

    % perform outlier detection on time difference of Code-minus-Carrier (CMC).
    % NaN samples (gaps in phase tracking, and the first epoch, which has no
    % predecessor) are treated as missing data by ISOUTLIER and never flagged.
    outliers = isoutlier(dcmc_dtime, opts.OutlierMethod);
    outliers(~isfinite(dcmc_dtime)) = false;

    % Collapse each run of consecutive flags to its leading sample, so that
    % every detected event contributes exactly one true value and a symmetric
    % dilation stays centred on the event.
    if opts.CollapseClusters
        outliers = outliers & ~[false; outliers(1:end-1)];
    end

    % --- Store in "obsSolutions"
    obsSolutions(ii).bool_phase_slip = reshape(outliers, 1, N);   % point flag
                                        % on the first corrupted epoch
%     Consider that dilatation is already implemented in the ROTI
%     computation function.
%     if opts.Dilate > 0
%         mask = movmax(double(outliers), [opts.Dilate opts.Dilate]) > 0;
%         obsSolutions(ii).bool_phase_slip_mask = reshape(mask, 1, N);
%     end

    if opts.StoreIndices
        obsSolutions(ii).dcmc_dtime = reshape(dcmc_dtime, 1, N);
    end

end

end % computeCMCphaseSlip