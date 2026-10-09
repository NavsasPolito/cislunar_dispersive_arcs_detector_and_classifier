function ROTI = rotiArcRestore(ROTI, obsSolutions, PRNid1)
%ROTIARCRRESTORE  Fills missing ROTI samples within arc-masked intervals.
%
% DESCRIPTION
% -------------------------------------------------------------------------
% This function restores missing ROTI observations within tracked arcs by
% interpolating over gaps using only the samples that belong to the expected
% arc mask. Samples outside the arc mask are left as NaN and are NOT filled.
%
% The intended use is to recover ROTI values to have contiguous arcs
%
% The function operates as follows:
%
%   1) A temporary ROTI vector (ROTI_tmp) is initialized as NaN.
%
%   2) Makima interpolation is performed using only samples that belong to
%      the arc mask (arcMask) AND are currently valid (non-NaN ROTI).
%      The interpolated values are evaluated at ALL mask positions.
%
%   3) Only the positions that:
%        - belong to the arc mask  (arcMask == true)
%        - are currently missing   (isnan(ROTI) == true)
%      receive the interpolated value. Existing valid samples are preserved.
%
% -------------------------------------------------------------------------
% INPUT
%
%   ROTI          [Nx1] double
%       Rate Of TEC Index vector. May contain NaN where ROTI could not be
%       computed (gaps shorter than the ROTI window, edge removal, etc.).
%
%   obsSolutions  struct
%       Observation solution structure for the satellite. Must contain:
%           .rxTime   [Nx1] receiver time vector (seconds)
%           .arcMask  [Nx1] logical mask of expected valid arc samples
%
%   PRNid1        integer
%       Index into obsSolutions corresponding to the primary frequency
%       of the satellite being processed.
%
% -------------------------------------------------------------------------
% OUTPUT
%
%   ROTI          [Nx1] double
%       ROTI vector with missing-but-masked samples filled by interpolation.
%       Samples outside the arc mask are untouched and remain NaN.
%
% -------------------------------------------------------------------------
% NOTES
%
%   - Makima interpolation is chosen for its robustness against overshoots
%     at arc boundaries compared to cubic spline interpolation.
%   - The function only FILLS gaps; it never overwrites existing valid data.
%   - Sufficient valid points must be present within the arc for interp1
%     to produce meaningful results (at least 2 non-NaN samples).
%   - This function must be called AFTER all segments of a satellite have
%     been processed, so that the full ROTI vector is available as anchor
%     data for interpolation.
%
% -------------------------------------------------------------------------
%   AUTHOR:
%       Francesco Fiorina, NavSAS - DET - Politecnico di Torino, 03/04/2026
%
%   VERSION:
%       1.0
%
%   COPYRIGHT:
%       NavSAS - DET - Politecnico di Torino, 2026. All rights reserved.
% -------------------------------------------------------------------------

% --- Extract time vector and arc mask -----------------------------------
time    = obsSolutions(PRNid1).rxTime;
arcMask = logical(obsSolutions(PRNid1).arcMask(:));   % [Nx1] logical

% --- Identify valid (non-NaN) samples within the arc --------------------
% These are the anchor points used to drive the interpolation
validInArc = arcMask & ~isnan(ROTI(:));

edges_arc = diff([false; arcMask(:); false]);
arc_start = find(edges_arc ==  1);
arc_end   = find(edges_arc == -1) - 1;
nArcs     = numel(arc_start);

for k = 1:nArcs
    segRange = arc_start(k):arc_end(k);

    valid_idx = find(validInArc(segRange));
    if isempty(valid_idx)
        arcMask(segRange) = false;
        continue
    end

    % Dual-frequency valid interval within the segment
    Start = valid_idx(1);
    Stop  = valid_idx(end);

    % Invalidate epochs outside the dual-frequency interval
    if Start > 1
        arcMask(segRange(1:Start-1)) = false;
    end
    if Stop < numel(segRange)
        arcMask(segRange(Stop+1:end)) = false;
    end
end

% Need at least 2 valid points to interpolate
if sum(validInArc) < 2
    return
end

% --- Interpolate over the full arc mask using valid anchor points --------
ROTI_tmp = nan(length(time), 1);
ROTI_tmp(arcMask) = interp1( time(validInArc), ...
                              ROTI(validInArc), ...
                              time(arcMask), ...
                              'makima' );

% --- Fill only the gaps: positions inside the arc that are currently NaN -
gapMask       = arcMask & isnan(ROTI(:));
ROTI(gapMask) = ROTI_tmp(gapMask);

end