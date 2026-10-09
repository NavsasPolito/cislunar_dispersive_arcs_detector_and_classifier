function trackingPoints = retrieveTrackingPoints(y, settings)
%RETRIEVETRACKINGPOINTS  Identify tracking-window points (start/end of
%                        contiguous visibility) from a NaN-masked time
%                        series.
%
% DESCRIPTION
% -------------------------------------------------------------------------
% Given a time series y in which NaN denotes "satellite not tracked at
% this epoch" (e.g. grazing altitude masked by L1+L5 availability), this
% function detects contiguous non-NaN segments, merges segments that are
% separated by a gap shorter than settings.MAX_TANGENT_GAP_S samples
% (treated as a transient track loss within a single physical pass),
% and returns one tracking-window record per merged group.
%
% Each tracking-window record contains the index and the y-value of the
% first valid sample of the group (start) and of the last valid sample
% (end). For example, when y is the geometric grazing altitude, these
% records represent the start and end of a single contiguous tracking
% pass of the satellite, regardless of which frequencies are tracked.
%
% INPUT
% -------------------------------------------------------------------------
% y        : double vector
%     Time series with NaN values marking epochs of no tracking.
%
% settings : struct with field
%     .MAX_TANGENT_GAP_S : maximum gap (in samples) between consecutive
%                          segments for them to be merged into a single
%                          tracking window.
%
% OUTPUT
% -------------------------------------------------------------------------
% trackingPoints : N-by-4 double matrix, one row per merged group:
%     Column 1 : startIdx  (sample index of group start)
%     Column 2 : endIdx    (sample index of group end)
%     Column 3 : startVal  (y value at startIdx)
%     Column 4 : endVal    (y value at endIdx)
%
%     Empty (0-by-4) when no valid tracking is present in y.
%
% -------------------------------------------------------------------------
%   AUTHOR:
%       Francesco Fiorina, NavSAS - DET - Politecnico di Torino
%
%   VERSION HISTORY:
%       1.0 - 23/04/2026 - Initial release. Handles 0, 1 and N segments
%             uniformly; pre-allocates output; validates required
%             settings field.
%
%   COPYRIGHT:
%       NavSAS - DET - Politecnico di Torino, 2026. All rights reserved.
% -------------------------------------------------------------------------

%% === Input validation ===
if ~isfield(settings, 'MAX_TANGENT_GAP_S')
    error('retrieveTrackingPoints:MissingSetting', ...
          'settings must contain the field MAX_TANGENT_GAP_S.');
end

%% === Segment the series into contiguous non-NaN arcs ===
y         = y(:);
y_valid   = ~isnan(y);
edges     = diff([false; y_valid; false]);
seg_start = find(edges ==  1);
seg_end   = find(edges == -1) - 1;
nSeg      = numel(seg_start);

%% === Early exit: no valid samples ===
if nSeg == 0
    trackingPoints = zeros(0, 4);
    return
end

%% === Group closely-spaced segments ===
% Two segments are merged when the gap between them is shorter than
% settings.MAX_TANGENT_GAP_S samples.
groups    = {};
g         = 1;
groups{g} = 1;
for s = 1 : nSeg-1
    gapToNext = seg_start(s+1) - seg_end(s) - 1;
    if gapToNext <= settings.MAX_TANGENT_GAP_S
        groups{g}(end+1) = s + 1;        % merge into current group
    else
        g = g + 1;
        groups{g} = s + 1;               % start a new group
    end
end
nGroups = numel(groups);

%% === Build output (one row per merged group) ===
trackingPoints = zeros(nGroups, 4);
for g = 1 : nGroups
    grp                  = groups{g};
    startIdx             = seg_start(grp(1));
    endIdx               = seg_end(grp(end));
    trackingPoints(g, 1) = startIdx;
    trackingPoints(g, 2) = endIdx;
    trackingPoints(g, 3) = y(startIdx);
    trackingPoints(g, 4) = y(endIdx);
end

end