function navClosest = selectClosestEphemeris(navUse, tQuery)
%SELECTCLOSESTEPHEMERIS  Select, for each PRN in a RINEX navigation
%                        timetable, the single ephemeris entry whose
%                        reference epoch is closest in time to tQuery.
%
% DESCRIPTION
% -------------------------------------------------------------------------
% MATLAB's gnssconstellation(tQuery, navUse) propagates every single
% ephemeris row in navUse, even if multiple rows are present for the
% same satellite (as is typical in broadcast RINEX, where a new
% ephemeris is issued every 2 h for GPS and every 10 min for Galileo).
% The function then returns one position/velocity entry per row, and
% downstream code must pick the correct one.
%
% This function pre-filters navUse so that, for each PRN, only the
% single ephemeris row whose reference time is closest to tQuery is
% retained. The resulting timetable can be passed to gnssconstellation
% directly and will produce one position/velocity entry per PRN.
%
% The selection is based on the Time column of the timetable (the
% ephemeris reference epoch, also known as t_oe for GPS / t_oc for
% clock parameters). No attempt is made to enforce a maximum age; the
% caller is responsible for checking that the selected ephemeris is
% recent enough to be valid (typical fit interval: 4 h for GPS,
% 3-4 h for Galileo).
%
% INPUT
% -------------------------------------------------------------------------
% navUse : timetable
%     RINEX navigation data as returned by rinexread(...).GPS or
%     rinexread(...).Galileo. Must contain the variables:
%         .Time        : datetime, ephemeris reference epoch
%         .SatelliteID : double,   PRN number
%
% tQuery : datetime (scalar)
%     Query epoch at which the ephemeris will be evaluated.
%
% OUTPUT
% -------------------------------------------------------------------------
% navClosest : timetable
%     Sub-timetable of navUse containing one row per unique PRN,
%     selected as the row whose Time is closest to tQuery.
%
% EXAMPLE
% -------------------------------------------------------------------------
%     rinex      = rinexread(settings.rnxPath);
%     navClosest = selectClosestEphemeris(rinex.GPS, tQuery);
%     [pos, vel, ids] = gnssconstellation(tQuery, navClosest);
%
% -------------------------------------------------------------------------
%   AUTHOR:
%       Francesco Fiorina, NavSAS - DET - Politecnico di Torino
%
%   VERSION HISTORY:
%       1.0 - 22/04/2026 - Initial release.
%
%   COPYRIGHT:
%       NavSAS - DET - Politecnico di Torino, 2026. All rights reserved.
% -------------------------------------------------------------------------

if isempty(navUse)
    navClosest = navUse;
    return
end

prns        = navUse.SatelliteID;
uniquePRNs  = unique(prns);
nPRN        = numel(uniquePRNs);
selectedRow = zeros(nPRN, 1);

for k = 1 : nPRN
    rowsPRN        = find(prns == uniquePRNs(k));
    [~, minIdx]    = min(abs(tQuery - navUse.Time(rowsPRN)));
    selectedRow(k) = rowsPRN(minIdx);
end

navClosest = navUse(selectedRow, :);

end