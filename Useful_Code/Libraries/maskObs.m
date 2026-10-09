function obsSolutions = maskObs(obsSolutions, settings)
%FILTEROBSBYQUALITY  Applies CN0 quality mask and extracts valid pseudorange arcs.
%
% DESCRIPTION
% -------------------------------------------------------------------------
% This function filters raw pseudorange and carrier phase measurements from
% dual-frequency GNSS observations based on a minimum C/N0 threshold.
% Observations below the threshold are replaced with NaN. An arc validity
% mask is also computed for each satellite pair, indicating epochs where
% at least one of the two frequencies has a valid pseudorange measurement.
%
% The arc mask is further refined so that, within each contiguous segment,
% only epochs falling inside the dual-frequency overlap interval are kept
% (i.e. from the latest first-valid to the earliest last-valid of L1/L5).
%
% The function operates on all dual-frequency satellite pairs identified
% by FindMultiFreqIndexes and applies constellation-selective processing
% based on the flags defined in the settings structure.
%
% -------------------------------------------------------------------------
% INPUT
%
% obsSolutions : struct array
%     Observation solution structure containing at least the fields:
%         .cn0       -> Carrier-to-noise density ratio [dB-Hz]
%         .rawP      -> Raw pseudorange measurements [m]
%         .rawCarrP  -> Raw carrier phase measurements [cycles]
%         .SV        -> Satellite identifier string
%
% settings : struct
%     Configuration structure containing at least:
%         .CN0_threshold  -> Minimum C/N0 threshold [dB-Hz]
%         .G              -> Enable GPS     (1 = yes, 0 = no)
%         .E              -> Enable Galileo (1 = yes, 0 = no)
%
% -------------------------------------------------------------------------
% OUTPUT
%
% obsSolutions : struct array
%     Updated structure with the following fields modified or added:
%
%         .rawP      -> Pseudorange with sub-threshold epochs set to NaN
%         .rawCarrP  -> Carrier phase with sub-threshold epochs set to NaN
%         .arcMask   -> Logical vector, true where both frequencies have
%                       overlapping valid pseudorange measurements
%
% -------------------------------------------------------------------------
% METHOD
%
% 1. Identify dual-frequency signal pairs via FindMultiFreqIndexes
% 2. For each valid pair, apply the CN0 threshold mask to rawP and rawCarrP
% 3. Compute the per-epoch arc validity mask as the element-wise OR of
%    the two masked pseudorange validity vectors
% 4. Trim each arc segment so that only epochs within the dual-frequency
%    overlap interval (latest start, earliest end) are retained
%
% -------------------------------------------------------------------------
%   AUTHOR:
%       Francesco Fiorina, NavSAS - DET - Politecnico di Torino, 03/04/2026
%
%   VERSION:
%       1.1
%
%   COPYRIGHT:
%       NavSAS - DET - Politecnico di Torino, 2026. All rights reserved.
% -------------------------------------------------------------------------

%--- Keep a list of PRNs present for each of the three bands
[prnIdx, obsSolutions] = FindMultiFreqIndexes(obsSolutions);

%--- Load Off-Boresight Angles
load(settings.obsPath,'tmpResults');

%% Select and restore the observables
for ii = 1:size(prnIdx, 2)

    %--- Check if at least two frequencies are present in the data
    if sum(~isnan(prnIdx(:,ii))) == 2

        valid1 = obsSolutions(prnIdx(1,ii)).cn0 > settings.CN0_threshold;
        valid2 = obsSolutions(prnIdx(2,ii)).cn0 > settings.CN0_threshold;

% Commented and disabled --> kept for future developments
%         if ~strcmp(settings.OP , 'OP77_0')
% 
%             %PRN1   = obsSolutions(prnIdx(1,ii)).PRN;
%             %PRN2 = obsSolutions(prnIdx(2,ii)).PRN;
%             
%             %const1 = extract(obsSolutions(prnIdx(1,ii)).SV,1);
%             %const2 = extract(obsSolutions(prnIdx(2,ii)).SV,1);
%             
%             %tag1 = sprintf('%s%02d',const1,PRN1);
%             %tag2 = sprintf('%s%02d',const2,PRN2);
%             
%             %th1 = settings.mainlobe_threshold.(const1).L1_E1;
%             %th2 = settings.mainlobe_threshold.(const2).L5_E5;
%             
%             %mainTh1 = (90 - tmpResults.timeseries.(tag1).el) < th1;
%             %mainTh2 = (90 - tmpResults.timeseries.(tag2).el) < th2;
%             
%             valid1 = valid1; %& mainTh1';
%             valid2 = valid2; %& mainTh2';
%         end

        P1      = nan(1,length(valid1));
        P2      = nan(1,length(valid2));
        CarrP1  = nan(1,length(valid1));
        CarrP2  = nan(1,length(valid2));

        cString = char(obsSolutions(prnIdx(1,ii)).SV);   % Extracting constellation
        clen    = length(cString);                       % SV+const name length
        if clen > 3                                      % Check whether is long name or short e.g. GPS21 or G21
            clen = 3;
        else
            clen = 1;
        end

        if settings.(cString(1:clen))

            %--- Apply CN0 and lobe mask to pseudorange measurements
            P1(valid1) = obsSolutions(prnIdx(1,ii)).rawP(valid1);
            P2(valid2) = obsSolutions(prnIdx(2,ii)).rawP(valid2);

            obsSolutions(prnIdx(1,ii)).rawP = P1;
            obsSolutions(prnIdx(2,ii)).rawP = P2;

            %--- Apply CN0 and lobe mask to carrier phase measurements
            CarrP1(valid1) = obsSolutions(prnIdx(1,ii)).CarrierPhase(valid1);%.rawCarrP(valid1);
            CarrP2(valid2) = obsSolutions(prnIdx(2,ii)).CarrierPhase(valid2);%.rawCarrP(valid2);

            obsSolutions(prnIdx(1,ii)).CarrierPhase = CarrP1; %.rawCarrP = CarrP1;
            obsSolutions(prnIdx(2,ii)).CarrierPhase = CarrP2; %.rawCarrP = CarrP2;

            %--- Arc mask: valid epoch if either frequency has a valid measurement
            arcMask = ~isnan(P1) | ~isnan(P2);

            %--- Trim arc edges: keep only epochs where both frequencies
            %    are simultaneously active (i.e. from the latest start to
            %    the earliest end of L1 and L5 within each segment)
            L1valid = ~isnan(P1);
            L5valid = ~isnan(P2);

            edges_arc = diff([false; arcMask(:); false]);
            arc_start = find(edges_arc ==  1);
            arc_end   = find(edges_arc == -1) - 1;
            nArcs     = numel(arc_start);

            for k = 1:nArcs
                segRange = arc_start(k):arc_end(k);

                L1_idx = find(L1valid(segRange));
                L5_idx = find(L5valid(segRange));

                % If either frequency has no valid sample, invalidate the whole segment
                if isempty(L1_idx) || isempty(L5_idx)
                    arcMask(segRange) = false;
                    continue
                end

                % Dual-frequency valid interval within the segment
                dualStart = max(L1_idx(1),   L5_idx(1));
                dualEnd   = min(L1_idx(end), L5_idx(end));

                % Invalidate epochs outside the dual-frequency interval
                if dualStart > 1
                    arcMask(segRange(1:dualStart-1)) = false;
                end
                if dualEnd < numel(segRange)
                    arcMask(segRange(dualEnd+1:end)) = false;
                end
            end

            obsSolutions(prnIdx(1,ii)).arcMask = arcMask;
            obsSolutions(prnIdx(2,ii)).arcMask = arcMask;

        end
    end
end

end