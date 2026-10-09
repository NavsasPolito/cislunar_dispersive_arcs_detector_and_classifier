function catalogue = eventIDv2_1(obsSolutions, settings)
% EVENTIDV2  Classify GNSS occultation events per satellite.
%
%   catalogue = eventIDv2(obsSolutions, OPid)
%
%   For each dual-frequency satellite, valid ROTI restored arcs are 
%   segmented and classified into one of four event types:
%
%     'ingress'     - Single arc where ROTI rises and signal is lost at peak.
%                     Max occurs at or near the END of the arc.
%                     Mean ROT > +MEAN_THRESHOLD (TEC rising on average).
%
%     'egress'      - Single arc where ROTI decays from a peak at the start.
%                     Max occurs at or near the BEGINNING of the arc.
%                     Mean ROT < -MEAN_THRESHOLD (TEC falling on average).
%
%     'tangent'     - One or more close arcs (gaps < MAX_TANGENT_GAP_S) that
%                     together form a U-shape: ROTI dips to a minimum and
%                     then recovers. Mean ROT ≈ 0 (symmetric pass).
%
%     'unclassified'- Arc does not cleanly match any pattern.
%
%   For a full occultation, ingress and egress appear as two separate arcs
%   separated by a long gap (> OCC_GAP_S). They are stored as two
%   independent catalogue entries.
%
%   Compared to eventIDv1, this version adds a discriminant based on the
%   mean Rate Of TEC (ROT) over the event window. Because ROT = d(STEC)/dt,
%   its sign encodes the direction of ionospheric change:
%       mean ROT > 0  →  TEC is increasing  →  ingress
%       mean ROT < 0  →  TEC is decreasing  →  egress
%       mean ROT ≈ 0  →  symmetric pass     →  tangent
%   This makes the classification more robust to noisy ROTI shapes.
%
%   INPUT
%     obsSolutions  - struct array from FindMultiFreqIndexes
%     OPid          - occultation profile identifier (scalar)
%
%   OUTPUT
%     catalogue     - struct array with fields:
%                       .OP                 occultation profile ID
%                       .SV                 space vehicle ID
%                       .PRN                PRN number
%                       .channel            frequency channel string
%                       .eventTypeGeom      geometric classification from grazing altitude
%                       .eventType          'ingress'|'egress'|'tangent'|'deep tangent'|'unclassified'
%                       .confidence         'high'|'low'|'very low'
%                       .eventStart         UTC datetime of (trimmed) event start
%                       .eventEnd           UTC datetime of (trimmed) event end
%                       .arcStart           UTC datetime of full arc start (before trim)
%                       .arcEnd             UTC datetime of full arc end (before trim)
%                       .grazingAltStart    grazing altitude [m] at event start
%                       .grazingAltEnd      grazing altitude [m] at event end
%                       .grazingRateMedian  median grazing altitude rate [m/s]
%                       .medianROT          median ROT (outlier-filtered) on initial window
%                       .medianROTtrim      median ROT (outlier-filtered) on final trimmed window
%                       .meanROT            mean ROT (outlier-filtered) on initial window
%                       .meanROTtrim        mean ROT (outlier-filtered) on final trimmed window
%                       .peakROTI           maximum ROTI value in the event
%                       .minROTI            minimum ROTI value in the event
%                       .peakSTEC           maximum leveled STEC in the event
%                       .minSTEC            minimum leveled STEC in the event
%                       .segStart           start index/indices of contributing arc(s)
%                       .segEnd             end index/indices of contributing arc(s)
%                       .peakIdx            index of peak ROTI in the full series
%                       .minIdx             index of minimum ROTI in the full series
%                       .trackingStartTime  UTC datetime of tracking window start
%                       .trackingEndTime    UTC datetime of tracking window end
%                       .trackingStartAlt   grazing altitude [m] at tracking start
%                       .trackingEndAlt     grazing altitude [m] at tracking end
%                       .trackingStartIdx   sample index of tracking window start
%                       .trackingEndIdx     sample index of tracking window end
%                       .override           0 not overridden, 1 overridden by look ahead logic, 2 overridden by look back logic
%
% -------------------------------------------------------------------------
%   AUTHOR:
%       Francesco Fiorina, NavSAS - DET - Politecnico di Torino
%
%   VERSION HISTORY:
%       1.0 - 17/03/2026 - Initial release (eventIDv1): shape-based
%             classification using peak position at arc edges and
%             U-shape detection for tangent events.
%       2.0 - 13/04/2026 - Added meanROT discriminant. Fixed uninitialized
%             trimStart/trimStop on empty arc. Added post-trim consistency
%             check with automatic recomputation of peak, min and meanROT.
%             All configuration parameters moved to initSettings.
%       2.1 - 14/04/2026 - Added look-back consistency check to prevent
%             two consecutive ingress or egress for the same satellite.
%             Fixed copy-paste error in gapToPrevGroup. Extended helper
%             with soft labels ('probable ingress/egress/tangent') for
%             ambiguous arcs; probable labels are stripped to base type
%             and stored with a separate .confidence field ('high'|'low').
%             Fixed PEAK_EDGE_TOL hardcoded inside helper overriding the
%             value passed from initSettings.
%       2.2 - 13/05/2026 - Switched ROT discriminant from mean to median
%             with quartile-based outlier removal (MEDIAN_THRESHOLD).
%             Catalogue now stores both initial-window and trimmed-window
%             values for medianROT and meanROT. Added 'deep tangent' event
%             type for passes with interior ROTI peak above both arc tails
%             (maxIsInterior && maxIsHigher) and low |medianROT|. Added
%             geometric validation block (eventTypeGeom) using grazing
%             altitude profile, rate, and DEEP_TANGENT_HEIGHT threshold.
%             Added restoredROTIclass switch: when false, segmentation
%             uses ROTI_restored mask but classification uses original
%             (non-interpolated) ROTI values. Look-ahead now computes
%             peak-at-edge on the full merged event window instead of
%             the first segment only. Look-back now checks SV identity
%             to prevent cross-satellite contamination. Added arcStart
%             and arcEnd fields storing the full (pre-trim) arc time
%             boundaries. Added tracking window fields (trackingStart/
%             EndTime, trackingStart/EndAlt, trackingStart/EndIdx) and
%             grazingRateMedian to catalogue. Sidelobe events
%             (grazingAltStart > 7 Mm) are discarded with a warning.
%             Added grazing altitude distribution boxplot and event/
%             tracking marker overlay on GrazingAltitudeArcs figure.
%       2.3 - 14/04/2026 - Added physical consistency checks (look-ahead
%             + look-back) implementing the ingress→[occultation]→egress
%             sequence constraint. Look-ahead now uses a safe trigger:
%             intervenes only if the current event has NOT already been
%             classified with high confidence as ingress or egress,
%             preventing degradation of correct main-loop classifications.
%             Look-back enforces asymmetric override on consecutive same-
%             type events (ingress→ingress overrides current; egress→
%             egress overrides previous) to preserve the physical order,
%             with same-SV guard to avoid cross-satellite comparisons.
%             Fixed copy-paste error in gapToPrev. Extended helper with
%             soft labels for ambiguous arcs; probable labels stripped
%             to base type with separate .confidence field ('high'|'low').
%             Helper now uses peakAtEnd/peakAtStart (consistent with main
%             loop) for high-confidence ingress/egress detection; old
%             endVal/startVal criterion kept commented for reference.
%             Added 'omitnan' to max() call for robustness on arcs with
%             NaN gaps. Removed unreachable 'deep tangent' branch (subset
%             of standard tangent), kept commented for future use. Fixed
%             PEAK_EDGE_TOL hardcoded inside helper overriding the value
%             passed from initSettings. Added override history.
%
%   COPYRIGHT:
%       NavSAS - DET - Politecnico di Torino, 2026. All rights reserved.
% -------------------------------------------------------------------------

% --- Command Window Log
fprintf('-> Event classification...\n')

% -------------------------------------------------------------------------
% CONFIGURATION  (all parameters loaded from settings)
% -------------------------------------------------------------------------
MAX_TANGENT_GAP_S   = settings.MAX_TANGENT_GAP_S;
OCC_GAP_S           = settings.OCC_GAP_S;
PEAK_EDGE_TOL       = settings.PEAK_EDGE_TOL;
MEDIAN_THRESHOLD    = settings.MEDIAN_THRESHOLD;
MIN_DURATION        = settings.MIN_DURATION;
DEEP_TANGENT_HEIGHT = settings.DEEP_TANGENT_HEIGHT; 

% -------------------------------------------------------------------------
% SETUP
% -------------------------------------------------------------------------
[prnIdx, ~] = FindMultiFreqIndexes(obsSolutions);

catalogue = struct('OP',                 {}, ...
                   'SV',                 {}, ...
                   'PRN',                {}, ...
                   'channel',            {}, ...
                   'eventTypeGeom',      {}, ...
                   'eventType',          {}, ...
                   'confidence',         {}, ...
                   'eventStart',         {}, ...
                   'eventEnd',           {}, ...
                   'arcStart',           {}, ...
                   'arcEnd',             {}, ...
                   'grazingAltStart',    {}, ... 
                   'grazingAltEnd',      {}, ...  
                   'grazingRateMedian',  {}, ...
                   'medianROT',          {}, ...
                   'medianROTtrim',      {}, ...
                   'meanROT',            {}, ...  
                   'meanROTtrim',        {}, ...
                   'peakROTI',           {}, ...
                   'minROTI',            {}, ...
                   'peakSTEC',           {}, ...
                   'minSTEC',            {}, ...
                   'segStart',           {}, ...
                   'segEnd',             {}, ...
                   'peakIdx',            {}, ...
                   'minIdx',             {}, ...
                   'trackingStartTime',  {}, ...   
                   'trackingEndTime',    {}, ...   
                   'trackingStartAlt',   {}, ...   
                   'trackingEndAlt',     {}, ...
                   'trackingStartIdx',   {}, ...
                   'trackingEndIdx',     {}, ...
                   'override',           {});      
idx = 1;

UTC_Time = datetime(gps_2_Utc([obsSolutions(1).weekNumber',obsSolutions(1).GPStime']));
% -------------------------------------------------------------------------
% LOOP OVER ALL DUAL-FREQUENCY SATELLITE PAIRS
% -------------------------------------------------------------------------
for ii = 1:size(prnIdx, 2)

    id1 = prnIdx(1, ii);
    id5 = prnIdx(2, ii);

    if isnan(id1 + id5)
        continue
    end

    if settings.restoredROTIclass
        ROTI           = obsSolutions(id1).ROTI_restored;
        ROTI_valid     = ~isnan(ROTI);
    else
        ROTI_restored  = obsSolutions(id1).ROTI_restored;
        ROTI_valid     = ~isnan(ROTI_restored);
    
        ROTI           = obsSolutions(id1).ROTI;
    end

    ROT                = obsSolutions(id1).ROT;

    % ---------------------------------------------------------------------
    % SEGMENT THE ROTI SERIES INTO CONTIGUOUS VALID ARCS
    % ---------------------------------------------------------------------
    edges     = diff([false; ROTI_valid(:); false]);
    seg_start = find(edges ==  1);
    seg_end   = find(edges == -1) - 1;
    nSeg      = numel(seg_start);

    if nSeg == 0
        continue
    end

    % ---------------------------------------------------------------------
    % COMPUTE PER-SEGMENT FEATURES
    % (gap to NEXT segment, peak position relative to arc edges)
    % ---------------------------------------------------------------------
    segFeats = struct();
    for s = 1:nSeg
        seg = ROTI(seg_start(s):seg_end(s));
        L   = numel(seg);

        [peakVal, peakRel] = max(seg);
        [minVal,  minRel ] = min(seg);

        tol = max(1, round(PEAK_EDGE_TOL * L));   % tolerance in samples

        segFeats(s).peakVal    = peakVal;
        segFeats(s).peakIdx    = seg_start(s) + peakRel - 1;                    % absolute index
        segFeats(s).minVal     = minVal;
        segFeats(s).minIdx     = seg_start(s) + minRel  - 1;                    % absolute index
        segFeats(s).peakAtEnd  = peakRel >= (L - tol);                          % max near end   → ingress
        segFeats(s).peakAtStart= peakRel <= (1 + tol);                          % max near start → egress
        segFeats(s).len        = L;

        % Gap to the NEXT segment (in samples = seconds at 1 Hz)
        if s < nSeg
            segFeats(s).gapToNext = seg_start(s+1) - seg_end(s) - 1;
        else
            segFeats(s).gapToNext = Inf;
        end
    end

    % ---------------------------------------------------------------------
    % GROUP SEGMENTS INTO CANDIDATE EVENTS
    %
    % Strategy:
    %   1. Greedily merge consecutive segments whose mutual gap is
    %      < MAX_TANGENT_GAP_S → candidate tangent groups.
    %   2. Any remaining pair of groups separated by a long gap
    %      (> OCC_GAP_S) is treated as an ingress+egress occultation.
    %   3. Isolated single segments are classified individually.
    % ---------------------------------------------------------------------

    % --- Step 1: group closely-spaced segments ---
    groups    = {};   % each cell holds a vector of segment indices
    g         = 1;
    groups{g} = 1;
    for s = 1:nSeg-1
        if segFeats(s).gapToNext <= MAX_TANGENT_GAP_S
            groups{g}(end+1) = s + 1;   % merge into current group
        else
            g = g + 1;
            groups{g} = s + 1;          % start a new group
        end
    end
    nGroups = numel(groups);

    % --- Step 2: classify each group ---

    for g = 1:nGroups
        sIdx = groups{g};               % segment indices in this group
        nS   = numel(sIdx);

        % Global start/end of the merged event (absolute sample indices)
        evStart = seg_start(sIdx(1));
        evEnd   = seg_end(  sIdx(end));

        if (evEnd - evStart) < MIN_DURATION
            continue
        end

        % Compute features on the initial window
        evROTI  = ROTI(evStart:evEnd);
        evROT   = ROT(evStart:evEnd);

        [peakROTI, peakRel] = max(evROTI, [], 'omitnan');
        [minROTI,  minRel ] = min(evROTI, [], 'omitnan');
        peakIdxAbs = evStart + peakRel - 1;
        minIdxAbs  = evStart + minRel  - 1;
        medianROT = median(evROT(~isoutlier(evROT,"quartiles")),'omitnan');
        meanROT   = mean(evROT(~isoutlier(evROT,"quartiles")),'omitnan');

        segStartVec = seg_start(sIdx);
        segEndVec   = seg_end(  sIdx);

        % -----------------------------------------------------------------
        % CLASSIFICATION RULES
        % -----------------------------------------------------------------
        if nS > 1
            % Multiple close segments → Tangent
            %   (the U-shape diagnostic: min is interior to the merged arc,
            %    not at either extremity)
            [evType,trimStart,trimStop] = classifyTangentOrUnclassified(evROTI, minRel, medianROT, minIdxAbs, MEDIAN_THRESHOLD, PEAK_EDGE_TOL);

            if (strcmp(evType, 'ingress') || strcmp(evType, 'probable ingress')) && ~isnan(trimStart)
                evStart = trimStart;
            elseif (strcmp(evType, 'egress') || strcmp(evType, 'probable egress')) && ~isnan(trimStop)
                evEnd = trimStop;
            elseif (strcmp(evType,'unclassified'))
                warning('[eventIDv2] %s — Unclassified event from epoch %d to epoch %d. Group %d. ',...
                    obsSolutions(id1).SV, evStart, evEnd, g);
            end

        else
            % Single segment: use peak position at arc edges
            s1 = sIdx(1);

            if segFeats(s1).peakAtEnd && ~segFeats(s1).peakAtStart && medianROT >= MEDIAN_THRESHOLD
                evType  = 'ingress';

                evStart = minIdxAbs;

            elseif segFeats(s1).peakAtStart && ~segFeats(s1).peakAtEnd && medianROT <= -MEDIAN_THRESHOLD
                evType = 'egress'; 

                evEnd   = minIdxAbs;

            elseif (segFeats(s1).peakAtStart && ~segFeats(s1).peakAtEnd && abs(medianROT) < MEDIAN_THRESHOLD) ||...
                   (segFeats(s1).peakAtEnd && ~segFeats(s1).peakAtStart && abs(medianROT) < MEDIAN_THRESHOLD)
                evType = 'tangent';

            else
                % Single arc but peak is interior → could be a short
                % tangent pass that had no data gaps
                % Peak is interior, or shape/ROT contradict each other
                [evType, trimStart, trimStop] = classifyTangentOrUnclassified( ...
                    evROTI, minRel, medianROT, minIdxAbs, MEDIAN_THRESHOLD, PEAK_EDGE_TOL);

                if (strcmp(evType, 'ingress') || strcmp(evType, 'probable ingress')) && ~isnan(trimStart)
                    evStart = trimStart;
                elseif (strcmp(evType, 'egress') || strcmp(evType, 'probable egress')) && ~isnan(trimStop)
                    evEnd = trimStop;
                elseif (strcmp(evType,'unclassified'))
                    warning('[eventIDv2] %s — Unclassified event from epoch %d to epoch %d. Group %d. ',...
                        obsSolutions(id1).SV, evStart, evEnd, g);
                end
            end
        end

        % -----------------------------------------------------------------
        % PHYSICAL CONSISTENCY CHECKS
        % -----------------------------------------------------------------
        % Before consolidating the classification, the just-classified
        % event is checked against the physical sequence constraint of an
        % occultation. By construction, a full occultation produces two
        % consecutive arcs separated by a long gap (≈ OCC_GAP_S, typically
        % ~30 min, while the satellite is behind the Earth):
        %
        %     ingress ──[long gap, occultation]──> egress
        %
        % Two manifestations of this constraint are enforced:
        %
        %   LOOK-AHEAD  (forward consistency)
        %     If a long gap follows the current group, the current group
        %     must be the leading edge of an occultation (ingress) or the
        %     trailing edge of the previous one (egress). Triggered only
        %     when the current group has NOT already been classified with
        %     high confidence as ingress or egress, to avoid degrading
        %     correct main-loop classifications. The shape (peak position)
        %     decides which of the two.
        %
        %   LOOK-BACK   (backward consistency)
        %     If the previous group (same satellite) is close in time
        %     (gap < OCC_GAP_S), the two groups belong to the same
        %     occultation event and their types must respect the physical
        %     sequence: an ingress cannot follow an ingress, an egress
        %     cannot follow an egress. Violations are corrected
        %     asymmetrically to preserve the ingress→egress order:
        %       ingress→ingress : current is overridden to egress
        %       egress→egress   : previous is overridden to ingress
        %                         (overriding current would yield the
        %                         impossible egress→ingress sequence)
        %
        % Tangent events are physically neutral (the satellite grazes the
        % ionosphere without being occulted), so they do not constrain
        % the sequence and are never targeted by these checks.
        % -----------------------------------------------------------------

        override = 0;   % reset per-group flag

        % --- Look-ahead (next group) ---
        if g < nGroups
            gapToNextGroup = seg_start(groups{g+1}(1)) - seg_end(sIdx(end)) - 1;

        % -----------------------------------------------------------------
        % SAFE TRIGGER RATIONALE
        % -----------------------------------------------------------------
        % A long gap after the current group is ambiguous: it can mean
        %   (a) we are inside an ongoing occultation (current group is
        %       ingress or egress, the gap is the satellite hiding behind
        %       the Earth), or
        %   (b) the current group is unrelated to an occultation
        %       (e.g. a real tangent pass that happens to be followed
        %       by an unrelated quiet interval before the next event).
        %
        % Distinguishing (a) from (b) requires the shape and medianROT,
        % which the main loop has already evaluated with strict criteria.
        % If the main loop returned a high-confidence label (ingress,
        % egress, or tangent), we trust it and skip the look-ahead.
        % Look-ahead intervenes only when the main loop was uncertain
        % (soft 'probable *' labels or 'unclassified'), where the
        % occultation evidence from the long gap is genuinely useful.
        % -----------------------------------------------------------------
        % Safe trigger: only intervene if the current event has NOT been
        % settled with high confidence by the main loop. A settled
        % classification means the main loop's strict criteria were met
        % (shape AND medianROT agreeing for ingress/egress; U-shape AND
        % near-zero medianROT for tangent). Soft labels ('probable *')
        % and 'unclassified' remain eligible for look-ahead correction.
        currentIsSettled = strcmp(evType, 'ingress') || ...
                           strcmp(evType, 'egress')  || ...
                           strcmp(evType, 'tangent');

            if gapToNextGroup > OCC_GAP_S && ~currentIsSettled

                L_ev   = numel(evROTI);
                tol_ev = max(1, round(PEAK_EDGE_TOL * L_ev));
                evPeakAtEnd   = peakRel >= (L_ev - tol_ev);
                evPeakAtStart = peakRel <= (1 + tol_ev);

                if evPeakAtEnd && medianROT >= MEDIAN_THRESHOLD
                    warning('[eventIDv2] %s — Group %d (was %s) overridden to probable ingress (long gap follows).', ...
                            obsSolutions(id1).SV, g, evType);
                    evType   = 'probable ingress';
                    evStart  = minIdxAbs;
                    override = 1;

                elseif evPeakAtStart && medianROT <= -MEDIAN_THRESHOLD
                    warning('[eventIDv2] %s — Group %d (was %s) overridden to probable egress (long gap follows).', ...
                            obsSolutions(id1).SV, g, evType);
                    evType   = 'probable egress';
                    evEnd    = minIdxAbs;
                    override = 1;
                end
            end
        end

        % --- OLD LOOK-AHEAD (kept for reference) ---
        % if g < nGroups
        %     gapToNextGroup = seg_start(groups{g+1}(1)) - seg_end(sIdx(end)) - 1;
        %     if gapToNextGroup > OCC_GAP_S && abs(medianROT) > MEDIAN_THRESHOLD
        %         L_ev   = numel(evROTI);
        %         tol_ev = max(1, round(PEAK_EDGE_TOL * L_ev));
        %         evPeakAtEnd   = peakRel >= (L_ev - tol_ev);
        %         evPeakAtStart = peakRel <= (1 + tol_ev);
        %         if evPeakAtEnd && medianROT >= MEDIAN_THRESHOLD && ~contains(evType,'ingress')
        %             if contains(evType,'egress')
        %                 warning('[eventIDv2] %s — Event from epoch %d to epoch %d (Group %d) overridden to ingress from egress.', ...
        %                         obsSolutions(id1).SV, evStart, evEnd, g);
        %             else
        %                 warning('[eventIDv2] %s — Event from epoch %d to epoch %d (Group %d) overridden to ingress from tangent.', ...
        %                         obsSolutions(id1).SV, evStart, evEnd, g);
        %             end
        %             evType      = 'probable ingress';
        %             evStart     = minIdxAbs;
        %             override    = 1;
        %         elseif evPeakAtStart && medianROT <= -MEDIAN_THRESHOLD && ~contains(evType,'egress')
        %             if contains(evType,'ingress')
        %                 warning('[eventIDv2] %s — Event from epoch %d to epoch %d (Group %d) overridden to egress from ingress.', ...
        %                         obsSolutions(id1).SV, evStart, evEnd, g);
        %             else
        %                 warning('[eventIDv2] %s — Event from epoch %d to epoch %d (Group %d) overridden to egress from tangent.', ...
        %                         obsSolutions(id1).SV, evStart, evEnd, g);
        %             end
        %             evType      = 'probable egress';
        %             evEnd       = minIdxAbs;
        %             override    = 1;
        %         end
        %     end
        % end

        % --- Look-back (previous group) ---
        if g > 1 && idx > 1
            prevEntry  = catalogue(idx-1);
            prevSegEnd = prevEntry.segEnd(end);
            gapToPrev  = seg_start(sIdx(1)) - prevSegEnd - 1;

            sameSV = strcmp(prevEntry.SV, obsSolutions(id1).SV);

            if gapToPrev < OCC_GAP_S && abs(medianROT) > MEDIAN_THRESHOLD && sameSV

                % Case 1: ingress → ingress
                %   Override CURRENT to egress (sequence becomes ingress → egress)
                if strcmp(evType, 'ingress') && strcmp(prevEntry.eventType, 'ingress')
                    evType   = 'probable egress';
                    evEnd    = minIdxAbs;
                    override = 2;
                    warning(['[eventIDv2] %s — two consecutive ingress detected. ', ...
                             'Current group %d overridden to probable egress.'], ...
                             obsSolutions(id1).SV, g);

                % Case 2: egress → egress
                %   Override PREVIOUS to ingress (sequence becomes ingress → egress)
                %   Cannot override current to ingress because egress → ingress
                %   would be physically inconsistent.
                elseif strcmp(evType, 'egress') && strcmp(prevEntry.eventType, 'egress')
                    catalogue(idx-1).eventType  = 'probable ingress';
                    catalogue(idx-1).confidence = 'low';
                    catalogue(idx-1).override   = 2;
                    warning(['[eventIDv2] %s — two consecutive egress detected. ', ...
                             'Previous group %d overridden to probable ingress.'], ...
                             obsSolutions(id1).SV, g-1);
                end
            end
        end

        % --- OLD LOOK-BACK (kept for reference) ---
        % if g > 1 && idx > 1
        %     prevEntry    = catalogue(idx-1);
        %     prevSegEnd   = prevEntry.segEnd(end);
        %     gapToPrev    = seg_start(sIdx(1)) - prevSegEnd - 1;
        %     if gapToPrev < OCC_GAP_S && abs(medianROT) > MEDIAN_THRESHOLD && strcmp(prevEntry.SV, obsSolutions(id1).SV)
        %         if strcmp(evType, 'ingress') && strcmp(prevEntry.eventType, 'ingress')
        %             evType = 'probable egress';
        %             evEnd  = minIdxAbs;
        %             warning(['[eventIDv2] %s — two consecutive ingress detected. ', ...
        %                      'Group %d overridden to egress from ingress.'], ...
        %                     obsSolutions(id1).SV, g);
        %         elseif strcmp(evType, 'egress') && strcmp(prevEntry.eventType, 'egress')
        %             evType  = 'probable ingress';
        %             evStart = minIdxAbs;
        %             warning(['[eventIDv2] %s — two consecutive egress detected. ', ...
        %                      'Group %d overridden to ingress from egress.'], ...
        %                     obsSolutions(id1).SV, g);
        %         end
        %     end
        % end
        
        % -----------------------------------------------------------------
        % POST-TRIM CONSISTENCY CHECK
        % Verify that peak and min indices fall within the trimmed window.
        % If not, recompute them on the trimmed interval and warn.
        % -----------------------------------------------------------------
        peakOutside = (peakIdxAbs < evStart) || (peakIdxAbs > evEnd);
        minOutside  = (minIdxAbs  < evStart) || (minIdxAbs  > evEnd);
        trimApplied = (evStart ~= seg_start(sIdx(1))) || (evEnd ~= seg_end(sIdx(end)));

        % Recompute on trimmed window if anything changed
        medianROT_trim = medianROT;
        meanROT_trim   = meanROT;

        if peakOutside || minOutside || trimApplied

            % Recompute all derived quantities on the trimmed window
            evROTI_trimmed = ROTI(evStart:evEnd);
            evROT_trimmed  = ROT(evStart:evEnd);  

            [peakROTI_trim, peakRel_trim] = max(evROTI_trimmed, [], 'omitnan');
            [minROTI_trim,  minRel_trim ] = min(evROTI_trimmed, [], 'omitnan');
            peakIdxAbs_trim = evStart + peakRel_trim - 1;
            minIdxAbs_trim  = evStart + minRel_trim  - 1;
            medianROT_trim  = median(evROT_trimmed(~isoutlier(evROT_trimmed,"quartiles")), 'omitnan');
            meanROT_trim    = mean(evROT_trimmed(~isoutlier(evROT_trimmed,"quartiles")), 'omitnan'); 

            if peakOutside
                warning(['[eventIDv2] %s — peak index fell outside trimmed window. ', ...
                    'Original: idx=%d (ROTI=%.4f), Trimmed: idx=%d (ROTI=%.4f). ', ...
                    'Offset: %d samples.'], ...
                    obsSolutions(id1).SV, ...
                    peakIdxAbs,      peakROTI, ...
                    peakIdxAbs_trim, peakROTI_trim, ...
                    abs(peakIdxAbs - peakIdxAbs_trim));
                peakROTI   = peakROTI_trim;
                peakIdxAbs = peakIdxAbs_trim;
            end

            if minOutside
                warning(['[eventIDv2] %s — min index fell outside trimmed window. ', ...
                    'Original: idx=%d (ROTI=%.4f), Trimmed: idx=%d (ROTI=%.4f). ', ...
                    'Offset: %d samples.'], ...
                    obsSolutions(id1).SV, ...
                    minIdxAbs,      minROTI, ...
                    minIdxAbs_trim, minROTI_trim, ...
                    abs(minIdxAbs - minIdxAbs_trim));
                minROTI   = minROTI_trim;
                minIdxAbs = minIdxAbs_trim;
            end
    
        end

        % Confidence level:
        if contains(evType, 'probable')
            confidence = 'low';
            evType     = strrep(evType, 'probable ', '');  % strip → 'ingress'/'egress'/'tangent'
        elseif contains(evType, 'unclassified')
            confidence = 'very low';
        else
            confidence = 'high';
        end

        % -----------------------------------------------------------------
        % EVENT TIME INFO
        % -----------------------------------------------------------------

        timeStart = UTC_Time(evStart);
        timeEnd   = UTC_Time(evEnd);

        timeArcStart = UTC_Time(seg_start(sIdx(1)));
        timeArcEnd   = UTC_Time(seg_end(sIdx(end)));

        % -----------------------------------------------------------------
        % GEOMETRIC VALIDATION
        % -----------------------------------------------------------------
        grazingAlt              = obsSolutions(id1).grazingAltitude(evStart:evEnd);
        grazingAltRate          = finite_diff_central(grazingAlt);
        grazingAltRateMedian    = median(grazingAltRate(~isoutlier(grazingAltRate,"quartiles")), 'omitnan');
        grazingAltStart         = grazingAlt(1);
        grazingAltEnd           = grazingAlt(end);
        grazingAltMin           = min(grazingAlt);

        if grazingAltStart > grazingAltEnd && grazingAltMin == grazingAltEnd && grazingAltRateMedian <= -1000
            evTypeGeom = 'ingress';
        elseif grazingAltStart < grazingAltEnd && grazingAltMin == grazingAltStart && grazingAltRateMedian >= 1000
            evTypeGeom = 'egress';
        elseif grazingAltMin <= grazingAltStart && grazingAltMin <= grazingAltEnd && grazingAltMin > DEEP_TANGENT_HEIGHT 
            evTypeGeom = 'tangent';
        elseif grazingAltMin <= grazingAltStart && grazingAltMin <= grazingAltEnd && grazingAltMin < DEEP_TANGENT_HEIGHT 
            evTypeGeom = 'tangent';%'deep tangent';
        else
            evTypeGeom = 'unclassified';
        end
        
        % Commented and disabled --> kept for future developments
        %if (grazingAltStart > settings.grazing_threshold && strcmp(settings.OP , 'OP77_0')) %grazingAltStart > settings.grazing_threshold || 
        %    warning('possible sidelobe event --> discarded and not classified')
        %    continue
        %end

        % -----------------------------------------------------------------
        % LOCATE THE TRACKING WINDOW CONTAINING THIS EVENT
        % Each event lies inside exactly one tracking window of the SV.
        % If no such window exists (should not happen, since events are
        % only built on valid arcs), record NaN.
        % -----------------------------------------------------------------
        
        if isfield(obsSolutions,'SVtrackingPoints')
            trackingPoints = obsSolutions(id1).SVtrackingPoints;
            trkMatch       = find(trackingPoints(:,1) <= evStart & ...
                                  trackingPoints(:,2) >= evEnd, 1);
        else
            trkMatch       = [];
        end
        
        if isempty(trkMatch) || ~isfield(obsSolutions,'SVtrackingPoints')
            trackingStartIdx  = NaN;
            trackingEndIdx    = NaN;
            trackingStartAlt  = NaN;
            trackingEndAlt    = NaN;

            trackingStartTime = NaN;
            trackingEndTime   = NaN;
        else
            trackingStartIdx  = trackingPoints(trkMatch, 1);
            trackingEndIdx    = trackingPoints(trkMatch, 2);
            trackingStartAlt  = trackingPoints(trkMatch, 3);
            trackingEndAlt    = trackingPoints(trkMatch, 4);

            trackingStartTime = UTC_Time(trackingStartIdx);
            trackingEndTime   = UTC_Time(trackingEndIdx);
        end

        % -----------------------------------------------------------------
        % STORE IN CATALOGUE
        % -----------------------------------------------------------------
        catalogue(idx).OP                 = settings.OP;
        catalogue(idx).SV                 = obsSolutions(id1).SV;
        catalogue(idx).PRN                = obsSolutions(id1).PRN;
        catalogue(idx).channel            = 'L1+L5';
        catalogue(idx).eventTypeGeom      = evTypeGeom;
        catalogue(idx).eventType          = evType;
        catalogue(idx).confidence         = confidence;
        catalogue(idx).eventStart         = timeStart;
        catalogue(idx).eventEnd           = timeEnd;
        catalogue(idx).arcStart           = timeArcStart;
        catalogue(idx).arcEnd             = timeArcEnd;
        catalogue(idx).grazingAltStart    = grazingAltStart;
        catalogue(idx).grazingAltEnd      = grazingAltEnd;
        catalogue(idx).grazingRateMedian  = grazingAltRateMedian;
        catalogue(idx).medianROT          = medianROT;
        catalogue(idx).medianROTtrim      = medianROT_trim; % final trimmed event
        catalogue(idx).meanROT            = meanROT;
        catalogue(idx).meanROTtrim        = meanROT_trim;   % final trimmed event
        catalogue(idx).peakROTI           = peakROTI;     
        catalogue(idx).minROTI            = minROTI;
        catalogue(idx).peakSTEC           = max(obsSolutions(id1).smoothedTEC(evStart:evEnd));
        catalogue(idx).minSTEC            = min(obsSolutions(id1).smoothedTEC(evStart:evEnd));
        catalogue(idx).segStart           = segStartVec;
        catalogue(idx).segEnd             = segEndVec;
        catalogue(idx).peakIdx            = peakIdxAbs;
        catalogue(idx).minIdx             = minIdxAbs;
        catalogue(idx).trackingStartTime  = trackingStartTime;
        catalogue(idx).trackingEndTime    = trackingEndTime;
        catalogue(idx).trackingStartAlt   = trackingStartAlt;
        catalogue(idx).trackingEndAlt     = trackingEndAlt;
        catalogue(idx).trackingStartIdx   = trackingStartIdx;
        catalogue(idx).trackingEndIdx     = trackingEndIdx;
        catalogue(idx).override           = override;
        idx = idx + 1; 

    end  % group loop
end      % satellite loop

%% -----------------------------------------------------------------------
% PLOT Grazing Altitude distribution at event start / end
% -------------------------------------------------------------------------
% Builds two boxplot panels from the event catalogue:
%   - Left  : aggregate over all PRNs, grouped by event type and
%             start/end.
%   - Right : separated by constellation (GPS / Galileo), same grouping.
%
% Each box shows the distribution of the grazing altitude observed at
% the first (evStart) or last (evEnd) sample of every event of a given
% type. Unclassified events are excluded.
% -------------------------------------------------------------------------
plotSwitch = settings.PlotGrazingAltitudeDistribution;
if plotSwitch

    % === Collect data from the catalogue ===
    evTypes   = {catalogue.eventType};
    svStrings = {catalogue.SV};                   % 'G##' / 'E##'
    const     = cellfun(@(s) s(1), svStrings);    % first char: 'G' or 'E'

    h_start   = [catalogue.grazingAltStart] / 1e3;   % [km]
    h_end     = [catalogue.grazingAltEnd]   / 1e3;   % [km]

    % Drop unclassified entries and NaN rows (e.g. failed sanity check)
    keep = ~strcmp(evTypes, 'unclassified') & ...
           ~isnan(h_start) & ~isnan(h_end);

    evTypes = evTypes(keep);
    const   = const(keep);
    h_start = h_start(keep);
    h_end   = h_end(keep);

    nEv = numel(evTypes);
    if nEv == 0
        warning('plotGrazingAltitudeDistribution: no valid events to plot.');
        return
    end

    % Stack start/end samples into one vector for boxplot grouping
    altitudes = [h_start(:); h_end(:)];
    types     = [evTypes(:); evTypes(:)];
    consts    = [const(:);   const(:)];
    moments = categorical([repmat("start", nEv, 1); repmat("end", nEv, 1)], ...
                      ["start", "end"], 'Ordinal', true);

    % === Figure ===
    h = figure();
    figurename = 'GrazingAltitudeDistribution';
    set(h, 'Name', figurename);

    Pix_SS = get(0, 'screensize');
    h.OuterPosition = [0, 0, Pix_SS(3), Pix_SS(4)/2];

    % === Left panel: aggregate over all PRNs ===
    %subplot(1, 2, 1);
    boxchart(categorical(types), altitudes, ...
             'GroupByColor', moments);

    % Reference lines
    yline(0,   '--k', 'Earth surface',           'Interpreter','latex', ...
          'LabelHorizontalAlignment','left', 'HandleVisibility','off');
    yline(80,  ':k',  'Ionosphere base (80 km)', 'Interpreter','latex', ...
          'LabelHorizontalAlignment','left', 'HandleVisibility','off');
    yline(350, '-.k', 'F2 peak (~350 km)',       'Interpreter','latex', ...
          'LabelHorizontalAlignment','left', 'HandleVisibility','off');

    xlabel('Event type',                  'Interpreter','latex');
    ylabel('Grazing altitude [km]',       'Interpreter','latex');
    legend('Location','best', 'Interpreter','latex');
    grid on
    set(gca, 'FontSize', 13, 'TickLabelInterpreter', 'latex'); 

    % === Overall title ===
    sgtitle('Grazing altitude distribution at event start vs end', ...
            'Interpreter','latex', 'FontSize', 15);
end

%% -----------------------------------------------------------------------
% OVERLAY classified event markers on the pre-existing arcs plot
% -------------------------------------------------------------------------
% If the "GrazingAltitudeArcs" figure was created upstream (by the
% grazing-altitude plotting block), decorate it here with GREEN filled
% squares at each classified event's start and RED filled squares at
% its end. Unclassified events and NaN grazing altitudes are skipped.
% Markers inherit the SV's DisplayName so that the legend's ItemHitFcn
% (action1) toggles the arc together with its markers.
% -------------------------------------------------------------------------
if isfield(settings,'PlotGrazingAltitudeArcs') && ...
        settings.PlotGrazingAltitudeArcs && ~isempty(catalogue)

    figurename = 'GrazingAltitudeArcs';
    hArcs      = findobj('Type','figure','Name',figurename);

    if isempty(hArcs)
        warning('eventIDv2_1:NoArcsFigure', ...
                'Figure "%s" not found: skipping event-marker overlay.', ...
                figurename);
    else
        figure(hArcs);
        hold on;

        % --- SVs actually drawn on the axes -----------------------
        ax          = gca;
        existingSVs = string({});
        if ~isempty(ax.Children)
            existingSVs = string(get(ax.Children, 'DisplayName'));
        end
        existingSVs = existingSVs(existingSVs ~= "");

        % --- Catalogue data ---------------------------------------
        catSV     = string({catalogue.SV});
        catType   = string({catalogue.eventType});
        catTStart = [catalogue.eventStart];
        catTEnd   = [catalogue.eventEnd];
        catHStart = [catalogue.grazingAltStart] / 1e3;
        catHEnd   = [catalogue.grazingAltEnd]   / 1e3;

        keepCat = (catType ~= "unclassified") & ...
                  ~isnan(catHStart) & ~isnan(catHEnd);

        catSV     = catSV(keepCat);
        catTStart = catTStart(keepCat);
        catTEnd   = catTEnd(keepCat);
        catHStart = catHStart(keepCat);
        catHEnd   = catHEnd(keepCat);

        % --- Overlay markers, one SV at a time --------------------
        for s = 1:numel(existingSVs)
            svName = char(existingSVs(s));     % <-- CHAR, matches arc DisplayName
            sel    = (catSV == existingSVs(s));
            if ~any(sel), continue; end

            plot(catTStart(sel), catHStart(sel), ...
                 's', 'MarkerSize', 8, ...
                 'MarkerFaceColor', [0 0.6 0], ...
                 'MarkerEdgeColor', 'k', ...
                 'LineStyle',       'none', ...
                 'DisplayName',     svName, ...         
                 'HandleVisibility','off');

            plot(catTEnd(sel), catHEnd(sel), ...
                 's', 'MarkerSize', 8, ...
                 'MarkerFaceColor', [0.85 0 0], ...
                 'MarkerEdgeColor', 'k', ...
                 'LineStyle',       'none', ...
                 'DisplayName',     svName, ...         
                 'HandleVisibility','off');
        end

        % --- Overlay tracking-window markers (stars), one SV at a time ---
        catTrkStartIdx = [catalogue.trackingStartIdx];
        catTrkEndIdx   = [catalogue.trackingEndIdx];
        catTrkStartAlt = [catalogue.trackingStartAlt] / 1e3;   % km
        catTrkEndAlt   = [catalogue.trackingEndAlt]   / 1e3;   % km

        catTrkStartTime = NaT(size(catTrkStartIdx));
        catTrkEndTime   = NaT(size(catTrkEndIdx));
        validTrk = ~isnan(catTrkStartIdx) & ~isnan(catTrkEndIdx);
        catTrkStartTime(validTrk) = UTC_Time(catTrkStartIdx(validTrk));
        catTrkEndTime(validTrk)   = UTC_Time(catTrkEndIdx(validTrk));

        catSVfull = string({catalogue.SV});
        keepTrk   = validTrk & ~isnan(catTrkStartAlt) & ~isnan(catTrkEndAlt);

        for s = 1:numel(existingSVs)
            svName = char(existingSVs(s));
            sel    = (catSVfull == existingSVs(s)) & keepTrk;
            if ~any(sel), continue; end

            % Deduplicate by tracking-window start to avoid stacking
            % multiple stars where several events share one tracking window
            [~, uniq] = unique(catTrkStartIdx(sel));
            selIdx    = find(sel);
            selIdx    = selIdx(uniq);

            plot(catTrkStartTime(selIdx), catTrkStartAlt(selIdx), ...
                 'p', 'MarkerSize', 12, ...
                 'MarkerFaceColor', [0 0.6 0], ...
                 'MarkerEdgeColor', 'k', ...
                 'LineStyle',       'none', ...
                 'DisplayName',     svName, ...
                 'HandleVisibility','off');

            plot(catTrkEndTime(selIdx), catTrkEndAlt(selIdx), ...
                 'p', 'MarkerSize', 12, ...
                 'MarkerFaceColor', [0.85 0 0], ...
                 'MarkerEdgeColor', 'k', ...
                 'LineStyle',       'none', ...
                 'DisplayName',     svName, ...
                 'HandleVisibility','off');
        end

        % --- Dummy handles for the legend entries -----------------
        plot(NaT, NaN, 's', 'MarkerFaceColor',[0 0.6 0], ...
             'MarkerEdgeColor','k', 'MarkerSize',8, 'LineStyle','none', ...
             'DisplayName','event start');
        plot(NaT, NaN, 's', 'MarkerFaceColor',[0.85 0 0], ...
             'MarkerEdgeColor','k', 'MarkerSize',8, 'LineStyle','none', ...
             'DisplayName','event end');
        plot(NaT, NaN, 'p', 'MarkerFaceColor',[0 0.6 0], ...
             'MarkerEdgeColor','k', 'MarkerSize',12, 'LineStyle','none', ...
             'DisplayName','tracking start');
        plot(NaT, NaN, 'p', 'MarkerFaceColor',[0.85 0 0], ...
             'MarkerEdgeColor','k', 'MarkerSize',12, 'LineStyle','none', ...
             'DisplayName','tracking end');

        % Refresh legend so the new entries appear
        legend('Location','eastoutside', 'NumColumns',2, ...
               'FontSize',8, 'Interpreter','latex');
        if exist('cleanLegend','file')
            cleanLegend(figurename);
        end
    end
end

%% --- Command Window Log
fprintf('Classification completed.\n')

end  % main function


% =========================================================================
% LOCAL HELPER: classify a single merged arc as tangent or unclassified
%
%   A tangent pass shows a U-shape: the minimum is interior (not at the
%   first or last valid sample), and the values at both ends are
%   meaningfully higher than the minimum.
% =========================================================================
function [evType, trimStart, trimStop] = classifyTangentOrUnclassified( ...
        arcROTI, minRel, medianROT, minIdxAbs, MEDIAN_THRESHOLD, PEAK_EDGE_TOL)

validIdx = find(~isnan(arcROTI));
if isempty(validIdx)
    evType    = 'unclassified';
    trimStart = NaN;
    trimStop  = NaN;
    return
end

L   = numel(validIdx);
tol = max(1, round(PEAK_EDGE_TOL * L));

% Max search robust to NaN gaps (index returned in arcROTI coordinates)
[maxVal, maxRel] = max(arcROTI, [], 'omitnan');

maxIsInterior = (maxRel > tol) && (maxRel < (numel(arcROTI) - tol));
minIsInterior = (minRel > tol) && (minRel < (numel(arcROTI) - tol));

% Edge-position diagnostics (consistent with main loop segFeats)
peakAtEnd   = maxRel >= (numel(arcROTI) - tol);
peakAtStart = maxRel <= (1 + tol);

% Values at the two ends of the valid data
startVal = arcROTI(validIdx(1));
endVal   = arcROTI(validIdx(end));
minVal   = arcROTI(minRel);

% U-shape and inverted-U checks
bothEndsHigher = (startVal > minVal) && (endVal > minVal);
maxIsHigher    = (startVal < maxVal) && (endVal < maxVal);

if minIsInterior && bothEndsHigher && abs(medianROT) < MEDIAN_THRESHOLD
    evType    = 'tangent';
    trimStart = NaN;
    trimStop  = NaN;

elseif minIsInterior && bothEndsHigher && maxIsInterior && maxIsHigher && abs(medianROT) < MEDIAN_THRESHOLD   % DEAD CODE: subset of previous branch, kept for future 'deep tangent' differentiation
    evType    = 'tangent';   % 'deep tangent';
    trimStart = NaN;
    trimStop  = NaN;

elseif peakAtEnd && medianROT >= MEDIAN_THRESHOLD   % was: endVal > startVal && medianROT >= MEDIAN_THRESHOLD
    evType    = 'ingress';
    trimStart = minIdxAbs;
    trimStop  = NaN;

elseif peakAtStart && medianROT <= -MEDIAN_THRESHOLD   % was: startVal > endVal && medianROT <= -MEDIAN_THRESHOLD
    evType    = 'egress';
    trimStart = NaN;
    trimStop  = minIdxAbs;

elseif abs(medianROT) < MEDIAN_THRESHOLD
    evType    = 'probable tangent';
    trimStart = NaN;
    trimStop  = NaN;

elseif medianROT >= MEDIAN_THRESHOLD
    evType    = 'probable ingress';
    trimStart = minIdxAbs;
    trimStop  = NaN;

elseif medianROT <= -MEDIAN_THRESHOLD
    evType    = 'probable egress';
    trimStart = NaN;
    trimStop  = minIdxAbs;

else
    evType    = 'unclassified';
    trimStart = NaN;
    trimStop  = NaN;
end

end