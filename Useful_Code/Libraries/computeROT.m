function [obsSolutions] = computeROT(obsSolutions, limit, epochDate, settings, varargin)
%COMPUTEROT Computes the Rate Of TEC (ROT) from STEC phase measurements.
%
% DESCRIPTION
% -------------------------------------------------------------------------
% This function estimates the Rate Of Total Electron Content (ROT) from
% Slant Total Electron Content (STEC) phase measurements for all satellites
% available in the input observation structure.
%
% The ROT is defined as the temporal derivative of STEC:
%
%   ROT(t) = d(STEC)/dt     [TECu/s]
%
% The derivative is computed using a central finite difference operator
% applied to contiguous valid arcs of STEC measurements. Data gaps caused
% by cycle slips, tracking loss, or invalid observations are automatically
% detected and handled by segmenting the STEC time series into valid arcs.
%
% Epochs affected by a cycle slip on either frequency are masked to NaN in
% the ROT output. The mask covers the slip epoch itself and the immediately
% following sample, since the central finite difference straddles two
% consecutive epochs.
%
% The computed ROT is stored in the field:
%
%   obsSolutions(i).ROT
%
% for both frequencies of each satellite pair.
%
% -------------------------------------------------------------------------
% INPUT
%
% obsSolutions : struct array
%     Observation solution structure containing at least the fields:
%         .STEC_phase        -> Slant TEC derived from carrier phase
%         .bool_phase_slip   -> Cycle slip flag (logical)
%         .rxTime            -> receiver time vector (seconds)
%         .SV                -> satellite identifier
%         .band              -> signal band label
%
% limit : [startIndex endIndex]
%     Index interval defining the time window used for plotting.
%
% epochDate : datetime array OR cell array
%     Time information used for plotting the time axis.
%
% settings : struct
%     Configuration structure containing at least:
%         .PlotROT  -> Enable ROT plot (1 = yes, 0 = no)
%
% varargin
%     Optional input:
%         figureNumber  -> figure handle for plotting.
%
% -------------------------------------------------------------------------
% OUTPUT
%
% obsSolutions : struct array
%     Updated structure containing the new field:
%
%         .ROT   -> Rate Of TEC [TECu/s], with cycle-slip-affected
%                   epochs set to NaN
%
% -------------------------------------------------------------------------
% METHOD
%
% 1. Identify dual-frequency signal pairs for each satellite
% 2. Extract STEC phase time series
% 3. Detect contiguous valid arcs (NaN gaps split segments)
% 4. Compute ROT using a central finite difference operator
% 5. Build a cycle-slip mask for each frequency (slip epoch and next sample)
%    and combine them with logical OR across both frequencies
% 6. Set ROT to NaN at masked epochs
% 7. Assign the result to both frequencies of the pair
%
% -------------------------------------------------------------------------
% REFERENCES
%
% The Rate Of TEC (ROT) is commonly used for ionospheric monitoring and
% scintillation analysis and is defined as the time derivative of the
% slant TEC derived from dual-frequency GNSS observations.
%
% -------------------------------------------------------------------------
%   AUTHOR:
%       Francesco Fiorina, NavSAS - DET - Politecnico di Torino, 10/03/2026
%
%   VERSION:
%       1.0 - 10/03/2026 - Initial release
%       1.1 - 11/05/2026 - Added cycle-slip masking (slip epoch and next
%                          sample set to NaN); settings struct replaces
%                          plotSwitch scalar input
%
%   COPYRIGHT:
%       NavSAS - DET - Politecnico di Torino, 2026. All rights reserved.
% -------------------------------------------------------------------------

% --- Command Window Log
fprintf('-> Computing Rate Of Tec - ROT...\n')

% --- Build time axis for plotting ---------------------------------------
timeIdx = floor(limit(1)):floor(limit(2));

% Convert epochDate to datetime if necessary
if isdatetime(epochDate)
    time = epochDate(timeIdx);
else
    time = datetime(cell2mat([epochDate(timeIdx,1) epochDate(timeIdx,2)]));
end

% Number of samples required around each point for central derivative
halfspan = 1; % 3-point central difference

% Identify indices corresponding to dual-frequency observations
[prnIdx,~] = FindMultiFreqIndexes(obsSolutions);

% ------------------------------------------------------------------------
% LOOP OVER ALL SATELLITES
% ------------------------------------------------------------------------
for ii = 1:size(prnIdx,2)

    % Indices of the two frequencies (e.g. L1-L5)
    id1 =  prnIdx(1,ii);
    id5 =  prnIdx(2,ii);

    % Skip satellites without valid dual-frequency pair
    if isnan(id1+id5)
        continue
    end

    % Extract STEC phase time series
    STEC_phase = obsSolutions(id1).STEC_phase;

    % Logical mask identifying valid samples
    STEC_valid = ~isnan(STEC_phase);

    % Receiver time vector
    t = obsSolutions(id1).rxTime;

    % Initialize ROT output vectors
    obsSolutions(id1).ROT = nan(size(STEC_valid));
    obsSolutions(id5).ROT = nan(size(STEC_valid));

    % --------------------------------------------------------------------
    % SEGMENTATION OF VALID ARCS
    % --------------------------------------------------------------------
    % Identify start and end indices of contiguous valid arcs

    edges = diff([false; STEC_valid(:); false]);
    seg_start = find(edges == 1);
    seg_end   = find(edges == -1) - 1;

    % --------------------------------------------------------------------
    % PROCESS EACH VALID SEGMENT
    % --------------------------------------------------------------------
    for s = 1:numel(seg_start)

        % Segment limits
        i1 = seg_start(s);
        i2 = seg_end(s);

        % Segment length
        L  = i2 - i1 + 1;

        % Check if segment is long enough for central derivative
        if L < 2*halfspan + 1
            warning(['Segment too short. ', ...
                'Found segment of length %d (indices %d:%d), ', ...
                'but at least %d consecutive valid samples are required for derivative. SEGMENT %d / %d IS SKIPPED for %s'], ...
                L, i1, i2, 2*halfspan + 1, s, numel(seg_start), obsSolutions(id1).SV);
            continue
        end

        % Extract STEC segment
        xs = STEC_phase(i1:i2);

        % ----------------------------------------------------------------
        % COMPUTE ROT
        % ----------------------------------------------------------------
        % Use central finite difference operator
        ROT_arc = finite_diff_central(xs',t(i1:i2));

        % Check for cycle slip
        % --- L1 slip mask (current sample OR previous OR next, no wrap-around) --- previously obsSolutions(id1).bool_phase_slip((i1):(i2)) | circshift(obsSolutions(id1).bool_phase_slip((i1):(i2)), 1);
        bool1 = obsSolutions(id1).bool_phase_slip((i1):(i2)); bool1=bool1(:)';
        shift1_right = [false, bool1(1:end-1)];   % shift di 1 a destra (lookback)
        shift1_left  = [bool1(2:end), false];     % shift di 1 a sinistra (lookahead)
        slip1        = bool1 | shift1_right | shift1_left;
        
        % --- L5 slip mask ---
        bool5 = obsSolutions(id5).bool_phase_slip((i1):(i2)); bool5=bool5(:)';
        shift5_right = [false, bool5(1:end-1)];
        shift5_left  = [bool5(2:end), false];
        slip5        = bool5 | shift5_right | shift5_left;

        slipMask = slip1 | slip5;
        
        % Apply slip mask to ROT
        ROT_arc(slipMask) = NaN;
        
        % Store masked ROT for both frequencies of the pair
        obsSolutions(id1).ROT((i1):(i2)) = ROT_arc;
        obsSolutions(id5).ROT((i1):(i2)) = ROT_arc;

    end
end


%% -----------------------------------------------------------------------
% PLOT ROT TIME SERIES
% -----------------------------------------------------------------------

plotSwitch = settings.PlotROT;

if plotSwitch

    % Create figure
    h = figure();
    figurename = 'RateOfTEC';
    set(h,'Name',figurename);
    hold on;

    % Plot ROT for each satellite
    for i = 1:size(prnIdx,2)

        % Skip invalid entries
        if isnan(prnIdx(1,i))
            continue
        end

        % Satellite ID
        svID = obsSolutions(prnIdx(1,i)).SV;

        % Satellite-specific color
        RGBcolor = PRNColors(svID);

        % Plot ROT
        plot(time,obsSolutions(prnIdx(1,i)).ROT(timeIdx), ...
            'Color',RGBcolor, ...
            'MarkerSize',2, ...
            'DisplayName',append(svID," ", ...
            obsSolutions(prnIdx(1,i)).band,"+", ...
            obsSolutions(prnIdx(2,i)).band), ...
            'LineWidth',2);
    end

    % Axis labels
    xlabel('Time','Interpreter','latex');
    ylabel('ROT [TECu/s]','Interpreter','latex');

    axis tight
    grid on

    % Tick vector for time axis
    ticksVector = getTicks(limit);

    % Plot subtitle showing time interval
    if isdatetime(epochDate)
        subtitle(strcat("Start Date: ",datestr(epochDate(limit(1))), ...
            "  |  End Date: ",datestr(epochDate(limit(2)))),'Interpreter','latex');
    else
        subtitle(strcat("Start Date: ",epochDate{1,1}, ...
            "  |  End Date: ",epochDate{ticksVector(6),1}),'Interpreter','latex');
    end

    % Formatting
    set(gca,'FontSize',13,'TickLabelInterpreter','latex');
    set(findobj(gcf,'type','legend'),'Interpreter','latex','FontSize',13);

    %--- Add listener for datatips
    dcm = datacursormode(gcf);   
    dcm.Enable = 'on';
    dcm.DisplayStyle = 'datatip';
    dcm.UpdateFcn = @customDataTip;

    % Legend
    hLeg = legend('Location','eastoutside','NumColumns',2,'FontSize',8,'Interpreter','latex');
    hLeg.ItemHitFcn = @action1;

    % Window placement
    title("RINEX Phase Obs - Rate of TEC",'Interpreter','latex');

    Pix_SS = get(0,'screensize');
    h.OuterPosition = [Pix_SS(3)/2 0 Pix_SS(3)/2 Pix_SS(4)/2];

    cleanLegend(figurename)
end


% --- Command Window Log
fprintf('Rate Of Tec computation ended.\n\n')

end