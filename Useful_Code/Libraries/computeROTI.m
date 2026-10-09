function [obsSolutions] = computeROTI(obsSolutions,limit, epochDate, settings,varargin)
%COMPUTEROTI Computes the Rate Of TEC Index (ROTI) from ROT measurements.
%
% DESCRIPTION
% -------------------------------------------------------------------------
% This function estimates the Rate Of TEC Index (ROTI) for all satellites
% present in the observation structure. ROTI is defined as the standard
% deviation of the Rate Of TEC (ROT) within a moving time window and is
% widely used to characterize ionospheric irregularities and scintillation.
%
% The algorithm performs the following steps:
%
%   1) Identify valid arcs of ROT observations (NaN gaps are excluded)
%   2) Segment the data to avoid discontinuities due to cycle slips or
%      tracking losses
%   3) Compute ROTI as the moving standard deviation of ROT
%   4) Assign the resulting ROTI values to the corresponding observations
%
% ROTI definition:
%
%       ROTI = std(ROT)
%
% where:
%
%       ROT   = d(STEC)/dt      [TECU/s]
%       ROTI  = standard deviation of ROT within a window
%
% The window length is defined in the settings structure and converted
% into samples using the receiver sampling interval.
%
% INPUT
% -------------------------------------------------------------------------
% obsSolutions : struct
%     Structure containing observation data and previously computed ROT.
%
% limit : [iStart iEnd]
%     Index interval used for plotting.
%
% epochDate : datetime array or cell array
%     Epoch time information used to generate the time axis.
%
% varargin
%     Optional input argument:
%         figureNumber
%
% OUTPUT
% -------------------------------------------------------------------------
% obsSolutions : struct
%     Updated observation structure containing the field:
%
%         obsSolutions(i).ROTI
%
%     which stores the Rate Of TEC Index for each observation epoch.
%
% NOTES
% -------------------------------------------------------------------------
% - ROTI is computed only over contiguous valid ROT arcs.
% - Short segments are skipped to ensure meaningful statistics.
% - The current implementation supports only fixed-length windows.
%
% REFERENCES
% -------------------------------------------------------------------------
% ROTI is commonly used in ionospheric research to detect plasma
% irregularities and scintillation effects in GNSS signals.
% 
% Windowing ROTI calculation based on:
% GOCE High-level Processing Facility
% GOCE Level 1 and Level 2 Reprocessing Campaign
% Technical Note
% -------------------------------------------------------------------------
%   AUTHOR:
%       Francesco Fiorina, NavSAS - DET - Politecnico di Torino, 10/03/2026
%
%   VERSION:
%       1.0
%
%   COPYRIGHT:
%       NavSAS - DET - Politecnico di Torino, 2026. All rights reserved.
% -------------------------------------------------------------------------

% --- Command Window Log
fprintf('-> Computing Rate Of Tec Index - ROTI...\n')

% --- Build time axis for plotting ---------------------------------------
timeIdx = floor(limit(1)):floor(limit(2));

% Convert epochDate to datetime if necessary
if isdatetime(epochDate)
    time = epochDate(timeIdx);
else
    time = datetime(cell2mat([epochDate(timeIdx,1) epochDate(timeIdx,2)]));
end

% Identify indices corresponding to dual-frequency observations
[prnIdx,~] = FindMultiFreqIndexes(obsSolutions);

% ------------------------------------------------------------------------
% ROTI WINDOW DEFINITION
% ------------------------------------------------------------------------
t = obsSolutions(1).rxTime; % Receiver time vector
t_samp= (t(2)-t(1));

switch settings.ROTI_window
    case 'variable'
        error('Method not implemented yet, chose -fixed-.')

    case 'fixed'
        
        window = round(settings.window_length/t_samp); %time window/t_sampling
end

% ------------------------------------------------------------------------
% LOOP OVER ALL SATELLITES
% ------------------------------------------------------------------------
for ii = 1:size(prnIdx,2)

    % Indices of the two frequencies (e.g. L1-L5)
    id1 =  prnIdx(1,ii);
    id5 =  prnIdx(2,ii);

    % Skip satellites without valid dual-frequency pair
    if isnan(id1)
        continue
    end

    % Extract ROT time series
    ROT = obsSolutions(id1).ROT;

    % Logical mask identifying valid samples
    ROT_valid = ~isnan(ROT);

    % Initialize ROT output vectors
    obsSolutions(id1).ROTI = nan(size(ROT_valid));
    obsSolutions(id5).ROTI = nan(size(ROT_valid));

    % --------------------------------------------------------------------
    % SEGMENTATION OF VALID ARCS
    % --------------------------------------------------------------------
    % Identify start and end indices of contiguous valid arcs

    edges = diff([false; ROT_valid(:); false]);
    seg_start = find(edges == 1);
    seg_end   = find(edges == -1) - 1;

%     % --------------------------------------------------------------------
%     % PROCESS EACH VALID SEGMENT - superseeded
%     % --------------------------------------------------------------------
%     for s = 1:numel(seg_start)
% 
%         % Segment limits
%         i1 = seg_start(s);
%         i2 = seg_end(s);
% 
%         % Segment length
%         L  = i2 - i1 + 1;
% 
%         % Check if segment is long enough for sliding window
%         if L < window * (1-settings.tol_window)
%             warning(['Segment too short. ', ...
%                 'Found ROT segment of length %d (indices %d:%d), ', ...
%                 'but at least %d consecutive valid samples are required for good statistics. SEGMENT %d / %d IS SKIPPED for %s'], ...
%                 L, i1, i2, window, s, numel(seg_start), obsSolutions(id1).SV);
%             continue
%         end
% 
%         % Extract STEC segment
%         ROTIwindow = movstd(ROT(i1:i2),window,'omitnan');
% 
%         % Edge effects removal
%         halfW = floor(window/2);
%         ROTIwindow(1:halfW) = NaN;
%         ROTIwindow(end-halfW+1:end) = NaN;
% 
%         % Store ROT for both frequencies of the pair
%         obsSolutions(id1).ROTI((i1):(i2)) = ROTIwindow;
%         obsSolutions(id5).ROTI((i1):(i2)) = ROTIwindow;
% 
%     end
    
    obsSolutions(id1).ROTI = movstd_nancheck(ROT,window,settings.tol_window);
    obsSolutions(id5).ROTI = movstd_nancheck(ROT,window,settings.tol_window);

    % --------------------------------------------------------------------
    % ARC RESTORATION
    % Fills gaps within the arc mask via interpolation.
    % Called AFTER all segments are processed so the full ROTI vector
    % is available as anchor data for interp1.
    % --------------------------------------------------------------------
    if settings.arcRestore
        ROTI_restored = rotiArcRestore(obsSolutions(id1).ROTI, obsSolutions, id1);
        obsSolutions(id1).ROTI_restored = ROTI_restored;
        obsSolutions(id5).ROTI_restored = ROTI_restored;
    end

end

   
%% -----------------------------------------------------------------------
% PLOT ROTI TIME SERIES
% -----------------------------------------------------------------------

plotSwitch = settings.PlotROTI;

if plotSwitch

    % Create figure
    h = figure();
    figurename = 'RateOfTECindex';
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

        % Plot ROTI
        plot(time,obsSolutions(prnIdx(1,i)).ROTI(timeIdx), ...
            'Color',RGBcolor, ...
            'MarkerSize',2, ...
            'DisplayName',append(svID," ", ...
            obsSolutions(prnIdx(1,i)).band,"+", ...
            obsSolutions(prnIdx(2,i)).band), ...
            'LineWidth',2);
    end

    % Axis labels
    xlabel('Time','Interpreter','latex');
    ylabel('ROTI [TECu/s]','Interpreter','latex');

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
    title("RINEX Phase Obs - Rate of TEC Index",'Interpreter','latex');

    Pix_SS = get(0,'screensize');
    h.OuterPosition = [Pix_SS(3)/2 0 Pix_SS(3)/2 Pix_SS(4)/2];

    cleanLegend(figurename)
end



plotSwitch = settings.PlotROTI_restored;

if plotSwitch

    % Create figure
    h = figure();
    figurename = 'RateOfTECindexRestored';
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
        plot(time,obsSolutions(prnIdx(1,i)).ROTI_restored(timeIdx), ...
            'Color',RGBcolor, ...
            'MarkerSize',2, ...
            'DisplayName',append(svID," ", ...
            obsSolutions(prnIdx(1,i)).band,"+", ...
            obsSolutions(prnIdx(2,i)).band), ...
            'LineWidth',2);
    end

    % Axis labels
    xlabel('Time','Interpreter','latex');
    ylabel('ROTI [TECu/s]','Interpreter','latex');

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
    title("RINEX Phase Obs - Rate of TEC Index Restored",'Interpreter','latex');

    Pix_SS = get(0,'screensize');
    h.OuterPosition = [Pix_SS(3)/2 0 Pix_SS(3)/2 Pix_SS(4)/2];

    cleanLegend(figurename)
end



% --- Command Window Log
fprintf('Rate Of Tec Index computation ended. \n\n')

end




%% monnezza
%         while ~end_of_seg
% 
%             i1w = i2w;
%             i2w = i1w + window -1;
% 
%             i2wfuture_check = i2w +window -1;
%             if i2wfuture_check > i2
%                 i2w = i2;
%             end
% 
%             if i2w > i2 || i2w == i2
%                 end_of_seg = true;
%                 i2w = i2;
%             end
%         
% 
%         % Segment length
%         L  = i2w - i1w + 1;
% 
%         % Check if segment is long enough for central derivative
%         if L < 60
%             warning(['Segment too short. ', ...
%                 'Found ROT segment of length %d (indices %d:%d), ', ...
%                 'but at least %d consecutive valid samples are required for good statistics. SEGMENT IS SKIPPED'], ...
%                 L, i1, i2, 60);
%             continue
%         end
% 
%         % Extract STEC segment
%         ROTIwindow = std(ROT(i1w:i2w))*ones(L,1);
% 
% 
%         % Store ROT for both frequencies of the pair
%         obsSolutions(id1).ROTI((i1w):(i2w)) = ROTIwindow;
%         obsSolutions(id5).ROTI((i1w):(i2w)) = ROTIwindow;
%         end