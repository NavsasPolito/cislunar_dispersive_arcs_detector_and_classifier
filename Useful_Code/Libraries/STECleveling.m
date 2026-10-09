function obsSolutions = STECleveling(obsSolutions,limit, epochDate, settings,varargin)
%SMOOTHTEC  Carrier-phase smoothing of slant TEC measurements
%
%   obsSolutions = SmoothTEC(obsSolutions)
%
%   This function applies a recursive carrier-phase leveling algorithm
%   to smooth slant TEC measurements. The method combines:
%
%       - Code-derived TEC  (STEC_code)
%       - Phase-derived TEC (STEC_phase)
%
%   The smoothing is based on the recursive averaging of the difference:
%
%       DeltaSTEC = STEC_code - STEC_phase
%
%   The offset is estimated as:
%
%       offset(t) = 1/N * [ (N-1)*offset(t-1) + DeltaSTEC(t) ]
%
%   where N is the number of valid observations since the last cycle slip.
%
%   If a cycle slip is detected (bool_phase_slip == true),
%   the filter is reinitialized.
%
%   NaN measurements do not update the filter but preserve time alignment.
%
%   INPUT:
%       obsSolutions  Struct array containing:
%           .STEC_code         Code-derived slant TEC
%           .STEC_phase        Phase-derived slant TEC
%           .bool_phase_slip   Cycle slip flag (logical)
%
%   OUTPUT:
%       obsSolutions  Updated struct with:
%           .smoothedTEC  Carrier-smoothed slant TEC
%
%   The output vectors preserve the original length and time alignment.
%
% -------------------------------------------------------------------------
% REFERENCES:
% "Wide Area Ionosphere Grid Modellingin the Auroral Region" Susan Skone
% pp. 34-35.
%
% -------------------------------------------------------------------------
%   AUTHOR:
%       Francesco Fiorina, NavSAS - DET - Politecnico di Torino, 04/03/2026
%
%   VERSION:
%       1.0
%
%   COPYRIGHT:
%       NavSAS - DET - Politecnico di Torino, 2026. All rights reserved.
% -------------------------------------------------------------------------

% --- Command Window Log
fprintf('-> Computing Smoothed Slant Total Electron Content - Smoothed STEC...\n')

%--- Channels and Epochs
numChannel = length(obsSolutions);                % Number of tracking channels (satellites/frequencies)
numEpochs  = length(obsSolutions(1).rxTime);      % Number of epochs

%--- Keep a list of PRNs present for each of the three bands
[prnIdx,obsSolutions] = FindMultiFreqIndexes(obsSolutions);

%--- Time axis
timeIdx = floor(limit(1)):floor(limit(2));
if isdatetime(epochDate)
    time = epochDate(timeIdx);
else
    time = datetime(cell2mat([epochDate(timeIdx,1) epochDate(timeIdx,2)]));
end

%--- t_gap tolerance for reinitialization
t_gap = settings.t_gap;

%% -------------------------------------------------------------------
% SMOOTHING CHANNEL BY CHANNEL
% --------------------------------------------------------------------
for ii = 1:numChannel                             % Loop over each satellite channel

    if isempty(obsSolutions(ii).STEC_phase)       % If no STEC estimate is present go to next channel
        continue
    end
    
    % Preallocate output vector for safety and speed
    obsSolutions(ii).smoothedTEC = NaN(1,numEpochs);
    
    % Compute code-phase TEC difference (DeltaSTEC)
    DeltaSTEC = obsSolutions(ii).STEC_code - ...
                obsSolutions(ii).STEC_phase;
    
    N = 0;                                        % Counter of valid samples since last reset

    newObs = 0;                                   % newObs epoch initialization
    lastObs = 0;                                  % lastObs epoch initialization

    firstOccurrence = find(~isnan(DeltaSTEC),1,'first'); % first occurence at which the stec is present
    if isempty(firstOccurrence)
        continue
    end
    
    for jj = 1:numEpochs                          % Loop over time
        
        % Check if phase measurement is valid
        if isnan(obsSolutions(ii).STEC_phase(jj))
            % If phase is invalid, smoothed TEC cannot be computed
            obsSolutions(ii).smoothedTEC(jj) = NaN;
        else
            newObs = jj;
        end
        
        % Detect filter reset conditions:
        % 1) First occurrence
        % 2) First VALID occurrence
        % 3) Cycle slip detected
        % 4) After t_gap hole
        if jj==1 || jj == firstOccurrence || obsSolutions(ii).bool_phase_slip(jj) || newObs-lastObs > t_gap
            
            if ~isnan(DeltaSTEC(jj))
                offset = DeltaSTEC(jj);            % Initialize offset
                N = 1;                             % Reset counter
            else
                offset = NaN;                      % No valid initialization
                N = 0;
            end
            
        else
            % Normal recursive update  
            if ~isnan(DeltaSTEC(jj))               % Update only if valid measurement
                N = N + 1;                         % Increase number of observations
                % Recursive averaging formula
                offset = ((N-1)*offset + DeltaSTEC(jj)) / N; 
            end
            % If DeltaSTEC is NaN:
            %   - Do not update N
            %   - Keep previous offset
        end
        
        % Compute smoothed TEC if offset is valid
        if ~isnan(offset)
            obsSolutions(ii).smoothedTEC(jj) = obsSolutions(ii).STEC_phase(jj) + offset;
        else
            obsSolutions(ii).smoothedTEC(jj) = NaN;
        end
        
        lastObs = newObs;
    end

end


%% -----------------------------------------------------------------------
% PLOT ROT TIME SERIES
% -----------------------------------------------------------------------

plotSwitch = settings.PlotSTECsmooth;

if plotSwitch

    h = figure();
    figurename = 'smoothedSlantTEC';
    set(h,'Name',figurename);
    hold on;


    for i = 1:size(prnIdx,2)
        %--- Find the right color
        if isnan(prnIdx(1,i))
            continue
        end
        svID = obsSolutions(prnIdx(1,i)).SV;
        RGBcolor = PRNColors(svID);
        plot(time,obsSolutions(prnIdx(1,i)).smoothedTEC(timeIdx),'Color',RGBcolor,'MarkerSize', 2,'DisplayName',append(svID," ",obsSolutions(prnIdx(1,i)).band,"+",obsSolutions(prnIdx(2,i)).band),'LineWidth',2);
    end

    xlabel('Time','Interpreter','latex');
    ylabel('Slant TEC [TECu]','Interpreter','latex');
    axis tight;
    grid on;
    ticksVector = getTicks(limit);
%     xticks(time(ismember(timeIdx,ticksVector))); %xticks(time(ticksVector));
%     xticklabels(string(timeofday(time(ismember(timeIdx,ticksVector))))); %xticklabels(string(timeofday(time(ticksVector))));
    if isdatetime(epochDate)
        subtitle(strcat("Start Date: ",datestr(epochDate(limit(1))),"  |  End Date: ",datestr(epochDate(limit(2)))),'Interpreter','latex');
    else
        subtitle(strcat("Start Date: ",epochDate{1,1},"  |  End Date: ",epochDate{ticksVector(6),1}),'Interpreter','latex');
    end

    %--- All elements in LaTeX
    set(gca, 'FontSize', 13, 'TickLabelInterpreter', 'latex');
    set(findobj(gcf,'type','legend'), 'Interpreter', 'latex', 'FontSize', 13);


    %--- Add listener for datatips
    dcm = datacursormode(gcf);   
    dcm.Enable = 'on';
    dcm.DisplayStyle = 'datatip';
    dcm.UpdateFcn = @customDataTip;

    %--- Legend
    hLeg = legend('Location','eastoutside','NumColumns',2,'FontSize',8,'Interpreter','latex');
    hLeg.ItemHitFcn = @action1;

    %--- Place on screen
    title("RINEX Obs - Smoothed Slant TEC",'Interpreter','latex');
    Pix_SS = get(0,'screensize');
    h.OuterPosition = [Pix_SS(3)/2 0 Pix_SS(3)/2 Pix_SS(4)/2];
end

% --- Command Window Log
fprintf('Smoothed STEC computation ended. \n\n')

end