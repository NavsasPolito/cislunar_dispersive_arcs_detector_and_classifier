function [msg, obsSolutions] = PlotTEC_code(directory, file, limit, epochDate, settings, varargin)
% This function generates Code slant TEC plot 
%
% Written by    Simone Zocca
% Extended by:  Francesco Fiorina


% Copyright (C) 2025 NavSAS Group, Politecnico di Torino
% Original code: Simone Zocca (LuNART, https://github.com/NavsasPolito/LuNART)
% Modified by: Francesco Fiorina, 2026
%
% This program is free software: you can redistribute it and/or modify
% it under the terms of the GNU General Public License as published by
% the Free Software Foundation, either version 3 of the License, or
% (at your option) any later version.
%
% This program is distributed in the hope that it will be useful,
% but WITHOUT ANY WARRANTY; without even the implied warranty of
% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
% GNU General Public License for more details.

% --- Command Window Log
fprintf('-> Computing Code Slant Total Electron Content - Code STEC...\n')

%--- Load files
if exist(directory,"dir")
    path = fullfile(directory,file);
    m = load(path);
else
    m.obsSolutions = file;
end

%--- Keep a list of PRNs present for each of the three bands
[prnIdx,m.obsSolutions] = FindMultiFreqIndexes(m.obsSolutions);

%--- Time axis
timeIdx = floor(limit(1)):floor(limit(2));
if isdatetime(epochDate)
    time = epochDate(timeIdx);
else
    time = datetime(cell2mat([epochDate(timeIdx,1) epochDate(timeIdx,2)]));
end

%% Generate plot
plotSwitch = settings.PlotSTECcode;

if plotSwitch
    h = figure();
    figurename = 'codeSlantTEC';
    set(h,'Name',figurename);
end

%--- Compute TEC
for i = 1:size(prnIdx,2)
    %--- Check is at least two frequencies are present in the data
    if sum(~isnan(prnIdx(:,i))) == 2

        cString = char(m.obsSolutions(prnIdx(1,i)).SV); % Extracting constellation
        clen    = length(cString);                      % SV+const name length
        if clen > 3                                     % Check whether is long name or short e.g. GPS21 or G21
            clen = 3;
        else 
            clen = 1;
        end

        if plotSwitch; hold on; end
        if settings.(cString(1:clen))
            %--- Coefficients
            f1 = m.obsSolutions(prnIdx(1,i)).freq;
            f2 = m.obsSolutions(prnIdx(2,i)).freq;
            cf1 = f1^2 / (f1^2 - f2^2);
            cf2 = f2^2 / (f1^2 - f2^2);
            %--- Measurements
            P1= m.obsSolutions(prnIdx(1,i)).rawP;
            P2= m.obsSolutions(prnIdx(2,i)).rawP;
            %--- Compute values
            ionoFree = (P1.*cf1) - (P2.*cf2);
            ionoError = P1 - ionoFree;
            STEC = ionoError.*( f1^2 / 40.3e16 );   % TEC = ( f1^2 * f2^2 / (40.3 * (f1^2  - f2^2) )).*(P1 - P2);
            
            m.obsSolutions(prnIdx(1,i)).STEC_code = STEC;
            m.obsSolutions(prnIdx(2,i)).STEC_code = STEC;
            m.obsSolutions(prnIdx(1,i)).maxSTEC_code = max(STEC(timeIdx));
            m.obsSolutions(prnIdx(2,i)).maxSTEC_code = max(STEC(timeIdx));
          

            %--- Find the right color
            svID = m.obsSolutions(prnIdx(1,i)).SV;
            RGBcolor = PRNColors(svID);
            if plotSwitch
                plot(time,STEC(timeIdx),'Color',RGBcolor,'MarkerSize', 2,'DisplayName',append(svID," ",m.obsSolutions(prnIdx(1,i)).band,"+",m.obsSolutions(prnIdx(2,i)).band),'LineWidth',1.1);
            end
        end
    end
end
obsSolutions = m.obsSolutions;

if plotSwitch
    xlabel('Time','Interpreter','latex');
    ylabel('Slant TEC [TECu]','Interpreter','latex');
    axis tight;
    grid on;
    ticksVector = getTicks(limit);
    %xticks(time(ismember(timeIdx,ticksVector))); %xticks(time(ticksVector));
    %xticklabels(string(timeofday(time(ismember(timeIdx,ticksVector))))); %xticklabels(string(timeofday(time(ticksVector))));
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
    title("RINEX Obs - Code Slant TEC",'Interpreter','latex');
    Pix_SS = get(0,'screensize');
    h.OuterPosition = [Pix_SS(3)/2 0 Pix_SS(3)/2 Pix_SS(4)/2];
end  
    %--- Message
    msg = "Code Slant TEC plot generated";
    % --- Command Window Log
    fprintf('Code STEC computation ended.\n\n')
end