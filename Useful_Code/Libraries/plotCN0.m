function plotCN0(obsSolutions, limit, epochDate, settings, varargin)
%PLOTCN0  Plots C/N0 time series for all dual-frequency satellite pairs.
%
% DESCRIPTION
% -------------------------------------------------------------------------
% This function plots the Carrier-to-Noise density ratio (C/N0) time series
% for each dual-frequency satellite pair found in the observation structure.
% Each satellite is drawn with its unique colour (via PRNColors). A
% horizontal dashed threshold line is overlaid using the value stored in
% settings.CN0_threshold.
%
% Formatting:
%   - LineWidth  : 1.5
%   - FontSize   : 14  (axes tick labels, axis labels, legend)
%   - Interpreter: LaTeX for all text elements
%
% -------------------------------------------------------------------------
% INPUT
%
%   obsSolutions : struct array
%       Observation solution structure. Must contain at least:
%           .cn0        [1×N] C/N0 vector [dB-Hz]
%           .SV         satellite identifier string
%           .band       signal band label
%           .rxTime     receiver time vector (seconds)
%
%   limit : [iStart iEnd]
%       Index interval defining the plotting window.
%
%   epochDate : datetime array OR cell array
%       Time information used to build the time axis.
%
%   settings : struct
%       Must contain:
%           .CN0_threshold   minimum C/N0 threshold [dB-Hz]
%
%   varargin
%       Optional:
%           varargin{1}  -> figure number (default: 800)
%
% -------------------------------------------------------------------------
% OUTPUT
%
%   Figure with C/N0 time series and threshold line.
%
% -------------------------------------------------------------------------
%   AUTHOR:
%       Francesco Fiorina, NavSAS - DET - Politecnico di Torino
%
%   VERSION:
%       1.0
%
%   COPYRIGHT:
%       NavSAS - DET - Politecnico di Torino, 2026. All rights reserved.
% -------------------------------------------------------------------------

if ~settings.PlotCN0
    return
end

% --- Command Window Log --------------------------------------------------
fprintf('-> Plotting C/N0 time series...\n')

% --- Time axis -----------------------------------------------------------
timeIdx = floor(limit(1)):floor(limit(2));

if isdatetime(epochDate)
    time = epochDate(timeIdx);
else
    time = datetime(cell2mat([epochDate(timeIdx,1) epochDate(timeIdx,2)]));
end

% --- Dual-frequency satellite index map ----------------------------------
[prnIdx, obsSolutions] = FindMultiFreqIndexes(obsSolutions);

% --- Create figure -------------------------------------------------------
h = figure();
set(h, 'Name', 'CN0_TimeSeries');
hold on;

% --- Plot C/N0 for each satellite ----------------------------------------
for i = 1:size(prnIdx, 2)

    % Skip incomplete pairs
    if isnan(prnIdx(1,i))
        continue
    end

    id1 = prnIdx(1,i);
    id5 = prnIdx(2,i);

    svID     = obsSolutions(id1).SV;
    RGBcolor = PRNColors(svID);

    dispName1 = append(svID, " ", obsSolutions(id1).band);
    dispName5 = append(svID, " ", obsSolutions(id5).band);

    % --- L1 (or equivalent first frequency) ---
    plot(time, obsSolutions(id1).cn0(timeIdx), ...
        'Color',       RGBcolor, ...
        'LineWidth',   1.5, ...
        'DisplayName', dispName1);

    % --- L5 (or equivalent second frequency): dashed ---
    plot(time, obsSolutions(id5).cn0(timeIdx), ...
        'Color',       RGBcolor, ...
        'LineStyle',   '--', ...
        'LineWidth',   1.5, ...
        'DisplayName', dispName5);
end

% --- Threshold line ------------------------------------------------------
yline(settings.CN0_threshold, ...
    'k:', ...
    'LineWidth', 2, ...
    'DisplayName', sprintf('Threshold: %g dB-Hz', settings.CN0_threshold), ...
    'Interpreter', 'latex', ...
    'LabelHorizontalAlignment', 'left', ...
    'LabelVerticalAlignment',   'bottom');

% --- Axis labels ---------------------------------------------------------
xlabel('Time',            'Interpreter', 'latex', 'FontSize', 14);
ylabel('C/N$_0$ [dB-Hz]', 'Interpreter', 'latex', 'FontSize', 14);

axis tight;
grid on;

% --- Subtitle (time window) ----------------------------------------------
ticksVector = getTicks(limit);

if isdatetime(epochDate)
    subtitle(strcat("Start Date: ", datestr(epochDate(limit(1))), ...
                    "  |  End Date: ",  datestr(epochDate(limit(2)))), ...
        'Interpreter', 'latex');
else
    subtitle(strcat("Start Date: ", epochDate{1,1}, ...
                    "  |  End Date: ",  epochDate{ticksVector(6),1}), ...
        'Interpreter', 'latex');
end

% --- Title ---------------------------------------------------------------
title('RINEX Obs -- C/N$_0$ Time Series', 'Interpreter', 'latex', 'FontSize', 14);

% --- Formatting ----------------------------------------------------------
set(gca, 'FontSize', 14, 'TickLabelInterpreter', 'latex');

% --- Data tips -----------------------------------------------------------
dcm = datacursormode(gcf);
dcm.Enable       = 'on';
dcm.DisplayStyle = 'datatip';
dcm.UpdateFcn    = @customDataTip;

% --- Legend --------------------------------------------------------------
hLeg = legend('Location', 'eastoutside', ...
              'NumColumns', 2, ...
              'FontSize',   10, ...
              'Interpreter','latex');
hLeg.ItemHitFcn = @action1;

% --- Window placement ----------------------------------------------------
Pix_SS = get(0, 'screensize');
h.OuterPosition = [Pix_SS(3)/2, 0, Pix_SS(3)/2, Pix_SS(4)/2];

% --- Command Window Log --------------------------------------------------
fprintf('C/N0 plot generated.\n\n')

end
