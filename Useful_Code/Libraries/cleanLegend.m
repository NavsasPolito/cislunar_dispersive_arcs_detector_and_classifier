function cleanLegend(figureName)
%CLEANLEGEND  Removes duplicate and hidden legend entries from a named figure.
%
% DESCRIPTION
% -------------------------------------------------------------------------
% This function cleans the legend of a figure identified by its Name
% property. It removes:
%
%   1) Entries whose DisplayName is empty or starts with an underscore '_'
%      (MATLAB convention for "hidden from legend").
%   2) Duplicate DisplayName entries — only the first occurrence is kept.
%
% -------------------------------------------------------------------------
% USAGE
%
%   cleanLegend(figureName)
%
%   where figureName is the string assigned to the figure's Name property,
%   e.g. 'RateOfTECindex' or 'CN0_TimeSeries'.
%
% -------------------------------------------------------------------------
% INPUT
%
%   figureName : char or string
%       The Name property of the target figure (set via
%       set(h,'Name', figureName) or figure('Name', figureName)).
%
% -------------------------------------------------------------------------
% EXAMPLE
%
%   cleanLegend('RateOfTECindex')
%   cleanLegend('CN0_TimeSeries')
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

% --- Input validation ----------------------------------------------------
if nargin < 1 || isempty(figureName)
    error('cleanLegend: figureName is required. Pass the Name property of the target figure.');
end
figureName = char(figureName);

% --- Find figure by Name -------------------------------------------------
allFigs = findall(0, 'Type', 'figure');
fig     = [];
for k = 1:numel(allFigs)
    if strcmp(allFigs(k).Name, figureName)
        fig = allFigs(k);
        break
    end
end

if isempty(fig)
    error('cleanLegend: no open figure with Name = ''%s'' was found.', figureName);
end

% --- Find legend(s) in the figure ----------------------------------------
legs = findobj(fig, 'Type', 'legend');
if isempty(legs)
    warning('cleanLegend: no legend found in figure ''%s''.', figureName);
    return
end

% --- Process each legend -------------------------------------------------
for iLeg = 1:numel(legs)

    ax = legs(iLeg).Axes;
    if isempty(ax)
        continue
    end

    kids = get(ax, 'Children');

    seenNames = {};

    for k = 1:numel(kids)
        obj = kids(k);

        if ~isprop(obj, 'DisplayName')
            continue
        end

        dname = strtrim(obj.DisplayName);

        % Rule 1: empty name or leading underscore -> hide
        if isempty(dname) || startsWith(dname, '_')
            obj.Annotation.LegendInformation.IconDisplayStyle = 'off';
            continue
        end

        % Rule 2: duplicate -> hide after first occurrence
        if any(strcmp(seenNames, dname))
            obj.Annotation.LegendInformation.IconDisplayStyle = 'off';
        else
            seenNames{end+1} = dname; %#ok<AGROW>
            obj.Annotation.LegendInformation.IconDisplayStyle = 'on';
        end
    end

    % Refresh legend preserving original position and column count
    prevLocation   = legs(iLeg).Location;
    prevNumColumns = legs(iLeg).NumColumns;
    prevFontSize   = legs(iLeg).FontSize;

    legend(ax, 'off');
    hLeg = legend(ax, ...
        'Location',    prevLocation, ...
        'NumColumns',  prevNumColumns, ...
        'FontSize',    prevFontSize, ...
        'Interpreter', 'latex');
    hLeg.ItemHitFcn = @action1;

end

fprintf('cleanLegend: legend cleaned for figure ''%s''.\n', figureName);

end
