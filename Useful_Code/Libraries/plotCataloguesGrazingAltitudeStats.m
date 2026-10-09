function plotCataloguesGrazingAltitudeStats(T, varargin)
%PLOTCATALOGUESGRAZINGALTITUDESTATS  Box-plot statistics of the grazing
%                                    altitude at event start and end
%                                    across all operations (OP).
%
% DESCRIPTION
% -------------------------------------------------------------------------
% Builds a two-panel figure from the aggregated event table produced by
% exportCatalogues.m:
%
%   Left  : aggregate distribution over all PRNs, grouped by event type
%           and split by moment (start/end).
%   Right : same distributions, further split by constellation
%           (GPS / Galileo).
%
% Each box shows the distribution of the grazing altitude observed at
% the first (eventStart) or last (eventEnd) sample of every event of a
% given type across every OP in the aggregated catalogue. Unclassified
% events and rows with NaN grazing altitude are excluded.
%
% INPUT
% -------------------------------------------------------------------------
% T : table OR char/string path to a .mat file containing a table
%     named 'T' (as saved by exportCatalogues.m).
%
%     Required columns:
%       .eventType       : cellstr or string
%       .constellation   : cellstr or string (e.g. 'G','E' or 'GPS','GAL')
%       .grazingAltStart : double [m]
%       .grazingAltEnd   : double [m]
%
% USAGE
% -------------------------------------------------------------------------
%   % from an already-loaded table
%   plotCataloguesGrazingAltitudeStats(T);
%
%   % from a saved .mat file
%   plotCataloguesGrazingAltitudeStats('Results/catalogue_all.mat');
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

%% === Accept either a table or a path to a .mat file ===
if ischar(T) || isstring(T)
    S = load(char(T), 'T');
    if ~isfield(S, 'T')
        error('plotCataloguesGrazingAltitudeStats:BadFile', ...
              'The file does not contain a variable named ''T''.');
    end
    T = S.T;
end

if ~istable(T)
    error('plotCataloguesGrazingAltitudeStats:BadInput', ...
          'Input must be a table or a .mat file path.');
end

required = {'eventType','constellation','grazingAltStart','grazingAltEnd'};
missing  = setdiff(required, T.Properties.VariableNames);
if ~isempty(missing)
    error('plotCataloguesGrazingAltitudeStats:MissingColumns', ...
          'Table is missing required columns: %s', strjoin(missing, ', '));
end

%% === Extract and sanitise data ===
evTypes = string(T.eventType);
consts  = string(T.constellation);

% Tolerate both 'G'/'E' short codes and 'GPS'/'GAL'/'Galileo' long names
consName = strings(size(consts));
consName(consts == "G" | consts == "GPS")                    = "GPS";
consName(consts == "E" | consts == "GAL" | consts == "Galileo") = "Galileo";
consName(consName == "")                                      = consts(consName == "");  % keep others as-is

h_start = T.grazingAltStart / 1e3;    % [km]
h_end   = T.grazingAltEnd   / 1e3;    % [km]

% Drop unclassified and NaN rows
keep = evTypes ~= "unclassified" & ~isnan(h_start) & ~isnan(h_end);
evTypes  = evTypes(keep);
consName = consName(keep);
h_start  = h_start(keep);
h_end    = h_end(keep);

nEv = numel(evTypes);
if nEv == 0
    warning('plotCataloguesGrazingAltitudeStats: no valid events to plot.');
    return
end

%% === Stack start/end samples into one vector for boxchart grouping ===
altitudes = [h_start(:);   h_end(:)];
types     = [evTypes(:);   evTypes(:)];
consStack = [consName(:);  consName(:)];

moments = categorical([repmat("start", nEv, 1); repmat("end", nEv, 1)], ...
                      ["start", "end"], 'Ordinal', true);

% Fix event-type category order so boxes always appear ingress, egress,
% tangent rather than alphabetically shuffled
typeOrder = ["ingress", "egress", "tangent"];
typeOrder = intersect(typeOrder, unique(types), 'stable');
typesCat  = categorical(types, typeOrder, ["ingress", "egress", "tangent"], 'Ordinal', true);

%% === Figure ===
h = figure();
figurename = 'CatalogueGrazingAltitudeStats';
set(h, 'Name', figurename);

Pix_SS = get(0, 'screensize');
h.OuterPosition = [0, 0, Pix_SS(3), Pix_SS(4)/2];

%% === Left panel: aggregate over all PRNs ===
% subplot(1, 2, 1) (uncomment if needed);
boxchart(typesCat, altitudes, 'GroupByColor', moments);

yline(0,   '--k', 'Earth surface',           'Interpreter','latex', ...
      'LabelHorizontalAlignment','left', 'HandleVisibility','off');
yline(80,  ':k',  'Ionosphere base (80 km)', 'Interpreter','latex', ...
      'LabelHorizontalAlignment','left', 'HandleVisibility','off');
yline(350, '-.k', 'F2 peak (~350 km)',       'Interpreter','latex', ...
      'LabelHorizontalAlignment','left', 'HandleVisibility','off');

xlabel('Event type',              'Interpreter','latex');
ylabel('Grazing altitude [km]',   'Interpreter','latex');
% title('All PRNs',                 'Interpreter','latex');
legend('Location','best',         'Interpreter','latex');
grid on
set(gca, 'FontSize', 13, 'TickLabelInterpreter','latex');

%% === Right panel: separated by constellation (uncomment if needed)===
% subplot(1, 2, 2);
% 
% % Combined x-axis label "type / constellation"
% xCombined = strcat(string(types), " / ", consStack);
% 
% % Build the full ordered category list so boxes are grouped first by
% % event type (ingress, egress, tangent) then by constellation within
% % each type
% uniqueConsts = unique(consStack, 'stable');
% fullOrder    = strings(0,1);
% for it = 1:numel(typeOrder)
%     for ic = 1:numel(uniqueConsts)
%         fullOrder(end+1,1) = strcat(typeOrder(it), " / ", uniqueConsts(ic)); %#ok<AGROW>
%     end
% end
% fullOrder = intersect(fullOrder, unique(xCombined), 'stable');
% xCat      = categorical(xCombined, fullOrder, 'Ordinal', true);
% 
% boxchart(xCat, altitudes, 'GroupByColor', moments);
% 
% yline(0,   '--k', 'Earth surface',           'Interpreter','latex', ...
%       'LabelHorizontalAlignment','left', 'HandleVisibility','off');
% yline(80,  ':k',  'Ionosphere base (80 km)', 'Interpreter','latex', ...
%       'LabelHorizontalAlignment','left', 'HandleVisibility','off');
% yline(350, '-.k', 'F2 peak (~350 km)',       'Interpreter','latex', ...
%       'LabelHorizontalAlignment','left', 'HandleVisibility','off');
% 
% xlabel('Event type / constellation',  'Interpreter','latex');
% ylabel('Grazing altitude [km]',       'Interpreter','latex');
% title('By constellation',             'Interpreter','latex');
% legend('Location','best',             'Interpreter','latex');
% grid on
% set(gca, 'FontSize', 13, 'TickLabelInterpreter','latex');
% xtickangle(30);

%% === Overall title ===
sgtitle(sprintf(['Grazing altitude distribution at event start vs end' ...
                 ' --- %d events across all OPs'], nEv), ...
        'Interpreter','latex', 'FontSize', 15);

end