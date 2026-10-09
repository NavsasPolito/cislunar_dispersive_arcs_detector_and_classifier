function plotEventStatistics(catalogue_path)
%PLOTEVENTSTATISTICS  Stacked bar plots of ionospheric event statistics.
%
%   plotEventStatistics()
%   plotEventStatistics(catalogue_path)
%
%   Reads catalogue_all.mat (output of exportCatalogues) and produces:
%     Figure 1 – events per OP stacked by event type (egress/ingress/tangent)
%                with vertical separators between mission phases.
%     Figure 2 – events per OP stacked by constellation (GAL/GPS)
%                with vertical separators between mission phases.
%
%   The phase order is: COMMISSIONING → TRANSIT → LLO → SURFACE
%
%   VERSION HISTORY
%     v1.0  2026-06-08  First release.
% -------------------------------------------------------------------------

%% 0. Defaults
if nargin < 1 || isempty(catalogue_path)
    catalogue_path = fullfile(pwd, 'catalogue_all.mat');
end
assert(isfile(catalogue_path), 'File not found: %s', catalogue_path);

%% 1. Load table
S = load(catalogue_path);
fn = fieldnames(S);
T  = S.(fn{1});   % works whether the variable is 'Tmat' or anything else
fprintf('Loaded %d events from %s\n', height(T), catalogue_path);

%% 2. Sanitise column types (they may be cell arrays of char or string arrays)
T.OP           = string(T.OP);
T.eventType    = lower(string(T.eventType));   % normalise case
T.constellation= upper(string(T.constellation));
T.missionPhase = upper(string(T.missionPhase));

%% 3. Ordered OP list (numeric sort already done by exportCatalogues)
%     Re-derive the numeric key to guarantee order.
allOPs  = unique(T.OP, 'stable');
opNums  = arrayfun(@opToNum, allOPs);
[~, si] = sort(opNums);
allOPs  = allOPs(si);
nOPs    = numel(allOPs);

%% 4. Phase order & phase→OP map
phaseOrder = ["COMMISSIONING","TRANSIT","LLO","SURFACE"];
% Assign each OP its phase (take the mode in case of inconsistencies)
opPhase = strings(nOPs, 1);
for k = 1:nOPs
    mask = T.OP == allOPs(k);
    ph   = T.missionPhase(mask);
    ph   = ph(strlength(ph) > 0);
    if isempty(ph), opPhase(k) = "UNKNOWN";
    else, [~,~,ic] = unique(ph); opPhase(k) = ph(mode(ic)); end
end

% Re-sort OPs by phase order first, then numeric within phase
phIdx = zeros(nOPs, 1);
for k = 1:nOPs
    idx = find(phaseOrder == opPhase(k), 1);
    if isempty(idx), idx = numel(phaseOrder) + 1; end
    phIdx(k) = idx;
end
[~, finalSort] = sortrows([phIdx, opNums(:)]);
allOPs  = allOPs(finalSort);
opPhase = opPhase(finalSort);
opNums  = opNums(finalSort);

%% 5. Event types & constellations
eventTypes = ["egress","ingress","tangent"];
consts     = ["GAL","GPS"];
nET = numel(eventTypes);
nC  = numel(consts);

%% 6. Build count matrices  [nOPs × nCategories]
countET = zeros(nOPs, nET);   % stacked by event type
countC  = zeros(nOPs, nC);    % stacked by constellation

for k = 1:nOPs
    mask = T.OP == allOPs(k);
    for e = 1:nET
        countET(k,e) = sum(mask & T.eventType == eventTypes(e));
    end
    for c = 1:nC
        countC(k,c)  = sum(mask & T.constellation == consts(c));
    end
end

%% 7. Phase separator positions
%     A separator is drawn BEFORE the first OP of each new phase.
[~, phaseStart] = unique(opPhase, 'stable');   % first occurrence of each phase
% Re-order according to phaseOrder
sepPos = [];   % x-positions of separators (between bars)
phaseLabPos = [];
phaseLabels = {};
currentPhase = "";
for k = 1:nOPs
    if opPhase(k) ~= currentPhase
        if k > 1, sepPos(end+1) = k - 0.5; end  %#ok<AGROW>
        phaseLabPos(end+1) = k;                   %#ok<AGROW>
        phaseLabels{end+1} = opPhase(k);          %#ok<AGROW>
        currentPhase = opPhase(k);
    end
end

%% 8. x-tick labels: OP number only (e.g. "OP38_0" → "38\_0")
xLabels = strrep(allOPs, "OP", "");   % strip leading "OP"
xLabels = strrep(xLabels, "_", "\_"); % escape underscore for TeX interpreter

%% 9. Colours
colET = [0.20 0.53 0.74;   % egress  – steel blue
         0.89 0.47 0.20;   % ingress – burnt orange
         0.42 0.68 0.35];  % tangent – muted green

colC  = [0.84 0.19 0.37;   % GAL – red/magenta
         0.20 0.53 0.74];  % GPS – steel blue

%% 10. ---- FIGURE 1: stacked by event type --------------------------------
fig1 = figure('Name','Events per OP – Event type', ...
    'Units','centimeters','Position',[2 10 28 9]);

ax1 = axes(fig1);
bh1 = bar(ax1, 1:nOPs, countET, 'stacked', 'BarWidth', 0.75);
for e = 1:nET
    bh1(e).FaceColor = colET(e,:);
    bh1(e).EdgeColor = 'none';
    bh1(e).DisplayName = char(eventTypes(e));
end

% Phase separators
yLim1 = ax1.YLim;
for s = 1:numel(sepPos)
    xline(ax1, sepPos(s), '--k', 'LineWidth', 0.8, 'HandleVisibility','off');
end

% Phase labels (above the bars)
yTop = max(sum(countET, 2)) * 1.05;
for p = 1:numel(phaseLabPos)
    % find end of this phase block
    if p < numel(phaseLabPos)
        xEnd = phaseLabPos(p+1) - 1;
    else
        xEnd = nOPs;
    end
    xMid = (phaseLabPos(p) + xEnd) / 2;
    text(ax1, xMid, yTop, phaseLabels{p}, ...
        'HorizontalAlignment','center','FontSize',8,'FontWeight','bold', ...
        'Interpreter','none');
end

ax1.XTick = 1:nOPs;
ax1.XTickLabel = xLabels;
ax1.XTickLabelRotation = 45;
ax1.TickLabelInterpreter = 'tex';
ax1.FontSize = 8;
ax1.YLabel.String = 'Number of events';
ax1.YLabel.FontSize = 10;
ax1.XLabel.String   = 'Operation ID';
ax1.XLabel.FontSize = 10;
ax1.Title.String    = 'Classified events per operation – breakdown by event type';
ax1.Title.FontSize  = 11;
grid(ax1,'on'); ax1.Layer = 'top';
legend(ax1, 'Location','northeast','FontSize',9,'Box','off');

ylim(ax1, [0, yTop * 1.12]);
fig1.Color = 'w';
ax1.Color  = [0.97 0.97 0.97];

%% 11. ---- FIGURE 2: stacked by constellation ----------------------------
fig2 = figure('Name','Events per OP – Constellation', ...
    'Units','centimeters','Position',[2 0 28 9]);

ax2 = axes(fig2);
bh2 = bar(ax2, 1:nOPs, countC, 'stacked', 'BarWidth', 0.75);
for c = 1:nC
    bh2(c).FaceColor = colC(c,:);
    bh2(c).EdgeColor = 'none';
    bh2(c).DisplayName = char(consts(c));
end

% Phase separators
for s = 1:numel(sepPos)
    xline(ax2, sepPos(s), '--k', 'LineWidth', 0.8, 'HandleVisibility','off');
end

% Phase labels
yTop2 = max(sum(countC, 2)) * 1.05;
for p = 1:numel(phaseLabPos)
    if p < numel(phaseLabPos), xEnd = phaseLabPos(p+1) - 1;
    else,                       xEnd = nOPs; end
    xMid = (phaseLabPos(p) + xEnd) / 2;
    text(ax2, xMid, yTop2, phaseLabels{p}, ...
        'HorizontalAlignment','center','FontSize',8,'FontWeight','bold', ...
        'Interpreter','none');
end

ax2.XTick = 1:nOPs;
ax2.XTickLabel = xLabels;
ax2.XTickLabelRotation = 45;
ax2.TickLabelInterpreter = 'tex';
ax2.FontSize = 8;
ax2.YLabel.String = 'Number of events';
ax2.YLabel.FontSize = 10;
ax2.XLabel.String   = 'Operation ID';
ax2.XLabel.FontSize = 10;
ax2.Title.String    = 'Classified events per operation – breakdown by constellation';
ax2.Title.FontSize  = 11;
grid(ax2,'on'); ax2.Layer = 'top';
legend(ax2, 'Location','northeast','FontSize',9,'Box','off');

ylim(ax2, [0, yTop2 * 1.12]);
fig2.Color = 'w';
ax2.Color  = [0.97 0.97 0.97];

end

% =========================================================================
%                         LOCAL HELPER
% =========================================================================
function n = opToNum(op)
op = strtrim(char(op));
op = regexprep(op, '^[Oo][Pp]', '');
parts = strsplit(op, '_');
major = str2double(parts{1});
if isnan(major), n = NaN; return; end
if numel(parts) >= 2 && ~isempty(parts{2})
    n = major + str2double(parts{2}) / 10^numel(parts{2});
else
    n = major;
end
end


% =========================================================================
% SNIPPET TO INSERT IN exportCatalogues.m
% =========================================================================
%
% WHERE: paste the function below AFTER the existing printStatistics()
%        function (around line 170 of exportCatalogues.m), and add ONE
%        line in the main body to call it.
%
% CALL SITE: in the main function body, after the line:
%
%     printStatistics(T);         ← already exists (step 8)
%
%   add immediately below it:
%
%     plotEventStatisticsInline(T);
%
% =========================================================================
