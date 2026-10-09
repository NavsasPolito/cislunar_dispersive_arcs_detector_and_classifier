function T = exportCatalogues(resultsRoot, outDir, missionCSV)
%EXPORTCATALOGUES  Aggregate per-OP catalogue.mat files and export tables.
%
%   T = exportCatalogues()
%   T = exportCatalogues(resultsRoot, outDir, missionXML)
%
%   WORKFLOW
%     1. Parse Mission_Scheduler.xml → OP-to-phase mapping
%     2. Scan OP*_* folders, discover ALL fields dynamically
%     3. GUI #1: select which OPs to include (by phase or individually)
%     4. GUI #2: select which columns per output format
%     5. Sort: OP (numeric) → constellation → PRN
%     6. Export .mat / .csv / .tex / .docx with merged cells
%     7. Print classification statistics to command window
% -------------------------------------------------------------------------
%   AUTHOR:
%       Francesco Fiorina, NavSAS - DET - Politecnico di Torino
%
%   VERSION HISTORY
%     v4.0  2026-05-12  Full rewrite: mission phases, dynamic event types,
%                        dual-GUI (OP + column), statistics summary.
%
%   COPYRIGHT:
%       NavSAS - DET - Politecnico di Torino, 2026. All rights reserved.
% -------------------------------------------------------------------------

%% 0. Defaults
if nargin<1||isempty(resultsRoot)
    resultsRoot = fullfile('/Users/francescofiorina/Library/CloudStorage', ...
        'OneDrive-PolitecnicodiTorino/FrancescoFiorina/PhD', ...
        'Works_OnGoing/Ionosphere_characterization_LuGRE/Results');
end
if nargin<2||isempty(outDir),    outDir = resultsRoot; end
if nargin<3||isempty(missionCSV)
    dataRoot = fullfile('/Users/francescofiorina/Library/CloudStorage', ...
        'OneDrive-PolitecnicodiTorino/FrancescoFiorina/PhD', ...
        'Works_OnGoing/Ionosphere_characterization_LuGRE/Datasets');
    missionCSV = fullfile(dataRoot, 'OPTABLE.csv');
    % If not found next to Results, try inside Useful_Code
    if ~isfile(missionCSV)
        warning('OP table not found at %s. Phase grouping disabled.', missionCSV);
    end
end

assert(isfolder(resultsRoot), 'Results root not found: %s', resultsRoot);
if ~isfolder(outDir), mkdir(outDir); end

%% 1. Parse OP table (CSV) → OP-to-phase map
opPhaseMap = containers.Map();  % 'OP23_0' → 'LLO'
phaseOPs   = struct();          % phaseOPs.LLO = {'OP23_0','OP24_0',...}
phases     = {};                % ordered phase names

% Single-letter phase codes used in the OP table → full mission phase names
phaseCode = containers.Map( {'C','T','L','S'}, ...
                            {'COMMISSIONING','TRANSIT','LLO','SURFACE'} );

if isfile(missionCSV)
    fprintf('Parsing OP table: %s\n', missionCSV);

    opts = detectImportOptions(missionCSV, 'Delimiter', ',');
    opts = setvartype(opts, 'char');        % keep '1_0' as text, not a number
    opts.ExtraColumnsRule = 'ignore';
    opts.EmptyLineRule    = 'skip';
    Tops = readtable(missionCSV, opts);

    % Locate the OP-ID and phase columns by header, ignoring case/spaces/units
    hdr      = regexprep(lower(opts.VariableNames), '[^a-z0-9]', '');
    colOP    = find(strcmp(hdr, 'opid'), 1);
    colPhase = find(strcmp(hdr, 'phase'), 1);
    assert(~isempty(colOP) && ~isempty(colPhase), ...
        'OP table must contain "OP ID" and "Phase" columns: %s', missionCSV);

    for r = 1:height(Tops)
        opID = strtrim(Tops{r, colOP}{1});              % e.g. '23_0'
        code = upper(strtrim(Tops{r, colPhase}{1}));    % e.g. 'L'
        if isempty(opID) || isempty(code), continue; end

        if isKey(phaseCode, code)
            pName = phaseCode(code);
        else
            pName = matlab.lang.makeValidName(code);    % already-expanded or new code
        end

        if ~isfield(phaseOPs, pName)        % first row of this phase
            phases{end+1}    = pName; %#ok<AGROW>
            phaseOPs.(pName) = {};
        end

        opFolder = ['OP' opID];
        opPhaseMap(opFolder)    = pName;
        phaseOPs.(pName){end+1} = opFolder;
    end
    fprintf('  Phases: %s\n', strjoin(phases, ', '));
    fprintf('  %d OPs mapped.\n', opPhaseMap.Count);
else
    warning('OP table not found at %s. Phase grouping disabled.', missionCSV);
end

%% 2. Discover OP folders and fields
d = dir(fullfile(resultsRoot, 'OP*_*'));
d = d([d.isdir]);
assert(~isempty(d), 'No OP*_* subfolders found in %s', resultsRoot);

fprintf('Found %d OP folders. Scanning fields...\n', numel(d));
allFieldNames = {};
validOPs = {};  % OP folder names that have a catalogue.mat

for k = 1:numel(d)
    catFile = fullfile(d(k).folder, d(k).name, 'catalogue.mat');
    if ~isfile(catFile), continue; end
    S = load(catFile); cat = pickCatalogue(S);
    if isempty(cat), continue; end
    allFieldNames = union(allFieldNames, fieldnames(cat), 'stable');
    validOPs{end+1} = d(k).name; %#ok<AGROW>
end
assert(~isempty(validOPs), 'No valid catalogue.mat files found.');
% Sort OPs numerically: OP1_0, OP2_0, ..., OP12_0, ..., OP78_0, OP78_1
opNums = cellfun(@opToNum, validOPs);
[~, sortIdx] = sort(opNums);
validOPs = validOPs(sortIdx);
fprintf('  %d catalogues, %d native fields.\n', numel(validOPs), numel(allFieldNames));

%% 3. GUI #1: OP selector (by phase + individual)
selectedOPs = opSelectorGUI(validOPs, phases, phaseOPs);
if isempty(selectedOPs)
    fprintf('Cancelled.\n'); T=[]; return
end
fprintf('Selected %d OPs.\n', numel(selectedOPs));

%% 4. Load selected catalogues
allRows = [];
for k = 1:numel(selectedOPs)
    opName  = selectedOPs{k};
    catFile = fullfile(resultsRoot, opName, 'catalogue.mat');
    S = load(catFile); cat = pickCatalogue(S);
    for j = 1:numel(cat), cat(j).OP = opName; end
    rows = catalogueToRows(cat, allFieldNames);
    allRows = [allRows; rows]; %#ok<AGROW>
    fprintf('  [ok] %s : %d events\n', opName, numel(cat));
end
T = struct2table(allRows);

%% 5. Add computed columns
T.OP_num = arrayfun(@(s) opToNum(char(s)), string(T.OP));
T.constellation = arrayfun(@(s) string(svToConst(char(s))), string(T.SV));
if all(ismember({'eventStart','eventEnd'}, T.Properties.VariableNames))
    T.duration_s = cellfun(@(a,b) durSecStr(a,b), T.eventStart, T.eventEnd);
end
if all(ismember({'segStart','segEnd'}, T.Properties.VariableNames))
    T.trackingDuration_s = cellfun(@trackDur, T.segStart, T.segEnd);
end

% Add mission phase column
T.missionPhase = repmat("",height(T),1);
for r = 1:height(T)
    op = char(T.OP(r));
    if isKey(opPhaseMap, op), T.missionPhase(r) = string(opPhaseMap(op)); end
end

% deltaSTEC
if all(ismember({'peakSTEC','minSTEC'}, T.Properties.VariableNames))
    T.deltaSTEC = cellfun(@(a,b) deltaCalc(a,b), ...
        num2cell(T.peakSTEC), num2cell(T.minSTEC));
end

%% 6. Sort
T = sortrows(T, {'OP_num','constellation','PRN'});
T.OP = arrayfun(@(x) string(numToOp(x)), T.OP_num);
T = movevars(T, 'OP', 'Before', 1);

%% 7. GUI #2: column selector per format
coreAlways = {'OP','constellation','SV','PRN'};
toggleable = setdiff(T.Properties.VariableNames, [coreAlways,{'OP_num'}], 'stable');
formats    = {'mat','csv','tex','docx'};
selections = columnSelectorGUI(toggleable, formats);
if isempty(selections)
    fprintf('Cancelled.\n'); T=[]; return
end

%% 8. Print statistics
printStatistics(T);
plotEventStatisticsInline(T);  

%% 9. Drop helper & export
T.OP_num = [];
basename  = fullfile(outDir, 'catalogue_all');
mergeCols = {'OP','constellation','SV','PRN'};

for f = 1:numel(formats)
    fmt  = formats{f};
    cols = unique([coreAlways, selections.(fmt)], 'stable');
    cols = intersect(cols, T.Properties.VariableNames, 'stable');
    Tsub = T(:, cols);
    switch fmt
        case 'mat'
            Tmat=Tsub; save([basename '.mat'],'Tmat'); %#ok<NASGU>
            fprintf('Saved %s.mat\n', basename);
        case 'csv'
            writetable(Tsub,[basename '.csv']);
            fprintf('Saved %s.csv\n', basename);
        case 'tex'
            writeLatex(Tsub,[basename '.tex'],mergeCols);
            fprintf('Saved %s.tex\n', basename);
        case 'docx'
            try
                writeDocx(Tsub,[basename '.docx'],mergeCols);
                fprintf('Saved %s.docx\n', basename);
            catch ME
                warning('DOCX failed: %s', ME.message);
            end
    end
end
fprintf('\nDone. %d events exported.\n', height(T));
end

% #########################################################################
%                   STATISTICS (dynamic event types)
% #########################################################################

% function printStatistics(T)
% fprintf('\n===== CLASSIFICATION STATISTICS =====\n');
% if ~ismember('eventType', T.Properties.VariableNames)
%     fprintf('  No eventType column found.\n'); return
% end
% 
% types = string(T.eventType);
% uTypes = unique(types(strlength(types)>0));
% nTotal = height(T);
% 
% fprintf('  Total events: %d\n', nTotal);
% fprintf('  %-20s  %5s  %6s\n', 'Event Type', 'Count', '%');
% fprintf('  %s\n', repmat('-',1,35));
% for i = 1:numel(uTypes)
%     n = sum(types == uTypes(i));
%     fprintf('  %-20s  %5d  %5.1f%%\n', uTypes(i), n, 100*n/nTotal);
% end
% 
% % Per-constellation breakdown
% if ismember('constellation', T.Properties.VariableNames)
%     consts = unique(string(T.constellation));
%     consts = consts(strlength(consts)>0);
%     fprintf('\n  Per-constellation breakdown:\n');
%     for c = 1:numel(consts)
%         mask = string(T.constellation)==consts(c);
%         nc = sum(mask);
%         fprintf('  %s (%d events):\n', consts(c), nc);
%         for i = 1:numel(uTypes)
%             n = sum(mask & types==uTypes(i));
%             if n>0, fprintf('    %-20s  %5d  %5.1f%%\n', uTypes(i), n, 100*n/nc); end
%         end
%     end
% end
% 
% % Per mission phase
% if ismember('missionPhase', T.Properties.VariableNames)
%     mps = unique(string(T.missionPhase));
%     mps = mps(strlength(mps)>0);
%     if ~isempty(mps)
%         fprintf('\n  Per-phase breakdown:\n');
%         for m = 1:numel(mps)
%             mask = string(T.missionPhase)==mps(m);
%             nm = sum(mask);
%             fprintf('  %s (%d events):\n', mps(m), nm);
%             for i = 1:numel(uTypes)
%                 n = sum(mask & types==uTypes(i));
%                 if n>0, fprintf('    %-20s  %5d  %5.1f%%\n', uTypes(i), n, 100*n/nm); end
%             end
%         end
%     end
% end
% fprintf('=====================================\n\n');
% end

function printStatistics(T)
fprintf('\n===== CLASSIFICATION STATISTICS =====\n');
if ~ismember('eventType', T.Properties.VariableNames)
    fprintf('  No eventType column found.\n'); return
end

types = string(T.eventType);
uTypes = unique(types(strlength(types)>0));
nTotal = height(T);

hasDur = ismember('duration_s', T.Properties.VariableNames);
if hasDur, dur = T.duration_s; end

% --- Overall ---
fprintf('  Total events: %d\n', nTotal);
if hasDur
    fprintf('  %-20s  %5s  %6s  %8s\n', 'Event Type', 'Count', '%', 'mean t[s]');
else
    fprintf('  %-20s  %5s  %6s\n', 'Event Type', 'Count', '%');
end
fprintf('  %s\n', repmat('-',1,50));
for i = 1:numel(uTypes)
    mask = types == uTypes(i);
    n = sum(mask);
    if hasDur
        mt = round(mean(dur(mask), 'omitnan'));
        fprintf('  %-20s  %5d  %5.1f%%  %8d\n', uTypes(i), n, 100*n/nTotal, mt);
    else
        fprintf('  %-20s  %5d  %5.1f%%\n', uTypes(i), n, 100*n/nTotal);
    end
end

% --- Per-constellation ---
if ismember('constellation', T.Properties.VariableNames)
    consts = unique(string(T.constellation));
    consts = consts(strlength(consts)>0);
    fprintf('\n  Per-constellation breakdown:\n');
    for c = 1:numel(consts)
        maskC = string(T.constellation)==consts(c);
        nc = sum(maskC);
        fprintf('  %s (%d events):\n', consts(c), nc);
        for i = 1:numel(uTypes)
            mask = maskC & types==uTypes(i);
            n = sum(mask);
            if n>0
                if hasDur
                    mt = round(mean(dur(mask), 'omitnan'));
                    fprintf('    %-20s  %5d  %5.1f%%  %8d\n', uTypes(i), n, 100*n/nc, mt);
                else
                    fprintf('    %-20s  %5d  %5.1f%%\n', uTypes(i), n, 100*n/nc);
                end
            end
        end
    end
end

% --- Per mission phase ---
if ismember('missionPhase', T.Properties.VariableNames)
    mps = unique(string(T.missionPhase));
    mps = mps(strlength(mps)>0);
    if ~isempty(mps)
        fprintf('\n  Per-phase breakdown:\n');
        for m = 1:numel(mps)
            maskP = string(T.missionPhase)==mps(m);
            nm = sum(maskP);
            fprintf('  %s (%d events):\n', mps(m), nm);
            for i = 1:numel(uTypes)
                mask = maskP & types==uTypes(i);
                n = sum(mask);
                if n>0
                    if hasDur
                        mt = round(mean(dur(mask), 'omitnan'));
                        fprintf('    %-20s  %5d  %5.1f%%  %8d\n', uTypes(i), n, 100*n/nm, mt);
                    else
                        fprintf('    %-20s  %5d  %5.1f%%\n', uTypes(i), n, 100*n/nm);
                    end
                end
            end
        end
    end
end
fprintf('=====================================\n\n');
end

function plotEventStatisticsInline(T)
%PLOTEVENTSTATISTICSINLINE  Called from exportCatalogues after table is built.
%  Identical logic to the standalone plotEventStatistics.m but receives
%  the table T directly instead of loading it from disk.

%% Sanitise
T.OP           = string(T.OP);
T.eventType    = lower(string(T.eventType));
T.constellation= upper(string(T.constellation));
T.missionPhase = upper(string(T.missionPhase));

%% Ordered OP list
allOPs  = unique(T.OP, 'stable');
opNums  = arrayfun(@opToNum, allOPs);   % opToNum already exists in exportCatalogues
[~, si] = sort(opNums);
allOPs  = allOPs(si);
nOPs    = numel(allOPs);

%% Phase per OP
phaseOrder = ["COMMISSIONING","TRANSIT","LLO","SURFACE"];
opPhase    = strings(nOPs,1);
for k = 1:nOPs
    mask = T.OP == allOPs(k);
    ph   = T.missionPhase(mask);
    ph   = ph(strlength(ph)>0);
    if isempty(ph), opPhase(k)="UNKNOWN"; else, [~,~,ic] = unique(ph); opPhase(k) = ph(mode(ic)); end
end

%% Sort OPs: phase order first, numeric within phase
phIdx = zeros(nOPs,1);
for k = 1:nOPs
    idx = find(phaseOrder == opPhase(k),1);
    if isempty(idx), idx = numel(phaseOrder)+1; end
    phIdx(k) = idx;
end
[~, fs]  = sortrows([phIdx, opNums(:)]);
allOPs   = allOPs(fs);
opPhase  = opPhase(fs);
opNums   = opNums(fs); %#ok<NASGU>

%% Count matrices
eventTypes = ["egress","ingress","tangent"];
consts     = ["GAL","GPS"];
nET = numel(eventTypes); nC = numel(consts);
countET = zeros(nOPs,nET); countC = zeros(nOPs,nC);
for k = 1:nOPs
    mask = T.OP == allOPs(k);
    for e = 1:nET, countET(k,e)=sum(mask & T.eventType==eventTypes(e)); end
    for c = 1:nC,  countC(k,c) =sum(mask & T.constellation==consts(c));  end
end

%% Phase separator positions & labels
sepPos=[]; phaseLabPos=[]; phaseLabels={}; currentPhase="";
for k = 1:nOPs
    if opPhase(k) ~= currentPhase
        if k>1, sepPos(end+1)=k-0.5; end         %#ok<AGROW>
        phaseLabPos(end+1)=k;                      %#ok<AGROW>
        phaseLabels{end+1}=char(opPhase(k));       %#ok<AGROW>
        currentPhase = opPhase(k);
    end
end

%% x-tick labels
xLabels = strrep(strrep(allOPs,"OP",""),"_","\_");

%% Colours
colET = [0.20 0.53 0.74; 0.89 0.47 0.20; 0.42 0.68 0.35];
colC  = [0.84 0.19 0.37; 0.20 0.53 0.74];

%% Helper: draw one figure
    %function drawFig(figNum, counts, catLabels, colours, titleStr)
        %fig = figure(figNum);
function drawFig(counts, catLabels, colours, titleStr)
    figNum =301;
        fig = figure;
        clf(fig);
        fig.Name  = titleStr;
        fig.Color = 'w';
        fig.Units = 'centimeters';
        fig.Position = [2, 10-(figNum-1)*10, 28, 9];

        ax = axes(fig);
        nCat = size(counts,2);
        bh   = bar(ax, 1:nOPs, counts, 'stacked', 'BarWidth', 0.75);
        for ci = 1:nCat
            bh(ci).FaceColor  = colours(ci,:);
            bh(ci).EdgeColor  = 'none';
            bh(ci).DisplayName= char(catLabels(ci));
        end

        yTop = max(sum(counts,2)) * 1.05;
        for s = 1:numel(sepPos)
            xline(ax, sepPos(s), '--k', 'LineWidth',0.8, 'HandleVisibility','off');
        end
        for p = 1:numel(phaseLabPos)
            if p<numel(phaseLabPos), xEnd=phaseLabPos(p+1)-1; else, xEnd=nOPs; end
            xMid=(phaseLabPos(p)+xEnd)/2;
            text(ax,xMid,yTop,phaseLabels{p}, ...
                'HorizontalAlignment','center','FontSize',8,'FontWeight','bold', ...
                'Interpreter','none');
        end

        ax.XTick               = 1:nOPs;
        ax.XTickLabel          = xLabels;
        ax.XTickLabelRotation  = 45;
        ax.TickLabelInterpreter= 'tex';
        ax.FontSize            = 8;
        ax.YLabel.String       = 'Number of events';
        ax.YLabel.FontSize     = 10;
        ax.XLabel.String       = 'Operation ID';
        ax.XLabel.FontSize     = 10;
        ax.Title.String        = titleStr;
        ax.Title.FontSize      = 11;
        ax.Color               = [0.97 0.97 0.97];
        ax.Layer               = 'top';
        grid(ax,'on');
        legend(ax,'Location','northeast','FontSize',9,'Box','off');
        ylim(ax,[0, yTop*1.15]);
    end

    % drawFig(301, countET, eventTypes, colET, ...
    %     'Classified events per operation – breakdown by event type');
    % drawFig(302, countC,  consts,     colC,  ...
    %     'Classified events per operation – breakdown by constellation');
    drawFig(countET, eventTypes, colET, ...
        'Classified events per operation – breakdown by event type');
    drawFig(countC,  consts,     colC,  ...
        'Classified events per operation – breakdown by constellation');
    
    %% deltaSTEC histograms
    if ismember('deltaSTEC', T.Properties.VariableNames)
        dSTEC  = T.deltaSTEC;
        validD = ~isnan(dSTEC);
        
        histSets = {"Total",   validD; ...
                    "egress",  validD & T.eventType=="egress"; ...
                    "ingress", validD & T.eventType=="ingress"; ...
                    "tangent", validD & T.eventType=="tangent"};
        
        colH = [0.50 0.50 0.50;
                0.20 0.53 0.74;
                0.89 0.47 0.20;
                0.42 0.68 0.35];
        
        for h = 1:size(histSets,1)
            lbl  = histSets{h,1};
            mask = histSets{h,2};
            vals = dSTEC(mask);
            if isempty(vals), continue; end
            
            fig = figure;
            fig.Name  = sprintf('deltaSTEC – %s', lbl);
            fig.Color = 'w';
            fig.Units = 'centimeters';
            fig.Position = [2+5*(h-1), 2, 16, 10];
            
            ax = axes(fig);
            histogram(ax, vals, 'BinMethod', 'fd', ...
                'FaceColor', colH(h,:), 'EdgeColor', 'w', 'FaceAlpha', 0.85);
            
            ax.TickLabelInterpreter = 'latex';
            ax.FontSize = 14;
            title(ax, sprintf('$\\Delta$STEC distribution -- %s ($N$=%d)', lbl, numel(vals)), ...
                'Interpreter', 'latex', 'FontSize', 14);
            xlabel(ax, '$\Delta$STEC [TECu]', 'Interpreter', 'latex', 'FontSize', 14);
            ylabel(ax, 'Count', 'Interpreter', 'latex', 'FontSize', 14);
            ax.Color = [0.97 0.97 0.97];
            ax.Layer = 'top';
            grid(ax, 'on');
            
            % Mean and median as legend entries
            mu  = mean(vals, 'omitnan');
            med = median(vals, 'omitnan');
            
            xline(ax, mu,  '-r',  'LineWidth', 1.2, 'HandleVisibility', 'on', ...
                'DisplayName', sprintf('Mean = %.1f TECu', mu));
            xline(ax, med, '--b', 'LineWidth', 1.2, 'HandleVisibility', 'on', ...
                'DisplayName', sprintf('Median = %.1f TECu', med));
            
            legend(ax, 'Interpreter', 'latex', 'FontSize', 12, ...
                'Location', 'northeast', 'Box', 'off');
        end
    end

end
% #########################################################################
%                   GUI #1: OP SELECTOR (by phase + individual)
% #########################################################################
function selectedOPs = opSelectorGUI(validOPs, phases, phaseOPs)

nOP = numel(validOPs);
rowH = 20;

% Left panel: phase switches + ALL; Right panel: individual OPs (scrollable)
leftW  = 220;
rightW = 250;
figH = max(400, nOP*rowH + 150);
figW   = leftW + rightW + 20;

fig = figure('Name','OP Selector','NumberTitle','off', ...
    'MenuBar','none','ToolBar','none','Resize','off', ...
    'Position',[300,150,figW,figH]);
fig.UserData.cancelled = false;
fig.CloseRequestFcn = @(src,~) closeCancelCb(src);

% --- Left panel: ALL + per-phase ---
yTop = figH - 30;
allSwitch = uicontrol(fig,'Style','checkbox','String','ALL OPs', ...
    'Value',1,'FontWeight','bold','FontSize',11, ...
    'Position',[15,yTop,190,22]);

phaseCb = {};
phaseNames = {};
y = yTop - 30;
for p = 1:numel(phases)
    pName = phases{p};
    if isfield(phaseOPs, pName)
        phaseNames{end+1} = pName; %#ok<AGROW>
        phaseCb{end+1} = uicontrol(fig,'Style','checkbox', ...
            'String',pName,'Value',1,'FontSize',10, ...
            'Position',[30,y,180,22]); %#ok<AGROW>
        y = y - 25;
    end
end

% --- Right panel: individual OPs (in a scrollable panel) ---
panelH = figH - 80 - 20;
scrollPanel = uipanel(fig,'Units','pixels', ...
    'Position',[leftW, 50, rightW, figH-80], ...
    'Title','Individual OPs');

opCb = gobjects(nOP,1);
for k = 1:nOP
    yPos = panelH - k*rowH;
    % Add phase tag if known
    lbl = validOPs{k};
    opCb(k) = uicontrol(scrollPanel,'Style','checkbox', ...
        'String',lbl,'Value',1, ...
        'Position',[10, yPos, 220, 20]);
end

% --- Assign callbacks after all created ---
allSwitch.Callback = @(src,~) toggleAllOPs(src, phaseCb, opCb, validOPs, phaseNames, phaseOPs);
for p = 1:numel(phaseCb)
    pName = phaseNames{p};
    phaseCb{p}.Callback = @(src,~) togglePhaseOPs(src, pName, phaseOPs, opCb, validOPs, allSwitch, phaseCb);
end
for k = 1:nOP
    opCb(k).Callback = @(~,~) syncPhaseAndAll(allSwitch, phaseCb, opCb, validOPs, phaseNames, phaseOPs);
end

% --- OK / Cancel ---
uicontrol(fig,'Style','pushbutton','String','OK', ...
    'Position',[figW/2-110,10,100,30],'Callback',@(~,~)uiresume(fig));
uicontrol(fig,'Style','pushbutton','String','Cancel', ...
    'Position',[figW/2+10,10,100,30],'Callback',@(~,~)closeCancelCb(fig));

uiwait(fig);
if ~isvalid(fig)||fig.UserData.cancelled
    selectedOPs=[]; if isvalid(fig), delete(fig); end; return
end

selectedOPs = {};
for k = 1:nOP
    if opCb(k).Value, selectedOPs{end+1} = validOPs{k}; end %#ok<AGROW>
end
delete(fig);
end

function toggleAllOPs(src, phaseCbs, opCbs, ~, ~, ~)
v = src.Value;
for i=1:numel(phaseCbs), phaseCbs{i}.Value=v; end
for i=1:numel(opCbs), opCbs(i).Value=v; end
end

function togglePhaseOPs(src, pName, phaseOPs, opCbs, validOPs, allCb, phaseCbs)
v = src.Value;
if isfield(phaseOPs, pName)
    members = phaseOPs.(pName);
    for k = 1:numel(validOPs)
        if ismember(validOPs{k}, members)
            opCbs(k).Value = v;
        end
    end
end
syncAllCb(allCb, phaseCbs, opCbs);
end

function syncPhaseAndAll(allCb, phaseCbs, opCbs, validOPs, phaseNames, phaseOPs)
% Sync phase checkboxes based on individual OP states
for p = 1:numel(phaseNames)
    pName = phaseNames{p};
    if isfield(phaseOPs, pName)
        members = phaseOPs.(pName);
        allOn = true;
        for k = 1:numel(validOPs)
            if ismember(validOPs{k}, members) && ~opCbs(k).Value
                allOn = false; break
            end
        end
        phaseCbs{p}.Value = allOn;
    end
end
syncAllCb(allCb, phaseCbs, opCbs);
end

function syncAllCb(allCb, ~, opCbs)
allCb.Value = all(arrayfun(@(c)c.Value, opCbs));
end

function closeCancelCb(src)
src.UserData.cancelled=true; uiresume(src);
end

% #########################################################################
%                   GUI #2: COLUMN SELECTOR
% #########################################################################
function selections = columnSelectorGUI(toggleable, formats)
nT=numel(toggleable); nF=numel(formats);
rowH=22; colW=200; padT=60;
figH=padT+nT*rowH+80; figW=colW*nF+40;

fig=figure('Name','Column Selector','NumberTitle','off', ...
    'MenuBar','none','ToolBar','none','Resize','off', ...
    'Position',[300,200,figW,figH]);
fig.UserData.cancelled=false;
fig.CloseRequestFcn=@(src,~)closeCancelCb(src);

cb=cell(1,nF); allCb=cell(1,nF);
for f=1:nF
    xOff=20+(f-1)*colW;
    uicontrol(fig,'Style','text','String',upper(formats{f}), ...
        'FontWeight','bold','FontSize',12, ...
        'Position',[xOff,figH-35,160,25],'HorizontalAlignment','center');
    allCb{f}=uicontrol(fig,'Style','checkbox','String','ALL', ...
        'Value',1,'FontWeight','bold','Position',[xOff,figH-55,160,20]);
    cb{f}=gobjects(nT,1);
    for c=1:nT
        cb{f}(c)=uicontrol(fig,'Style','checkbox','String',toggleable{c}, ...
            'Value',1,'Position',[xOff,figH-padT-c*rowH,180,20]);
    end
end
for f=1:nF
    allCb{f}.Callback=@(src,~)toggleAllCols(src,cb{f});
    for c=1:nT
        cb{f}(c).Callback=@(~,~)syncAllCol(allCb{f},cb{f});
    end
end

uicontrol(fig,'Style','pushbutton','String','OK', ...
    'Position',[figW/2-110,30,100,50],'Callback',@(~,~)uiresume(fig));
uicontrol(fig,'Style','pushbutton','String','Cancel', ...
    'Position',[figW/2+10,30,100,50],'Callback',@(~,~)closeCancelCb(fig));

uiwait(fig);
if ~isvalid(fig)||fig.UserData.cancelled
    selections=[]; if isvalid(fig),delete(fig);end; return
end
selections=struct();
for f=1:nF
    sel={};
    for c=1:nT, if cb{f}(c).Value, sel{end+1}=toggleable{c}; end; end %#ok<AGROW>
    selections.(formats{f})=sel;
end
delete(fig);
end

function toggleAllCols(src,cbs)
for i=1:numel(cbs),cbs(i).Value=src.Value;end
end
function syncAllCol(allCb,cbs)
allCb.Value=all(arrayfun(@(c)c.Value,cbs));
end

% #########################################################################
%                       CATALOGUE LOADING
% #########################################################################
function cat = pickCatalogue(S)
fn=fieldnames(S);
for c={'catalogue','event_catalogue','Catalogue','events'}
    if ismember(c{1},fn), cat=S.(c{1}); return; end
end
for i=1:numel(fn), if isstruct(S.(fn{i})), cat=S.(fn{i}); return; end; end
cat=[];
end

function rows = catalogueToRows(cat, masterFields)
rows=repmat(struct(),numel(cat),1);
for i=1:numel(cat)
    e=cat(i);
    for f=1:numel(masterFields)
        fn=masterFields{f};
        rows(i).(fn)=convertField(getf(e,fn,[]));
    end
end
end

function v = convertField(raw)
if isempty(raw), v=''; return; end
if isdatetime(raw)
    if isnat(raw),v='';else,v=char(datetime(raw,'Format','yyyy-MM-dd HH:mm:ss'));end; return
end
if isnumeric(raw)
    if isscalar(raw),v=double(raw);
    else, v=['[' strjoin(arrayfun(@(x)num2str(x),raw(:).','UniformOutput',false),',') ']'];
    end; return
end
try,v=char(string(raw));catch,v='';end
end

function v=getf(s,f,d), if isfield(s,f)&&~isempty(s.(f)),v=s.(f);else,v=d;end, end

function d=durSecStr(t1,t2)
if ischar(t1)&&ischar(t2)&&~isempty(t1)&&~isempty(t2)
    try,d=round(seconds(datetime(t2)-datetime(t1)));catch,d=NaN;end
else,d=NaN;end
end

function d=trackDur(segS,segE)
try
    if ischar(segS),segS=str2num(segS);end %#ok<ST2NM>
    if ischar(segE),segE=str2num(segE);end %#ok<ST2NM>
    if isnumeric(segS)&&isnumeric(segE)&&numel(segS)==numel(segE)
        d=sum(segE-segS);
    else,d=NaN;end
catch,d=NaN;end
end

function d=deltaCalc(a,b)
if isnumeric(a)&&isnumeric(b)&&~isnan(a)&&~isnan(b), d=a-b; else, d=NaN; end
end

% #########################################################################
%                        OP / CONSTELLATION HELPERS
% #########################################################################
function n=opToNum(op)
if iscell(op),op=op{1};end; op=strtrim(char(op));
op=regexprep(op,'^[Oo][Pp]',''); parts=strsplit(op,'_');
major=str2double(parts{1}); if isnan(major),n=NaN;return;end
if numel(parts)>=2&&~isempty(parts{2})
    n=major+str2double(parts{2})/10^numel(parts{2});
else,n=major;end
end

function s=numToOp(n)
if isnan(n),s='';return;end; major=floor(n); frac=n-major;
if frac==0, s=sprintf('OP%d_0',major);
else
    str=sprintf('%.10g',frac); dotIdx=strfind(str,'.');
    m=regexprep(str(dotIdx+1:end),'0+$','');
    if isempty(m),m='0';end; s=sprintf('OP%d_%s',major,m);
end
end

function c=svToConst(sv)
if iscell(sv),sv=sv{1};end; sv=char(sv);
if isempty(sv),c='';return;end
if numel(sv)>=3&&any(strcmpi(sv(1:3),{'GPS','GAL','GLO','BDS','QZS','IRN'}))
    c=upper(sv(1:3));return;end
switch upper(sv(1))
    case 'G',c='GPS';case 'E',c='GAL';case 'R',c='GLO';
    case 'C',c='BDS';case 'J',c='QZS';case 'I',c='IRN';
    otherwise,c=upper(sv(1));
end
end

% #########################################################################
%                        MERGE SPAN HELPERS
% #########################################################################
function spans=computeMergeSpans(vals)
N=numel(vals);spans=ones(N,1);r=1;
while r<=N, k=1;
    while r+k<=N&&strcmp(vals(r+k),vals(r)),k=k+1;end
    if k>1,spans(r)=k;spans(r+1:r+k-1)=0;end; r=r+k;
end; end

function spans=splitSpansAtPages(spans,pageRows)
N=numel(spans);r=1;
while r<=N, sp=spans(r);if sp<=1,r=r+1;continue;end
    rowOnPage=mod(r-1,pageRows);roomLeft=pageRows-rowOnPage;
    if sp>roomLeft&&roomLeft>0
        spans(r)=roomLeft;nS=r+roomLeft;nSp=sp-roomLeft;
        spans(r+1:nS-1)=0;spans(nS)=nSp;spans(nS+1:r+sp-1)=0;r=nS;
    else,r=r+sp;end
end; end

function sp=getSpan(col,row,mergeMap,spanArrays)
if isKey(mergeMap,col),sp=spanArrays{mergeMap(col)}(row);else,sp=1;end
end

function [mergeMap,spanArrays,mergeColIdx,firstDataCol]=buildMergeSpans(T,cols,mergeCols,pageRows)
mergeColIdx=[];
for i=1:numel(mergeCols)
    idx=find(strcmp(cols,mergeCols{i}));
    if ~isempty(idx),mergeColIdx(end+1)=idx;end %#ok<AGROW>
end
mergeColIdx=sort(mergeColIdx);
if isempty(mergeColIdx),firstDataCol=1;mergeMap=containers.Map();spanArrays={};return;end
firstDataCol=max(mergeColIdx)+1;
nRows=height(T);compositeKey=strings(nRows,1);
spanArrays={};mergeMap=containers.Map();
for m=1:numel(mergeColIdx)
    ci=mergeColIdx(m);colName=cols{ci};
    colVals=strings(nRows,1);
    for r=1:nRows,colVals(r)=string(formatCell(T.(colName)(r),colName));end
    compositeKey=compositeKey+"|"+colVals;
    sp=computeMergeSpans(compositeKey);
    if ~isinf(pageRows),sp=splitSpansAtPages(sp,pageRows);end
    spanArrays{end+1}=sp;mergeMap(colName)=numel(spanArrays); %#ok<AGROW>
end; end

% #########################################################################
%                        FORMATTING HELPERS
% #########################################################################
function s=formatCell(x,name)
if iscell(x),x=x{1};end
if isnumeric(x)&&isscalar(x)
    if isnan(x),s='--';
    elseif any(strcmp(name,{'PRN','duration_s','trackingDuration_s'})),s=sprintf('%d',round(x));
    else,s=sprintf('%.4g',x);end
elseif ischar(x)||isstring(x), s=char(x); if isempty(s),s='--';end
else,s='--';end
end

function s=texEscape(s)
s=strrep(s,'\','\textbackslash{}');
s=strrep(s,'&','\&');s=strrep(s,'%','\%');s=strrep(s,'$','\$');
s=strrep(s,'#','\#');s=strrep(s,'_','\_');s=strrep(s,'{','\{');
s=strrep(s,'}','\}');s=strrep(s,'~','\textasciitilde{}');s=strrep(s,'^','\textasciicircum{}');
end

function h=prettyHeader(col)
map=containers.Map( ...
    {'OP','constellation','SV','PRN','eventType','eventTypeGeom', ...
     'confidence','channel','eventStart','eventEnd','duration_s', ...
     'trackingDuration_s','missionPhase', ...
     'peakROTI','minROTI','meanROT','peakSTEC','minSTEC', 'deltaSTEC', ...
     'grazingAltStart','grazingAltEnd','segStart','segEnd','peakIdx','minIdx'}, ...
    {'OP','Const.','SV','PRN','Type','Type Geom.', ...
     'Conf.','Channel','Start (UTC)','End (UTC)','Dur. [s]', ...
     'Track. Dur. [s]','Phase', ...
     'peak ROTI [TECu]','min ROTI [TECu]','mean ROT [TECu/s]', ...
     'peak STEC [TECu]','min STEC [TECu]','$\Delta$STEC [TECu]', ...
     'Alt. Start [km]','Alt. End [km]','Seg. Start','Seg. End','Peak Idx','Min Idx'});
if isKey(map,col),h=map(col);else,h=col;end
end

function w=colWidth(col)
map=containers.Map( ...
    {'OP','constellation','SV','PRN','eventType','eventTypeGeom', ...
     'confidence','channel','eventStart','eventEnd','duration_s', ...
     'trackingDuration_s','missionPhase', ...
     'peakROTI','minROTI','meanROT','peakSTEC','minSTEC','deltaSTEC' ...
     'grazingAltStart','grazingAltEnd','segStart','segEnd','peakIdx','minIdx'}, ...
    {'1.1','0.9','0.7','0.5','1.2','1.2', ...
     '0.6','0.9','2.5','2.5','0.6', ...
     '0.7','1.3', ...
     '1.0','1.0','1.0','1.0','1.0','1.1', ...
     '1.1','1.1','1.0','1.0','0.7','0.7'});
if isKey(map,col),w=map(col);else,w='1.2';end
end

% #########################################################################
%                        LATEX EXPORT
% #########################################################################
function writeLatex(T,fname,~)
cols=T.Properties.VariableNames; nCols=numel(cols); nRows=height(T);
widthSpec=cellfun(@(c)sprintf('m{%scm}',colWidth(c)),cols,'UniformOutput',false);
colSpec=['|' strjoin(widthSpec,'|') '|'];
hdr=cellfun(@prettyHeader,cols,'UniformOutput',false);

fid=fopen(fname,'w'); cObj=onCleanup(@()fclose(fid)); %#ok<NASGU>
fprintf(fid,'%% Auto-generated by exportCatalogues.m v4\n');
fprintf(fid,'%% USAGE: \\input{catalogue_all.tex}\n');
fprintf(fid,'%%   \\usepackage{longtable, array, lscape}\n');
fprintf(fid,'%%   \\usepackage[table]{xcolor}\n\n');

fprintf(fid,'\\begin{landscape}\n{\\footnotesize\\centering\n');
fprintf(fid,'\\begin{longtable}{%s}\n',colSpec);
fprintf(fid,'\\caption{Aggregated event catalogue.}\\label{tab:catalogue}\\\\\n');

hdrCells=cellfun(@(h)sprintf('\\centering\\arraybackslash\\textbf{%s}',texEscape(h)),hdr,'UniformOutput',false);
hdrLine=strjoin(hdrCells,' & ');
for rep=1:2
    fprintf(fid,'\\hline\n\\rowcolor[gray]{0.85}\n');
    fprintf(fid,'%s \\\\\\hline\n',hdrLine);
    if rep==1,fprintf(fid,'\\endfirsthead\n');else,fprintf(fid,'\\endhead\n');end
end
fprintf(fid,'\\hline\n\\endfoot\n');

for r=1:nRows
    cellsRow=cell(1,nCols);
    for i=1:nCols
        cellsRow{i}=sprintf('\\centering %s',texEscape(formatCell(T.(cols{i})(r),cols{i})));
    end
    fprintf(fid,'%s \\\\\\hline\n',strjoin(cellsRow,' & '));
end
fprintf(fid,'\\end{longtable}\n}\\end{landscape}\n');
end

% #########################################################################
%                   DOCX EXPORT (MATLAB Report Generator)
% #########################################################################
function writeDocx(T,fname,~)
import mlreportgen.dom.*
cols=T.Properties.VariableNames; nCols=numel(cols); nRows=height(T);
hdr=cellfun(@prettyHeader,cols,'UniformOutput',false);
tS={FontFamily('Times New Roman'),FontSize('5pt')};

bodyData=cell(nRows+1,nCols);
for i=1:nCols, bodyData{1,i}=hdr{i}; end
for r=1:nRows, for i=1:nCols
    bodyData{r+1,i}=formatCell(T.(cols{i})(r),cols{i});
end; end

doc=Document(fname,'docx');
sp2=doc.CurrentPageLayout; sp2.PageSize.Orientation='landscape';
sp2.PageSize.Height='8.5in'; sp2.PageSize.Width='11in';
h1=Heading1('Aggregated Event Catalogue');
h1.Style={FontFamily('Times New Roman'),FontSize('14pt'),Bold(true)};
append(doc,h1);
p1=Paragraph(sprintf('Total events: %d',nRows));
p1.Style={FontFamily('Times New Roman'),FontSize('10pt')};
append(doc,p1);

tbl=Table(bodyData);
tbl.Style=[tS,{Border('solid','black','0.5pt'),RowSep('solid','black','0.5pt'), ...
    ColSep('solid','black','0.5pt'),Width('100%')}];

hR=tbl.Children(1);
for i=1:numel(hR.Children),hR.Children(i).Style=[hR.Children(i).Style,tS,{Bold(true)}];end
for r=1:nRows
    rO=tbl.Children(r+1);
    for i=1:nCols, rO.Children(i).Style=[rO.Children(i).Style,tS]; end
end

append(doc,tbl); close(doc);
end