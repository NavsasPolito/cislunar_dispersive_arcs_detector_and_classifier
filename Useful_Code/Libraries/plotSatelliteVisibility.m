function plotSatelliteVisibility(obsFile,settings)
% plotSatelliteVisibility  GNSS satellite visibility scatter (L1 vs L5)
%
%   plotSatelliteVisibility(obsFile)
%
%   Per ogni epoca in cui il satellite ha una misura valida (rawP non NaN)
%   disegna un pallino. Per lo stesso PRN vengono disegnate due righe di
%   pallini leggermente sfalsate in verticale: una per L1 (G1C / E1B) e
%   una per L5 (G5I / E5I), con colori diversi.
%
%   Il nome della figura viene generato automaticamente come
%       "visibility plot_<OPxx_y>"
%   estraendo dal path del file il token OP... (es. OP1_0) che deve
%   comparire come componente di cartella.
%
%       link = '/percorso/OP1_0/obsSolutions.mat';
%       plotSatelliteVisibility(link);
% -------------------------------------------------------------------------
%   AUTHOR:
%       Francesco Fiorina, NavSAS - DET - Politecnico di Torino, 13/04/2026
%
%   VERSION:
%       1.0
%
%   COPYRIGHT:
%       NavSAS - DET - Politecnico di Torino, 2026. All rights reserved.
% -------------------------------------------------------------------------
    
    if ~settings.PlotVisibility
        return
    end

    if nargin < 1 || isempty(obsFile)
        error('Devi passare il path del file obsSolutions.mat');
    end

    if isempty(settings.OP)
        opID = [];
    else
        opID = settings.OP;
    end

    if ischar(obsFile)
        % ---------- estrai nome operazione dal path ----------
        % cerca un token tipo OP<digits>(_<digits>)? tra separatori di path
        tok = regexp(obsFile, '[\\/](OP\d+(?:_\d+)?)[\\/]', 'tokens', 'once');
        if isempty(tok)
            opName = 'unknown';
        else
            opName = tok{1};
        end
    
        % ---------- load ----------
        S = load(obsFile);
        if ~isfield(S,'obsSolutions')
            error('Il file "%s" non contiene la variabile obsSolutions.', obsFile);
        end

        obs = S.obsSolutions;

    elseif isstruct(obsFile)
        obs = obsFile;
        opName = opID;
    else
        error('Format file not accepted')
    end
    
    N   = numel(obs);

    % ---------- metadati di ogni traccia ----------
    SVs      = strings(N,1);
    consName = strings(N,1);
    PRNs     = zeros(N,1);
    band     = strings(N,1);     % "L1" o "L5"

    for k = 1:N
        SVs(k)  = string(obs(k).SV);
        PRNs(k) = double(obs(k).PRN(1));

        switch extractBefore(SVs(k),2)
            case "G", consName(k) = "GPS";
            case "E", consName(k) = "Galileo";
            case "R", consName(k) = "GLONASS";
            case "C", consName(k) = "BeiDou";
            case "J", consName(k) = "QZSS";
            case "I", consName(k) = "IRNSS";
            otherwise, consName(k) = extractBefore(SVs(k),2);
        end

        ch = upper(string(obs(k).channel));
        if     contains(ch,"1"), band(k) = "L1";
        elseif contains(ch,"5"), band(k) = "L5";
        else,                    band(k) = ch;
        end
    end

    % ---------- asse Y: una riga per (costellazione, PRN) ----------
    satRows         = table(consName, PRNs, 'VariableNames',{'cons','PRN'});
    uniqueSats      = unique(satRows,'rows');
    uniqueSats      = sortrows(uniqueSats, {'cons','PRN'});
    uniqueSats.yRow = (1:height(uniqueSats))';

    yPos = zeros(N,1);
    for k = 1:N
        yPos(k) = uniqueSats.yRow( ...
            uniqueSats.cons == consName(k) & uniqueSats.PRN == PRNs(k));
    end

    % ---------- figura ----------
    figTitle = sprintf('visibility plot_%s', opName);
    fig = figure('Color','w','Position',[100 100 1200 700], ...
                 'Name', figTitle, 'NumberTitle','off');
    ax  = axes(fig); hold(ax,'on'); grid(ax,'on'); grid(ax,'minor'); box(ax,'on');
    ax.GridAlpha   = 0.25;
    ax.TickLabelInterpreter = 'latex';

    offL1 =  0.18;
    offL5 = -0.18;
    colL1 = [0.12 0.47 0.95];
    colL5 = [0.95 0.45 0.10];
    mSize = 22;

    leapSec  = 18;
    gpsEpoch = datetime(1980,1,6,0,0,0,'TimeZone','UTC');

    % prima passata: calcola i timestamps validi di ogni traccia
    tValidCell = cell(N,1);
    for k = 1:N
        e      = obs(k);
        rawP   = e.rawP(:);
        gpsSec = double(e.GPStime(:));
        wn     = double(e.weekNumber(:));
        nEp    = min(numel(rawP), numel(gpsSec));

        valid  = ~isnan(rawP(1:nEp)) & ~isnan(gpsSec(1:nEp));
        if ~any(valid)
            tValidCell{k} = NaT(0,1,'TimeZone','UTC');
            continue;
        end
        tUTC = gpsEpoch + seconds(wn.*7*86400 + gpsSec(1:nEp) - leapSec);
        tValidCell{k} = tUTC(valid);
    end

    % ---- banda gialla: dove L1 e L5 sono entrambi presenti ----
    % la disegno PRIMA degli scatter cosi' resta sotto
    colDual    = [1.0 0.92 0.20];   % giallo
    dualAlpha  = 0.30;              % trasparenza
    bandHeight = 0.70;              % altezza verticale della banda per satellite

    for j = 1:height(uniqueSats)
        idxL1 = find(band == "L1" & yPos == j);
        idxL5 = find(band == "L5" & yPos == j);
        if isempty(idxL1) || isempty(idxL5); continue; end

        tL1 = vertcat(tValidCell{idxL1});
        tL5 = vertcat(tValidCell{idxL5});
        if isempty(tL1) || isempty(tL5); continue; end

        % intersezione a livello di epoca (tolleranza 0.5 s)
        [~, iA] = ismembertol(datenum(tL1), datenum(tL5), ...
                              0.5/86400, 'DataScale', 1);

        tDual = tL1(iA > 0);
        if isempty(tDual); continue; end

        % segmenta in run contigui (gap > 2 s => nuovo segmento)
        tDual  = sort(tDual);
        gaps   = [true; seconds(diff(tDual)) > 2];
        segId  = cumsum(gaps);
        for s = 1:max(segId)
            ts = tDual(segId == s);
            if numel(ts) < 2; continue; end
            x0 = ts(1); x1 = ts(end);
            yc = j;
            patch(ax, [x0 x1 x1 x0], ...
                      [yc-bandHeight/2 yc-bandHeight/2 ...
                       yc+bandHeight/2 yc+bandHeight/2], ...
                  colDual, 'FaceAlpha', dualAlpha, 'EdgeColor','none', ...
                  'HandleVisibility','off');
        end
    end

    % seconda passata: scatter L1/L5 sopra la banda
    for k = 1:N
        tV = tValidCell{k};
        if isempty(tV); continue; end

        switch band(k)
            case "L1"
                yc  = yPos(k) + offL1;
                col = colL1;
            case "L5"
                yc  = yPos(k) + offL5;
                col = colL5;
            otherwise
                yc  = yPos(k);
                col = [0.5 0.5 0.5];
        end

        scatter(ax, tV, repmat(yc, numel(tV), 1), ...
                mSize, col, 'filled', 'Marker','s');
    end

    % ---- etichette asse Y: "GPS 10", "GAL 03", ... in LaTeX ----
    shortName = strings(height(uniqueSats),1);
    for j = 1:height(uniqueSats)
        switch uniqueSats.cons(j)
            case "GPS",     tag = "GPS";
            case "Galileo", tag = "GAL";
            case "GLONASS", tag = "GLO";
            case "BeiDou",  tag = "BDS";
            otherwise,      tag = uniqueSats.cons(j);
        end
        shortName(j) = sprintf('$\\mathrm{%s}\\ %02d$', tag, uniqueSats.PRN(j));
    end
    yticks(ax, 1:height(uniqueSats));
    yticklabels(ax, shortName);
    ylim(ax, [0.3, height(uniqueSats)+0.7]);
    ylabel(ax, 'Satellite (costellation + PRN)', 'Interpreter','latex','FontSize',13);

    xlabel(ax, 'UTC time', 'Interpreter','latex','FontSize',13);
    ax.XAxis.TickLabelFormat = 'HH:mm:ss';

    % separatori tra costellazioni
    [~, firstIdx] = unique(uniqueSats.cons, 'stable');
    for j = 2:numel(firstIdx)
        yline(ax, firstIdx(j)-0.5, 'k:', 'LineWidth', 1);
    end

    % titolo in LaTeX (escape degli underscore)
    opNameTex = strrep(opName, '_', '\_');
    title(ax, sprintf('Satellite visibility - L1 vs L5 (%s)', opNameTex), ...
          'Interpreter','latex','FontSize',13);

    % legenda
    hL1 = scatter(ax, NaT('TimeZone','UTC'), NaN, mSize, colL1, 'filled','Marker','s');
    hL5 = scatter(ax, NaT('TimeZone','UTC'), NaN, mSize, colL5, 'filled','Marker','s');
    hDF = patch(ax, NaN, NaN, colDual, 'FaceAlpha', dualAlpha, 'EdgeColor','none');
    legend([hL1 hL5 hDF], ...
           {'L1 (G1C / E1B)','L5 (G5I / E5I)','dual frequency'}, ...
           'Location','bestoutside', 'Interpreter','latex','FontSize',13);
end
