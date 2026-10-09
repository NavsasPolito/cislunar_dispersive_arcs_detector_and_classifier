function buildLuGRE_HDF5(rootDir, settings)
% BUILDLUGRE_HDF5  Assemble the LuGRE dispersive occultation-arc dataset (HDF5).
%
%   buildLuGRE_HDF5(rootDir, outFile, opts)
%
%   Builds a single self-describing HDF5 file from the aggregated catalogue
%   (catalogue_all.mat / catalogue_All.mat) plus, for every operation, its
%   per-OP obsSolutions.mat (time series) and spc_GroundTruth.mat (receiver
%   ECEF). One HDF5 group is written per catalogued ARC (formerly "event";
%   the two are unified as "arc" here).
%
%   Expected on-disk layout (relative to rootDir):
%     Results/Cataloghi_LuGRE_NavSAS_INGV/catalogue_all.mat        (aggregated table)
%     Results/Cataloghi_LuGRE_NavSAS_INGV/CATALOGUES/<OP>/obsSolutions.mat
%     Datasets/<OP>/RTP/SupportData/TRUTH/spc_GroundTruth.mat
%     Datasets/Mission_Scheduler.xml
%
%   INPUT
%     rootDir  - path to Ionosphere_characterization_LuGRE (default: pwd)
%     settings
%       .resultPath - output .h5 path (default: <rootDir>/LuGRE_dispersive_arcs.h5)
%       .opts       - struct of overrides (all optional):
%           .cataloguePath   explicit path to the aggregated catalogue .mat
%           .missionYear     year used for DOY->date conversion   (default 2025)
%           .missionStart    'dd-MM-yyyy' override (skip XML parsing)
%           .missionEnd      'dd-MM-yyyy' override
%           .DOI             DOI string                            (default placeholder)
%           .localTimeMode   'solar' (default) | 'civil'          (see NOTE 3)
%           .grazUnits       units label for grazing altitude      (default 'm')
%           .kpTable         table/struct with .time(datetime) & .Kp(double) for
%                         the activity flag; [] -> 'UNKNOWN'    (default [])
%           .idCounterScope  'perOP' (default) | 'global'
%
%   NOTES / OPEN POINTS (flagged for the author)
%     1. Mission_Scheduler.xml is parsed heuristically (DOY + year). The
%        detected window is printed; verify it, or set opts.missionStart/End.
%     2. Activity flag (Kp / NOAA G-scale) has no source column in the
%        catalogue. It is 'UNKNOWN' unless opts.kpTable is supplied.
%     3. "Local time": default is solar local time (LT = UTC + lon/15 h), the
%        physically meaningful quantity for ionospheric pierce points and
%        dependency-free. 'civil' mode needs a lat/lon->IANA-zone helper that
%        MATLAB does not provide natively (see localTimeSeries).
%     4. Grazing-altitude units are labelled 'm' (sample values ~7.8e6 m);
%        confirm before the Scientific Data deposition.
%     5. Ev.id uses the compact form Lu###@&&&$$$$ as specified here, NOT the
%        earlier hyphenated Lu-###-TYP-SV-OP form. Switch in makeEvID if needed.
%
% -------------------------------------------------------------------------
%   AUTHOR:
%       Francesco Fiorina, NavSAS - DET - Politecnico di Torino
%
%   VERSION HISTORY:
%       1.0 - 20/07/2026 - Initial release: HDF5 arc-catalogue builder with
%             file header (mission metadata, classification accuracy,
%             per-type summary statistics) and per-arc header + body groups.
%
%   COPYRIGHT:
%       NavSAS - DET - Politecnico di Torino, 2026. All rights reserved.
% -------------------------------------------------------------------------

% ===================== DEFAULTS =========================================
if nargin < 1 || isempty(rootDir),              rootDir = uigetdir(pwd, 'Select the Ionosphere_characterization_LuGRE path'); end
if nargin < 2 || isempty(settings.resultPath),  outFile = fullfile(rootDir,'Results','LuGRE_dispersive_arcs12.h5'); else, outFile = fullfile(settings.resultPath,'LuGRE_dispersive_arcs12.h5'); end
if isempty(settings.opts),        settings.opts = struct(); else, opts=settings.opts; end

def = struct('cataloguePath','', 'missionYear',2025, 'missionStart',"15-01-2025", ...
             'missionEnd',"17-03-2025", 'DOI','10.5281/zenodo.XXXXXXXX', ...
             'localTimeMode','solar', 'grazUnits','m', 'kpTable',[], ...
             'idCounterScope','perOP');
fn = fieldnames(def);
for k = 1:numel(fn)
    if ~isfield(opts,fn{k}) || isempty(opts.(fn{k}))
        opts.(fn{k}) = def.(fn{k});
    end
end

DATASET_NAME  = 'Lunar GNSS Receiver Experiment (LuGRE) dispersive occultation arcs';
MISSION_NAME  = 'Blue Ghost Mission 1 (BGM1)';
CONSTELLATION = 'GPS, GAL';
SAMPLING_HZ   = 1;

% Resolve catalogue path
catPath = opts.cataloguePath;
if isempty(catPath)
    cand = { fullfile(rootDir,'Results','catalogue_all.mat'), ...
             fullfile(rootDir,'Results','catalogue_all.mat') };
    catPath = cand{find(cellfun(@isfile,cand),1)};
    if isempty(catPath)
        error('[buildLuGRE] Aggregated catalogue not found. Set opts.cataloguePath.');
    end
end
schedPath = fullfile(rootDir,'Datasets','Mission_Scheduler.xml');

fprintf('============================================================\n');
fprintf('== buildLuGRE_HDF5 - LuGRE dispersive arc catalogue (HDF5) ==\n');
fprintf('============================================================\n');
fprintf('Catalogue : %s\n', catPath);
fprintf('Output    : %s\n\n', outFile);

% ===================== LOAD CATALOGUE ===================================
T = loadCatalogueTable(catPath);
h = height(T);
fprintf('Loaded catalogue: %d arcs, %d columns.\n', h, width(T));

% Pull every column once (missing columns are warned and filled)
OPc   = getStr(T,'OP',h);              SVc   = getStr(T,'SV',h);
phase = getStr(T,'missionPhase',h);    etype = getStr(T,'eventType',h);
egeom = getStr(T,'eventTypeGeom',h);   constl= getStr(T,'constellation',h);
tSAlt = getNum(T,'trackingStartAlt',h);tEAlt = getNum(T,'trackingEndAlt',h);
rotMed= getNum(T,'medianROTtrim',h);   pkROTI= getNum(T,'peakROTI',h);
mnROTI= getNum(T,'minROTI',h);         dStec = getNum(T,'deltaSTEC',h);
durS  = getNum(T,'duration_s',h);
tsT   = getDT(T,'trackingStartTime',h);teT   = getDT(T,'trackingEndTime',h);
evS   = getDT(T,'eventStart',h);       evE   = getDT(T,'eventEnd',h);

% ===================== FILE-HEADER STATISTICS ===========================
uTypes = unique(etype(strlength(etype)>0));

% Classification accuracy = correct arcs (eventType == eventTypeGeom) / total
correct = false(h,1);
for i = 1:h
    correct(i) = strcmpi(strtrim(char(etype(i))), strtrim(char(egeom(i)))) ...
                 && strlength(strtrim(etype(i)))>0;
end
nCorrect = sum(correct);
accPct   = 100 * nCorrect / max(h,1);
fprintf('Classification accuracy: %.2f%% (%d/%d correct).\n', accPct, nCorrect, h);

% Mission window from scheduler (heuristic, override-able)
[missionStart, missionEnd, schedInfo] = readMissionWindow( ...
    schedPath, opts.missionYear, opts.missionStart, opts.missionEnd);
fprintf('Mission window: %s -> %s   (%s)\n\n', missionStart, missionEnd, schedInfo);

% ===================== CREATE FILE + HEADER =============================
if isfile(outFile), delete(outFile); end

% h5num(outFile,'/header/dataset_name', DATASET_NAME);
% h5num(outFile,'/header/doi', settings.opts.DOI);
% h5num(outFile,'/header/mission_start_time', missionStart);
% h5num(outFile,'/header/mission_end_time', missionEnd);
% h5num(outFile,'/header/mission_name', MISSION_NAME);
% h5num(outFile,'/header/sampling_rate_Hz', SAMPLING_HZ);
% h5num(outFile,'/header/active_constellations', CONSTELLATION);
% h5num(outFile,'/header/n_arcs_total',   h);        % creates the file
% h5num(outFile,'/header/n_arcs_correct', nCorrect);
% h5num(outFile,'/header/classification_accuracy_pct', accPct);

% Per-type summary statistics under /header/statistics/<type>
for i = 1:numel(uTypes)
    ty   = uTypes(i);
    mask = etype == ty;
    n    = sum(mask);
    grp  = ['/header/statistics/' sanitize(char(ty))];
    h5num(outFile,[grp '/count'],              n);
    h5num(outFile,[grp '/percent'],            100*n/max(h,1));
    h5num(outFile,[grp '/mean_duration_s'],    mean(durS(mask),'omitnan'));
    h5num(outFile,[grp '/grazing_start_median'], median(tSAlt(mask),'omitnan'));
    h5num(outFile,[grp '/grazing_end_median'],   median(tEAlt(mask),'omitnan'));
    h5num(outFile,[grp '/deltaSTEC_median'],     median(dStec(mask),'omitnan'));
    h5writeatt(outFile,grp,'event_type',char(ty));
    h5writeatt(outFile,[grp '/grazing_start_median'],'units',opts.grazUnits);
    h5writeatt(outFile,[grp '/grazing_end_median'],  'units',opts.grazUnits);
    h5writeatt(outFile,[grp '/deltaSTEC_median'],    'units','TECU');
end

% Root (file header) attributes
h5writeatt(outFile,'/','dataset_name',           DATASET_NAME);
h5writeatt(outFile,'/','DOI',                     char(opts.DOI));
h5writeatt(outFile,'/','mission_name',            MISSION_NAME);
h5writeatt(outFile,'/','mission_start',           char(missionStart));
h5writeatt(outFile,'/','mission_end',             char(missionEnd));
h5writeatt(outFile,'/','sampling_rate_Hz',        SAMPLING_HZ);
h5writeatt(outFile,'/','active_constellations',   CONSTELLATION);
h5writeatt(outFile,'/','classification_accuracy_pct', accPct);
h5writeatt(outFile,'/','n_arcs_total',            h);
h5writeatt(outFile,'/','n_arcs_correct',          nCorrect);
h5writeatt(outFile,'/','created',                 char(datetime('now','Format','yyyy-MM-dd''T''HH:mm:ss')));

% ===================== PER-ARC LOOP (grouped by OP) =====================
ops = unique(OPc,'stable');
globalCount = 0;
nWritten    = 0;

for iop = 1:numel(ops)
    op    = ops(iop);
    rows  = find(OPc == op);            % catalogue rows for this OP (order preserved)
    opStr = char(op);

    % ---- Load this OP's obsSolutions + ground truth once ----
    obsPath   = fullfile(rootDir,'Results',opStr,'obsSolutions.mat');
    truthPath = fullfile(rootDir,'Datasets',opStr,'RTP','SupportData','TRUTH','spc_GroundTruth.mat');

    obs = loadObs(obsPath);
    GT  = loadGroundTruth(truthPath);

    if isempty(obs)
        warning('[buildLuGRE] %s: obsSolutions missing/empty - body series skipped for %d arcs.', ...
                opStr, numel(rows));
        utcNaive = datetime.empty(0,1); svAll = strings(1,0);
    else
        utcNaive = buildOpUTC(obs);               % unzoned UTC (wall time)
        svAll    = string([obs.SV]);              % 1 x Nchan
    end

    fprintf('[%2d/%2d] %-8s : %3d arc(s)', iop, numel(ops), opStr, numel(rows));

    perOPCount = 0;
    for r = rows(:)'
        perOPCount  = perOPCount + 1;
        globalCount = globalCount + 1;
        if strcmpi(opts.idCounterScope,'global'), cnt = globalCount; else, cnt = perOPCount; end

        evID = makeEvID(cnt, etype(r), SVc(r), op);
        base = ['/arcs/' evID];

        % ---- BODY: locate this SV in obsSolutions and window the series ----
        N = 0; idxWin = [];
        haveBody = false;
        if ~isempty(obs)
            id = find(svAll == string(SVc(r)), 1, 'first');
            if ~isempty(id)
                t0 = naiveUTC(tsT(r));  t1 = naiveUTC(teT(r));
                if ~isnat(t0) && ~isnat(t1)
                    win    = (utcNaive >= t0) & (utcNaive <= t1);
                    idxWin = find(win);
                    N      = numel(idxWin);
                    haveBody = N > 0;
                end
            else
                warning('[buildLuGRE] %s: SV %s not found in obsSolutions.', opStr, SVc(r));
            end
        end

        % n_samples marker guarantees the arc group exists for attributes
        h5num(outFile,[base '/body/n_samples'], N);

        grazMin = NaN;
        if haveBody
            gg = obs(id).grazingAltitude(:);
            grazMin = min(gg(idxWin),[],'omitnan');

            % 1-D series
            h5vec(outFile,[base '/body/roti'],           obs(id).ROTI(idxWin),        'TECU/s');
            h5vec(outFile,[base '/body/rot'],            obs(id).ROT(idxWin),         'TECU/s');
            h5vec(outFile,[base '/body/stec_levelled'],  obs(id).smoothedTEC(idxWin), 'TECU');
            h5vec(outFile,[base '/body/stec_phase'],     obs(id).STEC_phase(idxWin),  'TECU');
            h5vec(outFile,[base '/body/stec_code'],      obs(id).STEC_code(idxWin),   'TECU');
            h5vec(outFile,[base '/body/grazing_altitude'],obs(id).grazingAltitude(idxWin), opts.grazUnits);

            % Tangent-point -> geodetic latitude / longitude
            tp = obs(id).tangentPoint(idxWin,:);
            [lat,lon] = ecef2latlon(tp);
            h5vec(outFile,[base '/body/latitude'],  lat, 'deg');
            h5vec(outFile,[base '/body/longitude'], lon, 'deg');

            % ECEF series (SV from obs, receiver from ground truth)
            h5mat(outFile,[base '/body/sv_ecef'], obs(id).SVpos(idxWin,:), 'm');
            rx = groundTruthWindow(GT, idxWin);
            h5mat(outFile,[base '/body/rx_ecef'], rx, 'm');

            % Time stamps
            uw   = utcNaive(idxWin);
            uwz  = uw; 
            uwz.TimeZone = 'UTC';
            
            h5vecstr(outFile,[base '/body/utc'], string(uwz,'dd-MM-yyyy HH:mm:ss'));
            h5writeatt(outFile,[base '/body/utc'],'units','UTC dd-MM-yyyy HH:mm:ss');
            
            lst = localTimeSeries(uw, lon, opts.localTimeMode);
            h5vecstr(outFile,[base '/body/local_time'], string(lst,'dd-MM-yyyy HH:mm:ss'));
            h5writeatt(outFile,[base '/body/local_time'],'units','LT dd-MM-yyyy HH:mm:ss');
        end
        
        % ---- HEADER attributes on the arc group ----
        h5writeatt(outFile,base, 'op_id',         opStr);
        
        h5writeatt(outFile,base,'sv_id',         char(SVc(r)));
        h5writeatt(outFile,base,'mission_phase', char(phase(r)));


        Kp_flag= activityFlag(naiveUTC(tsT(r)), naiveUTC(teT(r)), opts.kpTable);
        h5writeatt(outFile,base,'activity_flag', char(Kp_flag(2)));
        h5writeatt(outFile,base,'Kp', char(Kp_flag(1)));
        %h5writeatt(outFile,[base,'ev_id'],         evID);
        %h5writeatt(outFile,[base,'ev_id'],         evID);

        % Arc block
        h5writeatt(outFile,base,'arc_grazing_start', tSAlt(r));
        h5writeatt(outFile,base,'arc_grazing_min',   grazMin);
        h5writeatt(outFile,base,'arc_grazing_end',   tEAlt(r));
        h5writeatt(outFile,base,'arc_grazing_units', opts.grazUnits);
        h5writeatt(outFile,base,'arc_utc_start',     isoStr(tsT(r)));
        h5writeatt(outFile,base,'arc_utc_end',       isoStr(teT(r)));

        % RO-arc block
        h5writeatt(outFile,base,'ro_rot_median', rotMed(r));
        h5writeatt(outFile,base,'ro_max_roti',   pkROTI(r));
        h5writeatt(outFile,base,'ro_min_roti',   mnROTI(r));
        h5writeatt(outFile,base,'ro_utc_start',  isoStr(evS(r)));
        h5writeatt(outFile,base,'ro_utc_end',    isoStr(evE(r)));
        h5writeatt(outFile,base,'ro_delta_stec', dStec(r));
        h5writeatt(outFile,base,'ro_type',       char(etype(r)));
        h5writeatt(outFile,base,'ro_status',     statusFlag(etype(r), egeom(r)));
        h5writeatt(outFile,base,'constellation', char(constl(r)));

        nWritten = nWritten + 1;
        
    end
    fprintf('  -> done.\n');
end

fprintf('\nWrote %d arc groups to %s\n', nWritten, outFile);
fprintf('============================================================\n');
end % ===================== END MAIN ====================================


% =========================================================================
%                              HELPERS
% =========================================================================

function T = loadCatalogueTable(p)
S  = load(p);
fn = fieldnames(S);
T  = [];
for k = 1:numel(fn)
    if istable(S.(fn{k})), T = S.(fn{k}); return; end
end
error('[buildLuGRE] No table variable found in %s', p);
end

function obs = loadObs(p)
obs = [];
if ~isfile(p), return; end
S  = load(p);
if isfield(S,'obsSolutions'), obs = S.obsSolutions; return; end
fn = fieldnames(S);                              % fallback: first struct array
for k = 1:numel(fn)
    if isstruct(S.(fn{k})) && numel(S.(fn{k}))>1, obs = S.(fn{k}); return; end
end
end

function GT = loadGroundTruth(p)
GT = [];
if ~isfile(p), return; end
S  = load(p);
if isfield(S,'GroundTruth'), GT = S.GroundTruth; return; end
fn = fieldnames(S);
for k = 1:numel(fn)
    v = S.(fn{k});
    if isstruct(v) && isfield(v,'pos_ecef'), GT = v; return; end
end
end

function M = gps_2_Utc(gps_data, leap_sec)
% gps_2_Utc Converts GPS week number and time of week to UTC [Y M D H Mi S]
%
% INPUTS:
%   gps_data: Nx2 matrix where column 1 is Week Number (wn) 
%             and column 2 is Time of Week in seconds (tow).
%   leap_sec: (Optional) Number of leap seconds GPS is ahead of UTC. 
%             Defaults to 18 (valid for dates after Dec 31, 2016).
%
% OUTPUT:
%   M: Nx6 matrix formatted as [Year, Month, Day, Hour, Minute, Second]

    % Set default leap seconds if not provided
    if nargin < 2
        leap_sec = 18; 
    end
    
    % Extract week number and time of week
    wn = gps_data(:, 1);
    tow = gps_data(:, 2);
    
    % Define the GPS starting epoch (Jan 6, 1980 00:00:00)
    gps_epoch = datetime(1980, 1, 6, 0, 0, 0);
    
    % Calculate total elapsed seconds (1 week = 604,800 seconds)
    % Subtract leap seconds to shift from GPS time to UTC time
    total_seconds = (wn .* 604800) + tow - leap_sec;
    
    % Generate the UTC datetime object array
    utc_time = gps_epoch + seconds(total_seconds);
    
    % Format the output into the requested Nx6 matrix
    M = [year(utc_time), month(utc_time), day(utc_time), ...
         hour(utc_time), minute(utc_time), second(utc_time)];
end

function utcNaive = buildOpUTC(obs)
% Build the operation UTC vector (unzoned wall time) from GPS time.
tow = double(obs(1).GPStime(:));
wn  = double(obs(1).weekNumber(:));
M   = gps_2_Utc([wn, tow]);          % [Y M D H Mi S]
utcNaive = datetime(M);              % unzoned; values are UTC
if ~isempty(utcNaive.TimeZone), utcNaive.TimeZone = ''; end
end

function d = naiveUTC(d)
% Normalise a catalogue datetime to unzoned UTC wall time for comparison.
if isdatetime(d)
    if ~isempty(d.TimeZone), d.TimeZone = ''; end
elseif isnumeric(d)
    d = datetime(d,'ConvertFrom','datenum');
end
end

function rx = groundTruthWindow(GT, idxWin)
rx = nan(numel(idxWin),3);
if isempty(GT) || ~isfield(GT,'pos_ecef'), return; end
P = GT.pos_ecef;
if size(P,1) >= max(idxWin)
    rx = P(idxWin,:);                % GT aligned 1:1 with the obs epoch grid
else
    warning('[buildLuGRE] Ground truth shorter than obs grid (%d < %d).', size(P,1), max(idxWin));
end
end

function [lat,lon] = ecef2latlon(tp)
% Geodetic lat/lon of ECEF points via ecef2lla (NaN-safe).
n   = size(tp,1);
lat = nan(n,1);  lon = nan(n,1);
good = all(isfinite(tp),2);
if any(good)
    lla = ecef2lla(tp(good,:));      % [lat lon alt] deg,deg,m
    lat(good) = lla(:,1);
    lon(good) = lla(:,2);
end
end

function lst = localTimeSeries(utcNaiveWin, lon, mode)
% Local time series.
%   'solar' (default): mean solar local time = UTC + lon/15  [hours, 0-24)
%   'civil'          : requires a lat/lon->IANA-zone helper (not native to
%                      MATLAB). Provide latlon2tz.m to enable; falls back to
%                      solar with a warning otherwise.
switch lower(mode)
    case 'civil'
        % latlon2tz.m  - not yet implemented future developments
        % if exist('latlon2tz','file') 
        %     % Per-sample civil local time (expensive; needs external tz data).
        %     lst = nan(numel(utcNaiveWin),1);
        %     for k = 1:numel(utcNaiveWin)
        %         try
        %             z = latlon2tz(NaN, lon(k));            
        %             t = utcNaiveWin(k); 
        %             t.TimeZone = 'UTC'; 
        %             t.TimeZone = z;
        %             lst(k) = hour(t) + minute(t)/60 + second(t)/3600;
        % 
        %         catch
        %             lst(k) = NaN;
        %         end
        %     end
        %     return
        % else
        %     warning('[buildLuGRE] localTimeMode=''civil'' needs latlon2tz.m; using solar LT.');
        % end
        error('[buildLuGRE] localTimeMode=''civil'' needs latlon2tz.m helper which must built. Until this helper does not exist use -solar- mode for local time')
end

lst                  = utcNaiveWin + (hours(lon/15));   % local
end

function L = typeLetter(et)
et = lower(strtrim(char(et)));
if     contains(et,'ingress'),  L = 'I';
elseif contains(et,'egress'),   L = 'E';
elseif contains(et,'deep'),     L = 'D';   % deep tangent (extension)
elseif contains(et,'tangent'),  L = 'T';
elseif contains(et,'unclass'),  L = 'U';   % unclassified (extension)
else,                            L = 'X';
end
end

function id = makeEvID(cnt, et, sv, op)
% Compact form:  Lu###@&&&$$$$   (### counter, @ type letter, &&& SV, $$$$ OP)
% For the earlier hyphenated form, replace the line below with:
%   id = sprintf('Lu-%03d-%s-%s-%s', cnt, typeCode3(et), char(sv), char(op));
id = sprintf('Lu%03d%s%s%s', cnt, typeLetter(et), char(sv), char(op));
end

function s = statusFlag(et, eg)
if strcmpi(strtrim(char(et)), strtrim(char(eg))) && strlength(strtrim(string(et)))>0
    s = 'Correct';
else
    s = 'Wrong';
end
end

function Kp = activityFlag(t0, t1, kp)
% NOAA G-scale from max Kp in the window, if a Kp table is supplied.
%flag = "UNKNOWN";
if isempty(kp) || isnat(t0) || isnat(t1), return; end
try
    
    m = (kp.time >= t0-hours(3)) & (kp.time <= t1+hours(3));
    
    if ~any(m) 
        return; 
    end
    kmax = max(kp.Kp(m));
    if     kmax >= 9, flag = "G5 (extreme)";
    elseif kmax >= 8, flag = "G4 (severe)";
    elseif kmax >= 7, flag = "G3 (strong)";
    elseif kmax >= 6, flag = "G2 (moderate)";
    elseif kmax >= 5, flag = "G1 (minor)";
    else,             flag = "quiet";
    end
catch
    flag = "UNKNOWN";
    kmax='UNKNOWN';
end
Kp=[kmax, flag];
disp(Kp)
end

function [sStr, eStr, info] = readMissionWindow(xmlPath, yearHint, sOv, eOv)
% Heuristic parse of Mission_Scheduler.xml: min/max day-of-year -> dd-MM-yyyy.
sStr = string(sOv);  eStr = string(eOv);  info = "";
if strlength(sStr)>0 && strlength(eStr)>0, info = "override"; return; end
if ~isfile(xmlPath)
    warning('[buildLuGRE] Mission_Scheduler.xml not found: %s', xmlPath);
    if strlength(sStr)==0, sStr = "UNKNOWN"; end
    if strlength(eStr)==0, eStr = "UNKNOWN"; end
    return
end
txt = fileread(xmlPath);

% Year
yk  = regexp(txt,'(?<!\d)(20\d\d)(?!\d)','tokens');
yrs = str2double(cellfun(@(c)c{1}, yk, 'uni', false));
yrs = yrs(yrs>=2024 & yrs<=2027);
if ~isempty(yrs), yr = mode(yrs); else, yr = yearHint; end

% Day-of-year candidates near likely keywords
doy = [];
kw  = {'DOY','doy','DayOfYear','dayOfYear','day_of_year','day-of-year'};
for k = 1:numel(kw)
    tk = regexp(txt,[kw{k} '\D{0,6}(\d{1,3})'],'tokens');
    for j = 1:numel(tk), doy(end+1) = str2double(tk{j}{1}); end %#ok<AGROW>
end
doy = doy(doy>=1 & doy<=366);

if isempty(doy)
    warning(['[buildLuGRE] Could not auto-detect DOY in scheduler. ', ...
             'Set opts.missionStart / opts.missionEnd (dd-MM-yyyy).']);
    if strlength(sStr)==0, sStr = "UNKNOWN"; end
    if strlength(eStr)==0, eStr = "UNKNOWN"; end
    return
end

d0 = datetime(yr,1,1) + caldays(min(doy)-1);  d0.Format = 'dd-MM-yyyy';
d1 = datetime(yr,1,1) + caldays(max(doy)-1);  d1.Format = 'dd-MM-yyyy';
if strlength(sStr)==0, sStr = string(d0); end
if strlength(eStr)==0, eStr = string(d1); end
info = sprintf('heuristic: year=%d, DOY=[%d,%d]', yr, min(doy), max(doy));
end

% ---- small typed accessors (missing columns warned + filled) ------------
function v = getStr(T,name,h)
if ismember(name,T.Properties.VariableNames), v = string(T.(name));
else, warning('[buildLuGRE] string column "%s" missing.',name); v = strings(h,1); end
end
function v = getNum(T,name,h)
if ismember(name,T.Properties.VariableNames), v = double(T.(name));
else, warning('[buildLuGRE] numeric column "%s" missing.',name); v = nan(h,1); end
end
function v = getDT(T,name,h)
if ismember(name,T.Properties.VariableNames)
    v = T.(name);
    if ~isdatetime(v), v = datetime(v); end %    if ~isdatetime(v), v = datetime(v,'ConvertFrom','datenum'); end
    if ~isempty(v.TimeZone), v.TimeZone = ''; end
else
    warning('[buildLuGRE] datetime column "%s" missing.',name); v = NaT(h,1);
end
end

% ---- HDF5 write helpers -------------------------------------------------
function h5num(f,name,val)
if isempty(val) || ~isscalar(val), val = NaN; end
h5create(f,name,[1 1],'Datatype','double');
h5write(f,name,double(val));
end
function h5vec(f,name,data,units)
data = double(data(:));
if isempty(data), return; end
h5create(f,name,[numel(data) 1],'Datatype','double');
h5write(f,name,data);
if nargin>=4 && ~isempty(units), h5writeatt(f,name,'units',units); end
end
function h5vecstr(f,name,data)
data = string(data(:));
if isempty(data), return; end
h5create(f,name,[numel(data) 1],'Datatype','string');
h5write(f,name,data);
end
function h5mat(f,name,data,units)
if isempty(data), return; end
h5create(f,name,size(data),'Datatype','double');
h5write(f,name,double(data));
if nargin>=4 && ~isempty(units), h5writeatt(f,name,'units',units); end
end

function s = isoStr(dt)
if ~isdatetime(dt) || isnat(dt), s = 'NaT'; return; end
dt.Format = 'yyyy-MM-dd''T''HH:mm:ss';
s = char(dt);
end

function s = sanitize(str)
s = regexprep(char(str),'[^A-Za-z0-9_]','_');
end

