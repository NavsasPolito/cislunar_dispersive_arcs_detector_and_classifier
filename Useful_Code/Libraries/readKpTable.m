function kp = readKpTable(txtPath, keepYear)
% READKPTABLE  Build the opts.kpTable input of buildLuGRE_HDF5 from the GFZ
%              "Kp_ap_since_1932.txt" distribution file.
%
%   kp = readKpTable(txtPath)
%   kp = readKpTable(txtPath, keepYear)
%
%   The whole file is parsed first (1932 -> present, ~276k three-hour
%   records); the year selection is applied afterwards, on the assembled
%   table. This keeps the reader independent of how the distribution file is
%   sliced and lets the same call be reused for any mission year.
%
%   INPUT
%     txtPath  - path to Kp_ap_since_1932.txt (GFZ Helmholtz Centre, Niemegk).
%                Format: 30 '#' header lines, then one blank-separated record
%                per three-hour UT interval with columns
%                  YYYY MM DD hh.h hh._m days days_m Kp ap D
%     keepYear - (optional) calendar year(s) to retain, matched on the START
%                of each interval (UT). Default 2025 (the BGM1 mission year).
%                Scalar (2025), vector ([2024 2025]) or 'all'/[] for the
%                whole file.
%
%   OUTPUT
%     kp - table, one row per three-hour Kp interval, sorted by time:
%            .time        datetime, START of the interval, naive UT (no
%                         TimeZone), Format 'dd-MM-yyyy HH:mm:ss'
%            
%            .Kp          double, planetary Kp index, 0 - 9 in thirds
%                         (0, 0.333, 0.667, 1, ...); NaN if missing (-1.000)
%            .ap          double, three-hourly ap amplitude; NaN if missing
%            .definitive  logical, true = definitive (D=1), false = preliminary
%
%          Only .time and .Kp are read by buildLuGRE_HDF5/activityFlag; the
%          remaining columns are carried for traceability and for the
%          interval-overlap test described in NOTE 3 below.
%
%   USAGE with buildLuGRE_HDF5
%     settings.opts.kpTable = readKpTable('Kp_ap_since_1932.txt');        % 2025
%     settings.opts.kpTable = readKpTable('Kp_ap_since_1932.txt', 'all'); % 1932 ->
%     buildLuGRE_HDF5(rootDir, settings);
%


DEFAULT_YEAR = 2025;               % BGM1 mission year

% ===================== INPUT CHECKS =====================================
if nargin < 1 || isempty(txtPath)
    [f,p] = uigetfile({'*.txt','GFZ Kp/ap file (*.txt)'}, 'Select Kp_ap_since_1932.txt');
    if isequal(f,0), error('[readKpTable] No file selected.'); end
    txtPath = fullfile(p,f);
end
if ~isfile(txtPath)
    error('[readKpTable] File not found: %s', txtPath);
end
if nargin < 2, keepYear = DEFAULT_YEAR; end

% ===================== READ THE WHOLE FILE ==============================
% CommentStyle '#' is used instead of HeaderLines so the reader survives a
% change in the number of header lines across GFZ file revisions.
fid = fopen(txtPath,'r');
if fid < 0, error('[readKpTable] Cannot open: %s', txtPath); end
try
    C = textscan(fid, repmat('%f',1,10), ...
                 'CommentStyle','#', 'MultipleDelimsAsOne',true, ...
                 'CollectOutput',true);
catch ME
    fclose(fid);  rethrow(ME);
end
fclose(fid);

M = C{1};
if isempty(M) || size(M,2) ~= 10
    error(['[readKpTable] Unexpected layout in %s: got %d column(s), ', ...
           'expected 10 (YYYY MM DD hh.h hh._m days days_m Kp ap D).'], ...
          txtPath, size(M,2));
end

% ===================== BUILD THE FULL TABLE =============================
t  = datetime(M(:,1), M(:,2), M(:,3)) + hours(M(:,4));   % naive UT, bin start
t.Format = 'dd-MM-yyyy HH:mm:ss';

Kp = M(:,8);   Kp(Kp < 0) = NaN;      % missing flagged -1.000 in the file
ap = M(:,9);   ap(ap < 0) = NaN;      % missing flagged -1
D  = M(:,10) == 1;

kp = table(t , Kp, ap, D, 'VariableNames', ...
           {'time','Kp','ap','definitive'});
kp = sortrows(kp,'time');
kp = kp(~isnat(kp.time),:);

nAll = height(kp);

% ===================== KEEP THE SELECTED YEAR(S) ========================
isTxt   = ischar(keepYear) || isstring(keepYear);
keepAll = isempty(keepYear) || (isTxt && any(strcmpi(string(keepYear),"all")));

if ~keepAll
    if isTxt
        keepYear = str2double(string(keepYear));
    end
    keepYear = keepYear(:)';
    if any(~isfinite(keepYear)) || any(keepYear ~= fix(keepYear))
        error('[readKpTable] keepYear must be a calendar year, a vector of years, or ''all''.');
    end
    kp = kp(ismember(year(kp.time), keepYear), :);      % matched on interval START
    if isempty(kp)
        error('[readKpTable] No Kp record found for year(s) %s in %s.', ...
              strjoin(string(keepYear),', '), txtPath);
    end
    
end


end
