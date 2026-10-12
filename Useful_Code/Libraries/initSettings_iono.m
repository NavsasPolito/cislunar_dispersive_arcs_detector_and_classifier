function settings = initSettings_iono(jj)
%INITSETTINGS Initializes configuration parameters for the GNSS
% post-processing pipeline.
%
% DESCRIPTION
% -------------------------------------------------------------------------
% This function creates and returns a structure containing all the user
% configuration parameters required by the GNSS post-processing toolbox.
%
% -------------------------------------------------------------------------
% OUTPUT
%
% settings : struct
%     Structure containing all configuration parameters used throughout
%     the processing pipeline.
%
% -------------------------------------------------------------------------
%   AUTHOR:
%       Francesco Fiorina, NavSAS - DET - Politecnico di Torino
%
%   VERSION HISTORY:
%       1.0 - 10/03/2026 - Initial release: dataset path, time window,
%             constellation selection, CN0 threshold, arc restoration,
%             leveling gap, and ROTI processing parameters.
%       2.0 - 13/04/2026 - Added eventIDv2 classification parameters:
%             MAX_TANGENT_GAP_S, OCC_GAP_S, PEAK_EDGE_TOL,
%             MEAN_THRESHOLD, MIN_DURATION.
%       3.0 - 20/05/2026 - Added thresholds for side lobe exclusion (or
%             mainlobe inclusion.
%       4.0 - 18/06/2026 - Added plot switch for S4c, S4c and sigmaphi
%             parameters.
%
%   COPYRIGHT:
%       NavSAS - DET - Politecnico di Torino, 2026. All rights reserved.
% -------------------------------------------------------------------------

%% --- Command Window log
fprintf('============================================================\n')
fprintf('==      GNSS observables ionospheric processing tool      ==\n')
fprintf('==                      NavSAS Group                      ==\n')
fprintf('==              written by Francesco Fiorina              ==\n')
fprintf('============================================================\n\n')


%% --- DATASET PATH
% Input Path to the folder containing the dataset
% (can alternatively be selected interactively with uigetdir)
% settings.inputPtah  = uigetdir('Select dataset folder');

OPs = { 'OP1_0',  'OP2_0',  'OP5_0',  'OP9_0',  'OP12_0', 'OP14_0', 'OP17_0', ...
        'OP18_0', 'OP21_0', 'OP22_0', 'OP23_0', 'OP27_0', 'OP37_0', 'OP38_0',...
        'OP40_0', 'OP74_0', 'OP76_0', 'OP77_0', 'OP77_1', 'OP78_1'};

settings.OP = string((OPs(jj)));

% --- Root folder: insert the folder containing Useful_Code, Datasets and Results.
settings.rootPath = '/absolute/path/to/root'; 

% --- Input paths
settings.datasetsPath = fullfile(settings.rootPath, 'Datasets');
settings.inputPath    = fullfile(settings.datasetsPath, settings.OP, 'RTP', 'SupportData');

% --- OPTABLE path
settings.opTablePath = fullfile(settings.datasetsPath,'OP TABLE.csv');

% --- kp table path
settings.kpTablePath = fullfile(settings.datasetsPath,'Kp_ap_since_1932.txt'); 

% Observables
settings.obsFilename = 'obsSolutions.mat';
settings.obsPath     = fullfile(settings.inputPath,'OBS', settings.obsFilename);

% Rinex
navDir = fullfile(settings.inputPath,'BROADCAST','NAV_multi');
rnxFile = dir(navDir);
rnxFile = rnxFile(~ismember({rnxFile.name},{'.','..'}));
settings.rnxFile = rnxFile;
settings.rnxPath = fullfile(navDir, rnxFile.name);

% Ground truth
settings.gtFile = 'spc_GroundTruth.mat';
settings.gtPath = fullfile(settings.inputPath,'TRUTH', settings.gtFile);

% Off-Boresight Angles
% settings.obaFile = 'offBoresightAngles.mat';
% settings.obaPath = fullfile(settings.inputPath,'SV', settings.obaFile);

% --- Output path
settings.resultPath = fullfile(settings.rootPath, 'Results');
settings.outputPath = fullfile(settings.resultPath, settings.OP);


if ~isfolder(settings.outputPath)
    [ok, msg] = mkdir(settings.outputPath);
    if ~ok
        error('initSettings:mkdirFailed', ...
            'Impossible creating %s: %s', settings.outputPath, msg);
    end
end

%% --- PROCESSING TIME WINDOW
settings.window = 'auto';           % 'manual' or 'auto' mode for window selection

switch settings.window
    case 'manual'
        settings.epochDate_start = [4 , 3];   % Start date of the analysis interval [day, month]
        settings.epochDate_end   = [4 , 3];   % End date of the analysis interval [day, month]
        settings.tStart_UTC      = [7 , 54];  % Start time (UTC) [hour, minute]
        settings.tEnd_UTC        = [20 , 31]; % End time (UTC) [hour, minute]
    case 'auto'
        tmp = load(settings.obsPath);  obsSolutions=tmp.obsSolutions;

        wn  = obsSolutions(1).weekNumber;
        tow = obsSolutions(1).GPStime;

        timeUTC = gps_2_Utc([wn', tow']);

        settings.epochDate_start = [timeUTC(1,3)   , timeUTC(1,2)  ];   % Start date of the analysis interval [day, month]
        settings.epochDate_end   = [timeUTC(end,3) , timeUTC(end,2)];   % End date of the analysis interval [day, month]
        settings.tStart_UTC      = [timeUTC(1,4)   , timeUTC(1,5)  ];   % Start time (UTC) [hour, minute]
        settings.tEnd_UTC        = [timeUTC(end,4) , timeUTC(end,5)];   % End time (UTC) [hour, minute]
end


%% --- GNSS CONSTELLATION SELECTION
settings.G = 1; % Enable GPS processing
settings.E = 1; % Enable Galileo processing

%% --- MEASUREMENT QUALITY THRESHOLDS
settings.CN0_threshold = 20; % Minimum C/N0 threshold (dB-Hz) for usable observations

%% --- MAINLOBE INCLUSION THRESHOLDS - not yet implemented future developments
% settings.mainlobe_threshold.G.L1_E1 = Inf;25; % [deg] Off-Boresight angle threshold for main lobe inclusion
% settings.mainlobe_threshold.G.L5_E5 = Inf;32; % [deg] Off-Boresight angle threshold for main lobe inclusion
% 
% settings.mainlobe_threshold.E.L1_E1 = Inf;25; % [deg] Off-Boresight angle threshold for main lobe inclusion
% settings.mainlobe_threshold.E.L5_E5 = Inf;27; % [deg] Off-Boresight angle threshold for main lobe inclusion
% 
% settings.grazing_threshold = Inf;         % [m] eg. 7e6 m --> if event starts over this threshold is considered as sidelobe

%% --- ROTI ARCS RESTORATION
settings.arcRestore = true; % Fills the ROTI arcs gaps with interpolated values

%% --- FILTER / PROCESSING PARAMETERS
settings.t_gap = 10; % Maximum allowed observation gap (seconds) before leveling filter reinitialization

%% --- ROTI PROCESSING PARAMETERS
settings.ROTI_window = 'fixed'; % length of the window over which is computed the ROTI 
                                   % - 'variable' : length varies between on each satellite (to be implemented)
                                   % - 'fixed'    : length is fixed at settings.window_length [s] (typical approach)
settings.window_length = 10*3;  % window length in seconds 
settings.tol_window    = 0.1;   % tolerance on missing samples 

%% --- EVENTID PARAMETERS
% Parameters for event classification
settings.restoredROTIclass = false;     % if true the classification elaborates the restored ROTI
                                        %    with a risk of classifying on the basis of an
                                        %    interpolation artifact due to restoration.

% --- Segmentation and grouping
settings.MAX_TANGENT_GAP_S = 5 * 60;    % [s] max gap between segments to still be
                                        %     considered part of the same tangent event
                                        %     (sampling is 1 Hz → gap in samples = gap in s)

settings.OCC_GAP_S         = 30 * 60;   % [s] min gap between two arcs to flag them as
                                        %     an ingress/egress occultation pair

settings.MIN_DURATION      = 100;       % [s] minimum event duration to be classified

% --- Classification thresholds
settings.PEAK_EDGE_TOL     = 0.05;      % [-] fraction of arc length defining the edge
                                        %     tolerance zone for peak detection
                                        %     (5% → ~3 s for a 60-s arc)

settings.MEDIAN_THRESHOLD    = 0.007;  % [TECu/s] threshold on |medianROT|:
                                        %     |medianROT| > threshold → ingress or egress
                                        %     |medianROT| < threshold → tangent

% --- Classification validation threshold
settings.DEEP_TANGENT_HEIGHT = 6e5;     % [m] threshold for defining a deep tangent event,
                                        %      hence sufficiently close to the F2 peak.

%% --- S4c PARAMMETERS - not yet implemented future developments
% settings.det_window             = 60;       % [s] window for detrending SI. 60 s at 1 Hz 
% settings.avg_window             = 60;       % [s] window for averaging detrended SI. 60 s at 1 Hz
%                                             % value in (0, 1] that represents the *maximum allowed
% settings.nan_threshold_S4c      = nan;      % fraction* of NaN values inside an averaging window
%                                             % CURRENTLY NOT SET --> CREATE MOVMEAN_NANCHECK

%% --- SIGMA_PHI PARAMETERS - not yet implemented future developments
% settings.sigma_phi_window       = 60;       % [s] window for averaging detrended phase. 60 s at 1 Hz
% settings.nan_threshold_sigmaPhi = 0;        % value in (0, 1] that represents the *maximum allowed
%                                             % fraction* of NaN values inside an averaginng window                                       

%% --- BENDING ANGLE PARAMETERS - not yet implemented future developments
% settings.d_phi         = 0.00001; % Correction increment in radians
% settings.h             = 1e-6;    % Correction increment
% settings.epsilon_fd    = 1e-3;    % Doppler difference tolerance
% settings.discard       = true;    % Discard bending angle when ray penetrates neutral atmosphere (tropo/strato). Ionospheric assumption invalid.
% settings.max_iter      = 1000;    % Iteration limit on gauss newton on angles iteration

%% --- PLOT SWITCHES
settings.PlotCN0                            = true;
% settings.PlotS4c                            = true; %S - not yet implemented future developments
settings.PlotSTECphase                      = true;
settings.PlotSTECcode                       = true;
settings.PlotSTECleveled                    = true;
settings.PlotROT                            = true;
settings.PlotROTI                           = true;
settings.PlotROTI_restored                  = true;
settings.PlotVisibility                     = true;
settings.PlotGrazingAltitude                = true;
settings.PlotGrazingAltitudeArcs            = true;
settings.PlotGrazingAltitudeDistribution    = true;

%% --- HDF5 file generation
settings.opts.cataloguePath   = fullfile(settings.resultPath,'catalogue_all.mat');
settings.opts.missionStart    = '15-01-2025-06:11:24';
settings.opts.missionEnd      = '17-03-2025-00:38:17';
settings.opts.DOI             = '10.5281/zenodo.23216092';
settings.opts.localTimeMode   = 'solar';  % 'civil'          (see NOTE 3 in buildLuGRE_HDF5.m)
settings.opts.grazUnits       = 'm';
settings.opts.kpTable         = readKpTable(settings.kpTablePath);
settings.opts.idCounterScope  = 'global'; % "Local time": default is solar local time (LT = UTC + lon/15 h), the physically meaningful 
                                          % quantity for ionospheric pierce points and dependency-free. 'civil' mode needs a lat/lon->IANA-zone 
                                          % helper that MATLAB does not provide natively and must be built (see localTimeSeries and NOTE 3 in buildLuGRE_HDF5.m).
                                        
end