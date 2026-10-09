for jj = 1:20
clearvars -except jj settings
close all
clc
dbstop if error

%% --- Selecting dataset
settings = initSettings_iono(jj);

load(settings.obsPath)      % Loading the dataset

%% --- Cycle slip check
obsSolutions = computeCMCphaseSlip(obsSolutions);

%% --- Window selection

% Find window limits
[idx_begin, idx_end, timeUTC] = findTimeWindow(obsSolutions, settings);

%% --- Arcs mask for restoration
obsSolutions = maskObs(obsSolutions, settings);

%% --- C/N0 analysis
plotCN0(obsSolutions,[idx_begin,idx_end],timeUTC,settings)

%% --- Computing and plotting code STEC estimates
[~,obsSolutions] = PlotTEC_code([],obsSolutions,[idx_begin,idx_end],timeUTC,settings);
[~,obsSolutions] = PlotTEC_phase([],obsSolutions,[idx_begin,idx_end],timeUTC,settings);

%% --- ROT 
obsSolutions = computeROT(obsSolutions,[idx_begin,idx_end],timeUTC,settings);

%% --- ROTI
obsSolutions = computeROTI(obsSolutions,[idx_begin,idx_end],timeUTC,settings);

%% --- STEC leveling
obsSolutions = STECleveling(obsSolutions,[idx_begin,idx_end],timeUTC,settings);

%% --- STEC calibration

%% --- Scintillation proxy index - left for future developments
% obsSolutions = computeS4c(obsSolutions,[idx_begin,idx_end],timeUTC,settings);

%% --- Grazing altitude validation
obsSolutions = grazingAltitude(obsSolutions, timeUTC, settings);

%% --- Event identification
event_catalogue = eventIDv2_1(obsSolutions, settings);

%% --- Satellite Visibility
plotSatelliteVisibility(obsSolutions,settings)

%% --- SAVE results
% --- Command Window Log
fprintf('\n-> Saving data...\n')
save(fullfile(settings.outputPath,'catalogue'),"event_catalogue")
disp('Saved: catalogue')
save(fullfile(settings.outputPath,'obsSolutions'),"obsSolutions")
disp('Saved: obsSolutions')
save(fullfile(settings.outputPath,'settings'),"settings")
disp('Saved: settings')
saveFigures(settings.outputPath)

end

%% --- ALL OP catalogue and statistics
completeCatalogueElaboration(jj,settings)
