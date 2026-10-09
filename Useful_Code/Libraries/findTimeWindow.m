function [idx_window_begin, idx_window_end, OP_time_dat] = findTimeWindow(obsSolutions, settings)
% FINDTIMEWINDOW
% Finds start and end indices of a UTC time window
% from GNSS observation structure.
%
% INPUT:
%   obsSolutions  -> struct containing GPStime and weekNumber
%   settings      -> struct with:
%       epochDate_start  [day, month]
%       epochDate_end    [day, month]
%       tStart_UTC       [hour, minute]
%       tEnd_UTC         [hour, minute]
%
% OUTPUT:
%   idx_window_begin
%   idx_window_end
%   OP_time_dat     -> UTC datetime vector
%
% EXAMPLE USAGE:
%       [idx_begin, idx_end, timeUTC] = findTimeWindow(obsSolutions, settings);
%
%   AUTHOR:
%       Fiorina Francesco, 
%       NavSAS - DET - Politecnico di Torino, 2/03/2026
%
%   VERSION:
%       1.0
%
%   COPYRIGHT:
%       NavSAS - DET - Politecnico di Torino, 2026. All rights reserved.

%% Extract settings
epochDate_start = settings.epochDate_start;
epochDate_end   = settings.epochDate_end;
tStart_UTC      = settings.tStart_UTC;
tEnd_UTC        = settings.tEnd_UTC;

%% Extract GPS Time
tow = obsSolutions(1).GPStime';      
wn  = obsSolutions(1).weekNumber';   

%% Convert GPS → UTC
OP_time_mat = gps_2_Utc([wn, tow]);  
OP_time_dat = datetime(OP_time_mat);

%% Find window limits

idx_window_begin = find( ...
    OP_time_mat(:,2)==epochDate_start(2) & ...
    OP_time_mat(:,3)==epochDate_start(1) & ...
    OP_time_mat(:,4)==tStart_UTC(1) & ...
    OP_time_mat(:,5)==tStart_UTC(2), ...
    1,'first');

idx_window_end = find( ...
    OP_time_mat(:,2)==epochDate_end(2) & ...
    OP_time_mat(:,3)==epochDate_end(1) & ...
    OP_time_mat(:,4)==tEnd_UTC(1) & ...
    OP_time_mat(:,5)==tEnd_UTC(2), ...
    1,'last');


%% Check consistency

if isempty(idx_window_begin)
    error('Start time not found in dataset.')
end

if isempty(idx_window_end)
    error('End time not found in dataset.')
end

if idx_window_begin >= idx_window_end
    error('Invalid time window: begin >= end.')
end

end
