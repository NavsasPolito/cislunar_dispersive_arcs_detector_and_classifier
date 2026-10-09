function obsSolutions = grazingAltitude(obsSolutions, timeUTC, settings)
%GRAZINGALTITUDE  Straight-line tangent altitude (SLTA) of every
%                 GNSS-to-receiver raypath over the processing time
%                 window, in WGS-84 geodetic coordinates.
%
% DESCRIPTION
% -------------------------------------------------------------------------
% For each epoch in timeUTC and for each GNSS satellite tracked in
% obsSolutions, this function:
%
%   1) Determines the signal time of flight dt = rawP/c from the
%      pseudorange of the first non-NaN frequency among L1/L5.
%      Epochs where neither frequency is tracked are skipped.
%   2) Propagates the GNSS SV position and velocity from the broadcast
%      ephemeris at the transmission time tTX = tRX - dt using
%      gnssconstellation (closest-ephemeris selection per PRN via
%      selectClosestEphemeris).
%   3) Applies the Sagnac correction, rotating the SV state vector
%      from ECEF(tTX) into ECEF(tRX) by omegaE * dt around the Z axis,
%      so that all geometry is computed in a single ECEF frame.
%   4) Computes the Straight-Line Tangent Altitude (SLTA) of the
%      SV-to-receiver line: the geodetic WGS-84 height of the point
%      on the signal line closest to the Earth's centre.
%   5) Stores the SV state, the grazing altitude, and the ECEF
%      coordinates of the tangent point into all obsSolutions rows
%      associated with that SV.
%
% GEOMETRY
% -------------------------------------------------------------------------
% Let r_rx be the receiver position and r_GNSS be the (Sagnac-corrected)
% GNSS transmitter position in ECEF. The signal line is
%
%       r(s) = r_rx + s * t_hat,     t_hat = (r_GNSS - r_rx)/|r_GNSS - r_rx|
%
% The Earth-centre distance of this line (the straight-line impact
% parameter) is obtained from the cross-product identity
%
%       |r_rx x r_t| = |r_rx| * |r_t| * sin(phi)
%
% with r_t = r_GNSS - r_rx and phi the angle between r_rx and r_t, so
%
%       altitude = |r_rx| * sin(phi) = |r_rx x r_t| / |r_t|
%
% (this scalar, representing the geocentric distance of the line, is
% retained internally for the sanity check below).
%
% The tangent point is the foot of the perpendicular from the origin
% onto the signal line:
%
%       s*            = -(r_rx . t_hat)
%       tangent_point = r_rx + s* * t_hat = r_rx - (r_rx . t_hat) * t_hat
%
% (note the minus sign: the foot-of-perpendicular parameter is
% -r_rx.t_hat, as required by the perpendicularity condition
% tangent_point . t_hat = 0).
%
% The grazing altitude returned is the WGS-84 geodetic height of the
% tangent point, obtained via ecef2geodetic. This is the standard
% definition of SLTA used in operational radio occultation
% (Springer Handbook Fig. 38.12). It accounts for the oblateness of
% the Earth, with up to ~14 km difference at the poles relative to a
% spherical-Earth approximation.
%
% LIGHT-TIME AND BIAS BUDGET
% -------------------------------------------------------------------------
% The time of flight dt = rawP/c is used directly, without solving the
% implicit equation t_tx = t_rx - |r_rx - r_GNSS(t_tx)|/c. The
% pseudorange rawP includes minor biases (satellite clock residual,
% ionospheric and tropospheric delays, code noise, possibly receiver
% clock offset) that map into errors on dt:
%   - satellite clock residual         ~1   m  -> ~3 ns -> ~10 mm SV shift
%   - ionospheric delay (in occ.)   ~1-20  m  -> ~3-70 ns -> ~0.3 m SV shift
%   - tropospheric delay         ~0 m (vacuum link)
%   - code noise                    ~1   m  -> ~3 ns -> ~10 mm SV shift
%   - receiver clock offset        depends on dataset
% With a receiver clock disciplined to better than ~1 ms, the resulting
% tangent-point error is at most a few metres, well below the
% kilometre-scale features of interest for occultation masking. For
% datasets with grossly uncalibrated receiver clocks, rawP should be
% corrected upstream.
%
% NOTE ON THE METHOD AT LUNAR ALTITUDES
% -------------------------------------------------------------------------
% This straight-line approach is preferred over the Doppler-based
% bending-angle inversion (Schreiner et al., 1999; Springer Handbook
% Chap. 38) for occultation masking and tagging at lunar-altitude
% receivers, for the following reasons:
%
%   - the occultation plane r_GNSS x r_rx is well-conditioned only for
%     LEO receivers; at lunar distances r_GNSS and r_rx can become
%     nearly collinear and the Doppler iteration becomes unstable;
%   - the Doppler signature of a given bending angle scales with the
%     receiver velocity projected in the occultation plane: ~7.5 km/s
%     for LEO, ~1 km/s for a lunar receiver, so the SNR on the
%     atmospheric Doppler drops by an order of magnitude;
%   - the ray separation from straight-line propagation is already
%     below the Fresnel zone (~1.5 km) at LEO (Schreiner et al., 1999,
%     Fig. 3 and Table 1); from lunar altitudes it is even smaller
%     relative to the total ray length.
%
% At lunar distances the light time can reach ~1 s, so both the
% transmission-time correction (~1 km of SV motion) and the Sagnac
% rotation (~2-3 km of frame rotation) are non-negligible and are
% explicitly applied here.
%
% INPUT
% -------------------------------------------------------------------------
% obsSolutions : struct array
%     Observation structure with one entry per tracked SV/frequency.
%     Each entry must have:
%         .SV    : char/string of the form 'G##' / 'E##'
%         .rawP  : [N_epochs x 1] pseudorange [m], NaN when not tracked
%
% timeUTC : datetime array [N_epochs x 1]
%     Epoch grid over which the grazing altitude is computed.
%
% settings : struct with fields
%     .rnxPath : path to the RINEX navigation file
%     .gtPath  : path to a .mat file containing a GroundTruth struct
%                with fields .pos_ecef [N_epochs x 3] [m]
%                          and .vel_ecef [N_epochs x 3] [m/s]
%     .G       : logical, enable GPS
%     .E       : logical, enable Galileo
%     .PlotGrazingAltitude     : logical, enable the bare arcs plot
%     .PlotGrazingAltitudeArcs : logical, enable the observable-masked
%                                arcs plot
%
% OUTPUT
% -------------------------------------------------------------------------
% obsSolutions : struct array
%     Updated structure. For each entry the following fields are
%     written (one row per epoch):
%         .SVpos            [N_epochs x 3]  SV ECEF position at tRX [m]
%                                           (Sagnac-corrected)
%         .SVvel            [N_epochs x 3]  SV ECEF velocity at tRX [m/s]
%                                           (Sagnac-corrected)
%         .grazingAltitude  [N_epochs x 1]  SLTA above WGS-84 [m]
%         .tangentPoint     [N_epochs x 3]  tangent-point ECEF [m]
%         .SVtrackingPoints [Mx4]           tracking-window records
%                                           (only when the arcs plot
%                                           is enabled and SV is dual-
%                                           frequency); see
%                                           retrieveTrackingPoints.m
%
% REFERENCES
% -------------------------------------------------------------------------
% - Springer Handbook of Global Navigation Satellite Systems (2017),
%   Chap. 38, Fig. 38.12 ("Straight Line Tangent point Altitude").
% - Schreiner W.S., Sokolovskiy S.V., Rocken C., Hunt D.C. (1999),
%   "Analysis and validation of GPS/MET radio occultation data in the
%   ionosphere", Radio Science, 34(4), 949-966.
% - IS-GPS-200, Section 20.3.3.4.3 (Earth rotation / Sagnac correction).
%
% -------------------------------------------------------------------------
%   AUTHOR:
%       Francesco Fiorina, NavSAS - DET - Politecnico di Torino
%
%   VERSION HISTORY:
%       1.0 - 21/04/2026  - Straight-line geometric SLTA, more robust
%                           at lunar altitudes. GroundTruth passed
%                           explicitly. Sign of the tangent-point
%                           projection corrected (foot of perpendicular
%                           requires s* = -r_rx . t_hat). Output
%                           fields pre-allocated. Geometric self-
%                           consistency check on the tangent point.
%       1.1 - 22/04/2026  - Closest-ephemeris pre-selection via
%                           selectClosestEphemeris before propagation.
%                           GroundTruth loaded from settings.gtPath.
%                           Two visualisation blocks added:
%                           PlotGrazingAltitude (raw) and
%                           PlotGrazingAltitudeArcs (observable-masked
%                           with tracking-point retrieval).
%       1.2 - 04/06/2026  - Light-time correction: SV state propagated
%                           at tTX = tRX - rawP/c instead of tRX.
%                           Sagnac rotation applied to bring r_GNSS
%                           and v_GNSS into the ECEF(tRX) frame.
%                           Epochs without any valid rawP (neither
%                           L1 nor L5) are skipped. Grazing altitude
%                           now returned as WGS-84 geodetic height
%                           of the tangent point (via ecef2geodetic),
%                           replacing the previous spherical-Earth
%                           approximation.
%
%   COPYRIGHT:
%       NavSAS - DET - Politecnico di Torino, 2026. All rights reserved.
% -------------------------------------------------------------------------

% --- Command Window Log
fprintf('-> Computing Grazing altitude and tangent point...\n')

%% === RINEX navigation data ===
rinex               = rinexread(settings.rnxPath);
epochs              = length(timeUTC);
earth_radius        = physconst('earthRadius');

wgs84 = wgs84Ellipsoid;
earth_semiminorAxis = wgs84.SemiminorAxis; deltaMinor = earth_semiminorAxis - earth_radius;
earth_semimajorAxis = wgs84.SemimajorAxis; deltaMajor = earth_semimajorAxis - earth_radius;


%% === Ground Truth ===
% GroundTruth : struct with fields
%     .pos_ecef : [N_epochs x 3] receiver ECEF position [m]
%     .vel_ecef : [N_epochs x 3] receiver ECEF velocity [m/s]
%                 (not used by this function but kept for interface
%                  consistency with computeBendingAngle)
tmp = load(settings.gtPath); GroundTruth = tmp.GroundTruth;

%% === Pre-allocate output fields ===
nSV = numel(obsSolutions);
for ii = 1 : nSV
    obsSolutions(ii).SVpos            = nan(epochs, 3);
    obsSolutions(ii).SVvel            = nan(epochs, 3);
    obsSolutions(ii).grazingAltitude  = nan(epochs, 1);
    obsSolutions(ii).tangentPoint     = nan(epochs, 3);
end

%% === Constellation loop ===
for const = 1 : 2

    if     const == 1 && settings.G
        navUse = rinex.GPS;
        tag    = 'G';
    elseif const == 2 && settings.E
        navUse = rinex.Galileo;
        tag    = 'E';
    else
        continue
    end

    %% === Epoch loop ===
%     for epoch = 1 : epochs
% 
%         tQuery = timeUTC(epoch);
% 
%         % Propagate the whole constellation at this epoch
%         navClosest = selectClosestEphemeris(navUse, tQuery);
%         [pos, vel, ids] = gnssconstellation(tQuery, navClosest);
% 
%         for ii = 1 : numel(ids)
%             prn = ids(ii);
% 
%             % Match this SV to one or more rows of obsSolutions (one
%             % per tracked frequency)
%             idx = find(strcat(tag, string(prn))== [obsSolutions.SV]);
%             if isempty(idx), continue; end
% 
%             r_GNSS = pos(ii, :);
%             v_GNSS = vel(ii, :);
%             r_rx   = GroundTruth.pos_ecef(epoch, :);
% 
%             % --- Straight-line grazing altitude ---------------------
%             % Vector along the ray, rx -> GNSS
%             r_tangent = r_GNSS - r_rx;
% 
%             % sin(phi) between r_rx and the ray, from the cross-product
%             % identity |a x b| = |a||b|sin(phi)
%             sinPhi = norm(cross(r_rx, r_tangent)) ...
%                      / (norm(r_rx) * norm(r_tangent));
% 
%             % Distance from Earth centre to the straight ray
%             altitude         = norm(r_rx) * sinPhi;
%             grazing_altitude = altitude - earth_radius;
% 
%             % --- Tangent point --------------------------------------
%             % Foot of perpendicular from the origin onto the line
%             % r(s) = r_rx + s * t_hat, i.e.
%             %   s* = -(r_rx . t_hat)
%             %   tangent_point = r_rx + s* * t_hat
%             r_tangent_hat = r_tangent / norm(r_tangent);
%             s_star        = -dot(r_rx, r_tangent_hat);
%             tangent_point = r_rx + s_star * r_tangent_hat;
% 
%             % --- Geometric sanity check -----------------------------
%             % By construction the tangent point is the point on the
%             % ray closest to Earth centre, so |tangent_point| must
%             % equal the straight-line impact parameter 'altitude'.
%             % The two are computed by independent algebraic paths
%             % (cross-product identity vs. projection onto t_hat), so
%             % a mismatch beyond floating-point rounding indicates a
%             % degenerate or ill-conditioned geometry.
%             tol = 1e-6 * altitude;   % ~1 ppm relative tolerance
%             if abs(norm(tangent_point) - altitude) > tol
%                 warning('grazingAltitude:InconsistentTangentPoint', ...
%                     ['Inconsistent tangent-point geometry at epoch ' ...
%                      '%d, SV %s%02d: |tangent_point| = %.3f m differs ' ...
%                      'from altitude = %.3f m by %.3e m. ' ...
%                      'Setting grazing_altitude and tangent_point to NaN.'], ...
%                     epoch, tag, prn, ...
%                     norm(tangent_point), altitude, ...
%                     abs(norm(tangent_point) - altitude));
%                 grazing_altitude = NaN;
%                 tangent_point    = [NaN NaN NaN];
%             end
% 
%             % --- Write to all matched obsSolutions rows -------------
%             for f = 1 : numel(idx)
%                 obsSolutions(idx(f)).SVpos(epoch, :)         = r_GNSS;
%                 obsSolutions(idx(f)).SVvel(epoch, :)         = v_GNSS;
%                 obsSolutions(idx(f)).grazingAltitude(epoch)  = grazing_altitude;
%                 obsSolutions(idx(f)).tangentPoint(epoch, :)  = tangent_point;
%             end
%         end
%     end
    for epoch = 1 : epochs
    
        tRX    = timeUTC(epoch);
        c      = physconst('LightSpeed');
        omegaE = 7.2921151467e-5;             % WGS-84 Earth rotation rate [rad/s]
    
        navClosest  = selectClosestEphemeris(navUse, tRX);
    
        % Visible PRN list (states at tRX used only to identify which SVs
        % have a valid ephemeris at this epoch)
        [~, ~, ids] = gnssconstellation(tRX, navClosest);
    
        for ii = 1 : numel(ids)
            prn = ids(ii);
    
            % Match this SV to its obsSolutions rows
            idx = find(strcat(tag, string(prn)) == [obsSolutions.SV]);
            if isempty(idx), continue; end
    
            % --- Pseudorange-based time of flight ----------------------
            % Use the first non-NaN rawP across L1/L5. If both are NaN at
            % this epoch the SV is not being tracked, so skip.
            rawP_epoch = NaN;
            for f = 1 : numel(idx)
                candidate = obsSolutions(idx(f)).rawP(epoch);
                if ~isnan(candidate)
                    rawP_epoch = candidate;
                    break
                end
            end
            if isnan(rawP_epoch), continue; end
    
            dt  = rawP_epoch / c;            % [s], includes minor biases
            tTX = tRX - seconds(dt);
    
            % --- Propagate this PRN at transmission time --------------
            navOnePRN = navClosest(navClosest.SatelliteID == prn, :);
            if isempty(navOnePRN), continue; end
            [p, v, ~] = gnssconstellation(tTX, navOnePRN);
    
            % --- Sagnac: rotate from ECEF(tTX) to ECEF(tRX) -----------
            theta_E = omegaE * dt;
            R       = [ cos(theta_E)  sin(theta_E)  0
                       -sin(theta_E)  cos(theta_E)  0
                                  0             0   1 ];
            r_GNSS  = (R * p(1, :).').';
            v_GNSS  = (R * v(1, :).').';
    
            % --- Geometry ---------------------------------------------
            r_rx      = GroundTruth.pos_ecef(epoch, :);
            % Vector along the ray, rx -> GNSS
            r_tangent = r_GNSS - r_rx;
    
            % sin(phi) between r_rx and the ray, from the cross-product
            % identity |a x b| = |a||b|sin(phi)
            sinPhi    = norm(cross(r_rx, r_tangent)) ...
                        / (norm(r_rx) * norm(r_tangent));

            % Distance from Earth centre to the straight ray    
            altitude                = norm(r_rx) * sinPhi;
            grazing_altitude_radius = altitude - earth_radius;
  
            % --- Tangent point --------------------------------------
            % Foot of perpendicular from the origin onto the line
            % r(s) = r_rx + s * t_hat, i.e.
            %   s* = -(r_rx . t_hat)
            %   tangent_point = r_rx + s* * t_hat
            r_tangent_hat = r_tangent / norm(r_tangent);
            s_star        = -dot(r_rx, r_tangent_hat);
            tangent_point = r_rx + s_star * r_tangent_hat;

            % Tangent point in geodetic coordinates (WGS-84 ellipsoid)
            lla = ecef2lla(tangent_point);          % [lat_deg, lon_deg, h_m]
            grazing_altitude = lla(3);              % geodetic height above WGS-84

            % TOOLBOX AGNOSTIC VERSION
            % wgs84 = wgs84Ellipsoid('meter');
            % [~, ~, grazing_altitude] = ecef2geodetic(wgs84, ...
            %    tangent_point(1), tangent_point(2), tangent_point(3));
    
            % --- Geometric sanity check -----------------------------
            % By construction the tangent point is the point on the
            % ray closest to Earth centre, so |tangent_point| must
            % equal the straight-line impact parameter 'altitude'.
            % The two are computed by independent algebraic paths
            % (cross-product identity vs. projection onto t_hat), so
            % a mismatch beyond floating-point rounding indicates a
            % degenerate or ill-conditioned geometry.
            tol = 1e-6 * altitude;
            if abs(norm(tangent_point) - altitude) > tol
                warning('grazingAltitude:InconsistentTangentPoint', ...
                    ['Inconsistent tangent-point geometry at epoch %d, ' ...
                     'SV %s%02d: |tangent_point| = %.3f m differs from ' ...
                     'altitude = %.3f m by %.3e m. Setting NaN.'], ...
                    epoch, tag, prn, ...
                    norm(tangent_point), altitude, ...
                    abs(norm(tangent_point) - altitude));
                grazing_altitude = NaN;
                tangent_point    = [NaN NaN NaN];
            end
    
            % --- Write to all matched obsSolutions rows ----------------
            for f = 1 : numel(idx)
                obsSolutions(idx(f)).SVpos(epoch, :)                        = r_GNSS;
                obsSolutions(idx(f)).SVvel(epoch, :)                        = v_GNSS;
                obsSolutions(idx(f)).grazingAltitude(epoch)                 = grazing_altitude;
                obsSolutions(idx(f)).grazingAltitude_Rearth_mean(epoch)     = grazing_altitude_radius;
                obsSolutions(idx(f)).tangentPoint(epoch, :)                 = tangent_point;
            end
        end
    end
end
%% -----------------------------------------------------------------------
% PLOT Grazing Altitudes
% -----------------------------------------------------------------------
plotSwitch = settings.PlotGrazingAltitude;
if plotSwitch

    % Create figure
    h = figure();
    figurename = 'GrazingAltitude';
    set(h, 'Name', figurename);
    hold on;

    % Track satellites already drawn on this axes to avoid duplicates
    % (each SV appears once per tracked frequency in obsSolutions)
    ax          = gca;
    existingSVs = string({});
    if ~isempty(ax.Children)
        existingSVs = string(get(ax.Children, 'DisplayName'));
    end
    plottedSVs = existingSVs;

    % Plot grazing altitude for each satellite (km for readability)
    for i = 1:length(obsSolutions)
        svID = obsSolutions(i).SV;

        % Skip if this SV is already on the plot
        if any(plottedSVs == string(svID))
            continue
        end

        RGBcolor = PRNColors(svID);

        plot(timeUTC, obsSolutions(i).grazingAltitude/1e3, ...
            'Color',       RGBcolor, ...
            'DisplayName', svID, ...
            'LineWidth',   2);

        plottedSVs(end+1) = string(svID); 
    end

    % Reference lines: Earth surface and ionosphere base
    yline(0,  '--k', 'Earth surface', 'Interpreter','latex', ...
          'LabelHorizontalAlignment','left','HandleVisibility','off');
    yline(deltaMinor/1e3,  '--k', 'Earth Semiminor Axis', 'Interpreter','latex', ...
          'LabelHorizontalAlignment','left','HandleVisibility','off');
    yline(deltaMajor/1e3,  '--k', 'Earth Semimajor Axis', 'Interpreter','latex', ...
          'LabelHorizontalAlignment','left','HandleVisibility','off');
    yline(80, ':k',  'Ionosphere base (80 km)', 'Interpreter','latex', ...
          'LabelHorizontalAlignment','left','HandleVisibility','off');
    yline(350, '-.k', 'F2 peak (~350 km)',       'Interpreter','latex', ...
          'LabelHorizontalAlignment','left','HandleVisibility','off');

    % Axis labels
    xlabel('Time', 'Interpreter','latex');
    ylabel('Grazing altitude [km]', 'Interpreter','latex');
    axis tight
    grid on

    % Subtitle with time interval (built directly from timeUTC)
    subtitle(strcat("Start: ", datestr(timeUTC(1)), ...
                    "  |  End: ", datestr(timeUTC(end))), ...
             'Interpreter','latex');

    % Formatting
    set(gca, 'FontSize', 13, 'TickLabelInterpreter', 'latex');

    % Datatips
    dcm              = datacursormode(gcf);
    dcm.Enable       = 'on';
    dcm.DisplayStyle = 'datatip';
    dcm.UpdateFcn    = @customDataTip;

    % Legend
    hLeg = legend('Location','eastoutside', 'NumColumns',2, ...
                  'FontSize',8, 'Interpreter','latex');
    hLeg.ItemHitFcn = @action1;

    % Title and window placement
    title('Grazing altitude - geometric', 'Interpreter','latex');
    Pix_SS = get(0,'screensize');
    h.OuterPosition = [Pix_SS(3)/2, 0, Pix_SS(3)/2, Pix_SS(4)/2];

   % cleanLegend(figurename)
end

%% -----------------------------------------------------------------------
% PLOT Grazing Altitudes Arcs
% -----------------------------------------------------------------------
plotSwitch = settings.PlotGrazingAltitudeArcs;
if plotSwitch

    % Create figure
    h = figure();
    figurename = 'GrazingAltitudeArcs';
    set(h, 'Name', figurename);
    hold on;

    % Track satellites already drawn on this axes to avoid duplicates
    % (each SV appears once per tracked frequency in obsSolutions)
    ax          = gca;
    existingSVs = string({});
    if ~isempty(ax.Children)
        existingSVs = string(get(ax.Children, 'DisplayName'));
    end
    plottedSVs = existingSVs;

    % Plot grazing altitude for each satellite (km for readability),
    % masked to the epochs where ROTI is defined
    for i = 1:length(obsSolutions)
        svID = obsSolutions(i).SV;

        % Skip if SV does not exists in double frequency
        if numel(find(svID == [obsSolutions.SV])) < 2
            continue
        end

        % Skip if this SV is already on the plot
        if any(plottedSVs == string(svID))
            continue
        end

        % Find all obsSolutions rows for this SV
        sameSVrows = find(string(svID) == [obsSolutions.SV]);

        % Build a "tracking-valid" mask: keep epochs where at least one of
        % the SV's frequencies has valid pseudorange (rawP). Use & to require
        % BOTH missing before NaNing; use | to require ALL present.
        if numel(sameSVrows) >= 2
            bothMissing = isnan(obsSolutions(sameSVrows(1)).rawP) & ...
                isnan(obsSolutions(sameSVrows(end)).rawP);
        else
            bothMissing = isnan(obsSolutions(sameSVrows(1)).rawP);
        end

        y = obsSolutions(i).grazingAltitude / 1e3;
        y(bothMissing) = NaN;

        % Skip SVs with no valid ROTI at all
        if all(isnan(y))
            continue
        end

        % Retrieve tracking point start and end
        trackingPoints = retrieveTrackingPoints(y*1e3,settings);
        if numel(sameSVrows) >= 2
            obsSolutions(sameSVrows(1)).SVtrackingPoints = trackingPoints;
            obsSolutions(sameSVrows(end)).SVtrackingPoints = trackingPoints; 
        else
            obsSolutions(sameSVrows(1)).SVtrackingPoints = trackingPoints;
        end

        % Plot

        RGBcolor = PRNColors(svID);

        plot(timeUTC, y, ...
            'Color',       RGBcolor, ...
            'DisplayName', svID, ...
            'LineWidth',   2);

        plottedSVs(end+1) = string(svID);
    end

    % Reference lines: Earth surface, ionosphere base, F2 peak
    yline(0,   '--k', 'Earth surface',           'Interpreter','latex', ...
          'LabelHorizontalAlignment','left','HandleVisibility','off');
    yline(deltaMinor/1e3,  '--k', 'Earth Semiminor Axis', 'Interpreter','latex', ...
          'LabelHorizontalAlignment','left','HandleVisibility','off');
    yline(deltaMajor/1e3,  '--k', 'Earth Semimajor Axis', 'Interpreter','latex', ...
          'LabelHorizontalAlignment','left','HandleVisibility','off');
    yline(80,  ':k',  'Ionosphere base (80 km)', 'Interpreter','latex', ...
          'LabelHorizontalAlignment','left','HandleVisibility','off');
    yline(350, '-.k', 'F2 peak (~350 km)',       'Interpreter','latex', ...
          'LabelHorizontalAlignment','left','HandleVisibility','off');

    % Axis labels
    xlabel('Time', 'Interpreter','latex');
    ylabel('Grazing altitude [km]', 'Interpreter','latex');
    axis tight
    grid on

    % Subtitle with time interval (built directly from timeUTC)
    subtitle(strcat("Start: ", datestr(timeUTC(1)), ...
                    "  |  End: ", datestr(timeUTC(end))), ...
             'Interpreter','latex');

    % Formatting
    set(gca, 'FontSize', 13, 'TickLabelInterpreter', 'latex');

    % Datatips
    dcm              = datacursormode(gcf);
    dcm.Enable       = 'on';
    dcm.DisplayStyle = 'datatip';
    dcm.UpdateFcn    = @customDataTip;

    % Legend
    hLeg = legend('Location','eastoutside', 'NumColumns',2, ...
                  'FontSize',8, 'Interpreter','latex');
    hLeg.ItemHitFcn = @action1;

    % Title and window placement
    title('Grazing altitude (Observables valid arcs only)', 'Interpreter','latex');
    Pix_SS = get(0,'screensize');
    h.OuterPosition = [Pix_SS(3)/2, 0, Pix_SS(3)/2, Pix_SS(4)/2];

    cleanLegend(figurename)
end

% --- Command Window Log
fprintf('Grazing altitude and tangent point computation ended. \n\n')


end