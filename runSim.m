%% AE 6505 Final Project - Patrick Barry
% File to run full simulation with extended Kalman filter

clear; clc; close all;
rng(700);
%% Set Options

% Dr. Gunter, if you're reviewing this code, I'd set this the below to true, it's
% pretty cool even if the animation takes > 1 minute to fully draw out
showSetupPlot = false; % set to true to see cool (but long-animated) plot of the missile trajectory and constellation design!
customConstellationDesign = false; % set to true to specifically manipulate constellation parameters in "Build constellation" section
constellationChoice = 3; % either 1, 2, or 3 - irrelevant if customConstellationDesign is true
    % constellationChoice = 1: Space-Based Infrared System (SBIRS)
    % constellationChoice = 2: Space Development Agency Tranche 1 Constellation
    % constellationChoice = 3: Self-made MEO constellation Design

%% Load missile trajectory information
load('missileTraj.mat')
% Loads the following:
% - xs_misile (6xN), time history of r, v, of ballistic missile
% - ts_missile (N), time record for xs_missile
% - phaseOfFlight (N), cell array denoting which phase of flight the
% missile is in, either boost, midcourse, or reentry
% - betas (N), apparent brightness coefficient for IR sensors that missile
% gives off in each stage of flight
% - r_launchSite_ECI, location of launch site in ECI frame
% - r_target_ECI, location of target in ECI frame

switch constellationChoice
    case 1
        load("constellation1.mat")
    case 2
        load("constellation2.mat")
    case 3
        load("constellation3.mat")
    otherwise
        fprintf("Custom constellation design selected, make sure you set parameters under ""Build Constellation"" section.\n");
end



%% Constants

n = 6; % number of states
m = 4; % number of measurements per satellite

RE = 6371000; % radius of Earth, m
dt = ts_missile(2) - ts_missile(1); % align missile sim timestep with this sim timestep

minObservableMissileHeight = 100*1000;

measCovs.az = 0.5/100; % rad
measCovs.el = 0.5/100; % rad
measCovs.azrate = 0.5/1000; % rad/s
measCovs.elrate = 0.5/1000; % rad/s

processCovs.xvel = 0;
processCovs.yvel = 0;
processCovs.zvel = 0;
processCovs.xaccel = 5;
processCovs.yaccel = 5;
processCovs.zaccel = 5;

Qc = diag([processCovs.xvel^2, processCovs.yvel^2, processCovs.zvel^2, ...
    processCovs.xaccel^2, processCovs.yaccel^2, processCovs.zaccel^2]);

Q = Qc*dt;

% For plotting Earth
[xe, ye, ze] = sphere(50); % arg is resolution
xe = RE*xe;
ye = RE*ye;
ze = RE*ze;

%% Constellation Design

if customConstellationDesign
    SatelliteArray = [];
    satID = 1;
    
    % Constellation 1 Information
    constellationName = 'MEO_1';
    N_constellation = 4;
    a = 10000*1000 + RE; % semimajor axis, m 
    e = 0.0; % eccentricity
    inc = -60; % inclination (deg)
    O = 0;
    w = 0;
    fs = spreadSatellitesEvenly(N_constellation); 
    
    for i = 1:N_constellation
        sat = Satellite(satID, constellationName, a, e, inc, O, w, fs(i));
        SatelliteArray = [SatelliteArray; sat];
    
        satID = satID + 1;
    
    end
    
    % Constellation 2 Information
    constellationName = 'MEO_2';
    N_constellation = 4;
    a = 10000*1000 + 6371000; % semimajor axis, m
    e = 0.0; % eccentricity
    inc = 0; % inclination (deg)
    O = 0;
    w = 0;
    fs = spreadSatellitesEvenly(N_constellation); 
    
    for i = 1:N_constellation
        sat = Satellite(satID, constellationName, a, e, inc, O, w, fs(i));
        SatelliteArray = [SatelliteArray; sat];
    
        satID = satID + 1;
    
    end

    % Constellation 3 Information
    constellationName = 'MEO_3';
    N_constellation = 4;
    a = 10000*1000 + 6371000; % semimajor axis, m
    e = 0.0; % eccentricity
    inc = 60; % inclination (deg)
    O = 0;
    w = 0;
    fs = spreadSatellitesEvenly(N_constellation); 
    
    for i = 1:N_constellation
        sat = Satellite(satID, constellationName, a, e, inc, O, w, fs(i));
        SatelliteArray = [SatelliteArray; sat];
    
        satID = satID + 1;
    
    end


    max_viewdist = 10000*1000 + 6371000;
end
fprintf("Constellation initialized!\n")

%% Show simulation setup plot

if showSetupPlot
    figure(1); clf
    h_launchSite = scatter3(r_launchSite_ECI(1), r_launchSite_ECI(2), r_launchSite_ECI(3), 80, 'filled');
    hold on
    h_target = scatter3(r_target_ECI(1), r_target_ECI(2), r_target_ECI(3), 80, 'filled');
    h_missile = scatter3(xs_missile(1,1), xs_missile(2,1), xs_missile(3,1), 80, 'filled');
    h_sats = [];
    for i = 1 :length(SatelliteArray)
        h_sats(i) = scatter3(SatelliteArray(i).x(1), SatelliteArray(i).x(2), SatelliteArray(i).x(3), 80, 'magenta','filled');
    end
    h_Earth = surf(xe, ye, ze, ...
        'FaceColor', [0.2 0.6 1], ...
        'EdgeColor', 'none', ...
        'FaceAlpha', 0.5); % semi-transparent
    h_traj_missile = plot3(NaN, NaN, NaN, 'r-');
    h_sat_trajs = [];
    for i = 1:length(SatelliteArray)
        h_sat_trajs(i) = plot3(NaN, NaN, NaN, 'g-');
    end
    axis equal
    xlim([-max_viewdist max_viewdist])
    ylim([-max_viewdist max_viewdist])
    zlim([-max_viewdist max_viewdist])
    view([45, 20])
    grid on
    xlabel('ECI x (m)', 'FontSize', 20);
    ylabel('ECI y (m)', 'FontSize', 20);
    zlabel('ECI z (m)', 'FontSize', 20);
    title('Simulation Setup', 'FontSize', 20)
    legend([h_launchSite, h_target, h_missile, h_sats], {'Launch Site', 'Target', 'Missile', 'Satellites'})
    
    for k = 1:length(xs_missile)
        set(h_missile, 'XData', xs_missile(1,k), 'YData', xs_missile(2,k), 'ZData', xs_missile(3,k));
        set(h_traj_missile, 'XData', xs_missile(1,1:k), 'YData', xs_missile(2,1:k), 'ZData', xs_missile(3,1:k));
    
        for i = 1:length(SatelliteArray)
            SatelliteArray(i).propagateSatellite(dt);
    
            set(h_sats(i), 'XData', SatelliteArray(i).x(1), 'YData', SatelliteArray(i).x(2), 'ZData', SatelliteArray(i).x(3))
            set(h_sat_trajs(i), ...
                'XData', [get(h_sat_trajs(i),'XData'), SatelliteArray(i).x(1)], ...
                'YData', [get(h_sat_trajs(i),'YData'), SatelliteArray(i).x(2)], ...
                'ZData', [get(h_sat_trajs(i),'ZData'), SatelliteArray(i).x(3)]);
        end
    
        drawnow;
    end
end

%% Run the simulation
fprintf("\nRunning simulation!\n")
attemptEstimation = false; % flag for whether Kalman filter should be run based on missile altitude
estimationBegan = false; % flag for whether estimation has occurred
timeEstimationBegan = nan;
missileStillVisible = true;
timeEstimationEnded = nan;

% Allocate for KF matrices
xs_estimated = nan(6, length(ts_missile));
Ps_estimated = nan(length(ts_missile), 6, 6);


% Allocate for metrics
NObserving = zeros(length(ts_missile), 1); % Number of observing satellites at any given timestep
isObservable = zeros(length(ts_missile), 1); % True/False of whether missile is observable based on its altitude


for i = 1:length(ts_missile) % Timescale is based off of missile time, since this whole sim revolves around missile flight

    % Update missile truth state
    x_missile = xs_missile(:,i);
    beta = betas(i);
    flightPhase = phaseOfFlight{i};

    % Update satellite(s) truth state and observability of missile
    NObs_i = 0;
    for j = 1:length(SatelliteArray)
        SatelliteArray(j).propagateSatellite(dt); % TODO: ADD NOISE TO SATELLITE PROPAGATION?
        SatelliteArray(j).checkObservability(x_missile(1:3), beta);
        if SatelliteArray(j).canObserve
            NObs_i = NObs_i + 1;
        end
    end
    NObserving(i) = NObs_i;

    % Determine if estimation filter should be run based on satellite
    % altitude
    if norm(x_missile(1:3)) > RE + minObservableMissileHeight
        attemptEstimation = true;
        isObservable(i) = 1;
    else
        attemptEstimation = false;
        NObserving(i) = 0;
    end

    % State estimation doesn't begin until satellite is past observable
    % altitude and at least one satellite observes it
    if ~attemptEstimation && ~estimationBegan % missile cannot be viewed and estimation has not yet begun

        % Hold estimated state at arbitary zero while covariance is set to
        % be very high
        xs_estimated(:,i) = nan;
        Ps_estimated(i,:,:) = nan;

    elseif attemptEstimation && NObserving(i) > 0 && ~estimationBegan % 1st time missile can be viewed and is being viewed
        estimationBegan = true;
        timeEstimationBegan = ts_missile(i);

        % Filter initialization, or the first time the filter actually
        % estimates

        % Get simulated measurement yk (what each satellite sees)
        yk = formSimulatedMeasurement(SatelliteArray, measCovs, x_missile);

        % Initialize EKF based on line of sight (NObs = 1) or triangulation
        % (NObs >= 2)
        [x0, P0] = initializeEKF(yk, SatelliteArray);
        xs_estimated(:, i) = x0;
        Ps_estimated(i,:,:) = P0;

    elseif attemptEstimation && NObserving(i) > 0 && estimationBegan % missile can be viewed and is being viewed with prior info
        % Classic EKF

        % Time update
        A = computeA(xs_estimated(:,i-1));
        Phi = computeSTM(A, dt, 2); % compute STM to 2nd order
        x_prior = propagateMissile(xs_estimated(:,i-1), dt);
        P_prior = Phi*squeeze(Ps_estimated(i-1,:,:))*Phi' + Q;

        % Get simulated measurement yk (what each satellite sees)
        yk = formSimulatedMeasurement(SatelliteArray, measCovs, x_missile);

        % Measurement update
        H = formH(SatelliteArray, x_prior);
        R = formR(measCovs, NObserving(i), m);
        K = P_prior*H'*inv(H*P_prior*H' + R);
        innov = yk - formPredictedMeasurement(SatelliteArray, x_prior);
        for k = 1:NObserving(i) % Ensure angles lie in the +/- pi range to avoid massive residuals (analogous to frisbee problem)
            az_idx = 1 + (k-1)*4;
            el_idx = 2 + (k-1)*4;
            innov(az_idx) = wrapToPi(innov(az_idx));
            innov(el_idx) = wrapToPi(innov(el_idx));
        end
        xs_estimated(:,i) = x_prior + K*innov;
        Ps_estimated(i,:,:) = (eye(n) - K*H)*P_prior*(eye(n)-K*H)' + K*R*K';

    elseif attemptEstimation && NObserving(i) == 0 && estimationBegan % missile can be viewed but no satellites can view it, and estimation already began --> just propagate state/cov
        timeEstimationEnded = ts_missile(i) - dt;
        % Pure propagation

        A = computeA(xs_estimated(:,i-1));
        Phi = computeSTM(A, dt, 2); % compute STM to 2nd order
        xs_estimated(:,i) = propagateMissile(xs_estimated(:,i-1), dt);
        Ps_estimated(i, :, :) = Phi*squeeze(Ps_estimated(i-1,:,:))*Phi' + Q;

    elseif ~attemptEstimation & estimationBegan %% missile cannot be viewed but estimation already began --> propagate with reentry dynamics
        
        if missileStillVisible
            timeEstimationEnded = ts_missile(i);
            missileStillVisible = false;
        end
        
        % Pure propagation

        A = computeAReentry(xs_estimated(:,i-1));
        Phi = computeSTM(A, dt, 2); % compute STM to 2nd order
        xs_estimated(:,i) = propagateMissileReentry(xs_estimated(:,i-1), dt);
        Ps_estimated(i, :, :) = Phi*squeeze(Ps_estimated(i-1,:,:))*Phi' + Q;
    end

end

%% Plots/Results

timeBoostEnd = nan;
timeReentryBegin = nan;
for i = 1:length(phaseOfFlight)
    if strcmp(phaseOfFlight{i}, 'MIDCOURSE') && isnan(timeBoostEnd)
        timeBoostEnd = ts_missile(i-1);
    elseif strcmp(phaseOfFlight{i}, 'REENTRY') && isnan(timeReentryBegin)
        timeReentryBegin = ts_missile(i);
    end
end

threeSigmas = processCovarianceToSigmas(Ps_estimated, ts_missile);

% Number of satellites observing at each time, and whether
% missile can be observed
figure(2); clf
plot(1:timeEstimationBegan-1, NObserving(1:timeEstimationBegan-1), 'r', 'LineWidth', 4)
hold on
plot(timeEstimationBegan:timeEstimationEnded+1, NObserving(timeEstimationBegan:timeEstimationEnded+1), 'g', 'LineWidth', 4)
plot(timeEstimationEnded+1:length(ts_missile), NObserving(timeEstimationEnded+1:end), 'r', 'LineWidth', 4)
ylabel('Number of Satellites Observing Missile')
xlabel('Time (s)')
legend('Missile Below Visible Altitude', 'Missile Above Missile Altitude', 'Location', 'best')
grid on
box on
title('Missile Observability During Flight')

% Missile Estimated and True Trajectory plot (position) with phases of
% flight and when estimation began/ended
figure(3); clf
subplot(3,2,1)
plot(ts_missile, xs_missile(1,:), 'LineWidth', 2)
hold on
plot(ts_missile, xs_estimated(1,:), 'LineWidth', 2)
ylabel('ECI x Position (m)')
xlabel('Time (s)')
plot([timeEstimationEnded timeEstimationEnded], [-1e10 1e10], 'k--')
plot([timeBoostEnd timeBoostEnd], [-1e10 1e10], 'magenta--')
plot([timeReentryBegin timeReentryBegin], [-1e10 1e10], 'magenta--')
plot([timeEstimationBegan timeEstimationBegan], [-1e10 1e10], 'k--')
legend('Actual', 'Estimated', 'Estimation Bounds', 'Midcourse Phase Bounds', 'Location', 'best')
ylim_low = min([xs_estimated(1,:), xs_missile(1,:)]);
ylim_high = max([xs_estimated(1,:), xs_missile(1,:)]);
ylim([ylim_low ylim_high])
grid on
box on

subplot(3,2,2)
plot(ts_missile, xs_missile(4,:), 'LineWidth', 2)
hold on
plot(ts_missile, xs_estimated(4,:), 'LineWidth', 2)
ylabel('ECI x Velocity (m/s)')
xlabel('Time (s)')
plot([timeEstimationEnded timeEstimationEnded], [-1e10 1e10], 'k--')
plot([timeBoostEnd timeBoostEnd], [-1e10 1e10], 'magenta--')
plot([timeReentryBegin timeReentryBegin], [-1e10 1e10], 'magenta--')
plot([timeEstimationBegan timeEstimationBegan], [-1e10 1e10], 'k--')
legend('Actual', 'Estimated', 'Estimation Bounds', 'Midcourse Phase Bounds', 'Location', 'best')
ylim_low = min([xs_estimated(4,:), xs_missile(4,:)]);
ylim_high = max([xs_estimated(4,:), xs_missile(4,:)]);
ylim([ylim_low ylim_high])
grid on
box on


subplot(3,2,3)
plot(ts_missile, xs_missile(2,:), 'LineWidth', 2)
hold on
plot(ts_missile, xs_estimated(2,:), 'LineWidth', 2)
ylabel('ECI y Position (m)')
xlabel('Time (s)')
plot([timeEstimationEnded timeEstimationEnded], [-1e10 1e10], 'k--')
plot([timeBoostEnd timeBoostEnd], [-1e10 1e10], 'magenta--')
plot([timeReentryBegin timeReentryBegin], [-1e10 1e10], 'magenta--')
plot([timeEstimationBegan timeEstimationBegan], [-1e10 1e10], 'k--')
legend('Actual', 'Estimated', 'Estimation Bounds', 'Midcourse Phase Bounds', 'Location', 'best')
ylim_low = min([xs_estimated(2,:), xs_missile(2,:)]);
ylim_high = max([xs_estimated(2,:), xs_missile(2,:)]);
ylim([ylim_low ylim_high])
grid on
box on

subplot(3,2,4)
plot(ts_missile, xs_missile(5,:), 'LineWidth', 2)
hold on
plot(ts_missile, xs_estimated(5,:), 'LineWidth', 2)
ylabel('ECI y Velocity (m/s)')
xlabel('Time (s)')
legend('Estimated', 'Actual')
plot([timeEstimationEnded timeEstimationEnded], [-1e10 1e10], 'k--')
plot([timeBoostEnd timeBoostEnd], [-1e10 1e10], 'magenta--')
plot([timeReentryBegin timeReentryBegin], [-1e10 1e10], 'magenta--')
plot([timeEstimationBegan timeEstimationBegan], [-1e10 1e10], 'k--')
legend('Actual', 'Estimated', 'Estimation Bounds', 'Midcourse Phase Bounds', 'Location', 'best')
ylim_low = min([xs_estimated(5,:), xs_missile(5,:)]);
ylim_high = max([xs_estimated(5,:), xs_missile(5,:)]);
ylim([ylim_low ylim_high])
grid on
box on

subplot(3,2,5)
plot(ts_missile, xs_missile(3,:), 'LineWidth', 2)
hold on
plot(ts_missile, xs_estimated(3,:), 'LineWidth', 2)
ylabel('ECI z Position (m)')
xlabel('Time (s)')
plot([timeEstimationEnded timeEstimationEnded], [-1e10 1e10], 'k--')
plot([timeBoostEnd timeBoostEnd], [-1e10 1e10], 'magenta--')
plot([timeReentryBegin timeReentryBegin], [-1e10 1e10], 'magenta--')
plot([timeEstimationBegan timeEstimationBegan], [-1e10 1e10], 'k--')
legend('Actual', 'Estimated', 'Estimation Bounds', 'Midcourse Phase Bounds', 'Location', 'best')
ylim_low = min([xs_estimated(3,:), xs_missile(3,:)]);
ylim_high = max([xs_estimated(3,:), xs_missile(3,:)]);
ylim([ylim_low ylim_high])
grid on
box on

subplot(3,2,6)
plot(ts_missile, xs_missile(6,:), 'LineWidth', 2)
hold on
plot(ts_missile, xs_estimated(6,:), 'LineWidth', 2)
ylabel('ECI z Velocity (m/s)')
xlabel('Time (s)')
plot([timeEstimationEnded timeEstimationEnded], [-1e10 1e10], 'k--')
plot([timeBoostEnd timeBoostEnd], [-1e10 1e10], 'magenta--')
plot([timeReentryBegin timeReentryBegin], [-1e10 1e10], 'magenta--')
plot([timeEstimationBegan timeEstimationBegan], [-1e10 1e10], 'k--')
legend('Actual', 'Estimated', 'Estimation Bounds', 'Midcourse Phase Bounds', 'Location', 'best')
ylim_low = min([xs_estimated(6,:), xs_missile(6,:)]);
ylim_high = max([xs_estimated(6,:), xs_missile(6,:)]);
ylim([ylim_low ylim_high])
grid on
box on

sgtitle('True and Estimated Trajectories')

%  Estimation error with covariance bounds
figure(4); clf

subplot(3,2,1)
plot(ts_missile, xs_missile(1,:) - xs_estimated(1,:))
hold on
plot(ts_missile, threeSigmas(1,:), 'k-')
plot([timeEstimationBegan timeEstimationBegan], [-1e10 1e10], 'k--')
plot(ts_missile, -threeSigmas(1,:), 'k-')
plot([timeEstimationEnded timeEstimationEnded], [-1e10 1e10], 'k--')
ylabel('ECI x Position Error (m)')
xlabel('Time (s)')
legend('Error', '3\sigma Bound', 'Estimation Bounds', 'Location', 'best')
ylim_low = min((-4/3)*threeSigmas(1,:));
ylim_high = max((4/3)*threeSigmas(1,:));
ylim([ylim_low ylim_high])
grid on
box on


subplot(3,2,2)
plot(ts_missile, xs_missile(4,:) - xs_estimated(4,:))
hold on
plot(ts_missile, threeSigmas(4,:), 'k-')
plot([timeEstimationBegan timeEstimationBegan], [-1e10 1e10], 'k--')
plot(ts_missile, -threeSigmas(4,:), 'k-')
plot([timeEstimationEnded timeEstimationEnded], [-1e10 1e10], 'k--')
ylabel('ECI x Velocity Error (m/s)')
xlabel('Time (s)')
legend('Error', '3\sigma Bound', 'Estimation Bounds', 'Location', 'best')
ylim_low = min((-4/3)*threeSigmas(4,:));
ylim_high = max((4/3)*threeSigmas(4,:));
ylim([ylim_low ylim_high])
grid on
box on

subplot(3,2,3)
plot(ts_missile, xs_missile(2,:) - xs_estimated(2,:))
hold on
plot(ts_missile, threeSigmas(2,:), 'k-')
plot([timeEstimationBegan timeEstimationBegan], [-1e10 1e10], 'k--')
plot(ts_missile, -threeSigmas(2,:), 'k-')
plot([timeEstimationEnded timeEstimationEnded], [-1e10 1e10], 'k--')
ylabel('ECI y Position Error (m/s)')
xlabel('Time (s)')
legend('Error', '3\sigma Bound', 'Estimation Bounds', 'Location', 'best')
ylim_low = min((-4/3)*threeSigmas(2,:));
ylim_high = max((4/3)*threeSigmas(2,:));
ylim([ylim_low ylim_high])
grid on
box on

subplot(3,2,4)
plot(ts_missile, xs_missile(5,:) - xs_estimated(5,:))
hold on
plot(ts_missile, threeSigmas(5,:), 'k-')
plot([timeEstimationBegan timeEstimationBegan], [-1e10 1e10], 'k--')
plot(ts_missile, -threeSigmas(5,:), 'k-')
plot([timeEstimationEnded timeEstimationEnded], [-1e10 1e10], 'k--')
ylabel('ECI y Velocity Error (m/s)')
xlabel('Time (s)')
legend('Error', '3\sigma Bound', 'Estimation Bounds', 'Location', 'best')
ylim_low = min((-4/3)*threeSigmas(5,:));
ylim_high = max((4/3)*threeSigmas(5,:));
ylim([ylim_low ylim_high])
grid on
box on

subplot(3,2,5)
plot(ts_missile, xs_missile(3,:) - xs_estimated(3,:))
hold on
plot(ts_missile, threeSigmas(3,:), 'k-')
plot([timeEstimationBegan timeEstimationBegan], [-1e10 1e10], 'k--')
plot(ts_missile, -threeSigmas(3,:), 'k-')
plot([timeEstimationEnded timeEstimationEnded], [-1e10 1e10], 'k--')
ylabel('ECI z Position Error (m)')
xlabel('Time (s)')
legend('Error', '3\sigma Bound', 'Estimation Bounds', 'Location', 'best')
ylim_low = min((-4/3)*threeSigmas(3,:));
ylim_high = max((4/3)*threeSigmas(3,:));
ylim([ylim_low ylim_high])
grid on
box on

subplot(3,2,6)
plot(ts_missile, xs_missile(6,:) - xs_estimated(6,:))
hold on
plot(ts_missile, threeSigmas(6,:), 'k-')
plot([timeEstimationBegan timeEstimationBegan], [-1e10 1e10], 'k--')
plot(ts_missile, -threeSigmas(6,:), 'k-')
plot([timeEstimationEnded timeEstimationEnded], [-1e10 1e10], 'k--')
ylabel('ECI z Velocity Error (m)')
xlabel('Time (s)')
legend('Error', '3\sigma Bound', 'Estimation Bounds', 'Location', 'best')
ylim_low = min((-4/3)*threeSigmas(6,:));
ylim_high = max((4/3)*threeSigmas(6,:));
ylim([ylim_low ylim_high])
grid on
box on

sgtitle('Estimation Errors with Three-Sigma Bounds')

[missDist, crossrangeError, CEP] = getMissStatistics(xs_estimated, Ps_estimated, xs_missile, ts_missile);
threeSigmasEstimation = threeSigmas(:, timeEstimationBegan+1:timeEstimationEnded);

% Compute impact in LLA (actual, estimated)
actualImpactLLA = ECI_2_LLA(xs_missile(1:3,end), ts_missile(end));

for i = timeEstimationEnded:dt:ts_missile(end)
    if norm(xs_estimated(1:3, i)) - RE < 0
        break
    end
end
estimatedImpactLLA = ECI_2_LLA(xs_estimated(1:3, i), ts_missile(i));

latError = actualImpactLLA(1) - estimatedImpactLLA(1);
lonError = actualImpactLLA(2) - estimatedImpactLLA(2);

estErrorSigmas = std(threeSigmasEstimation, 0, 2);
estErrorSigma_px = estErrorSigmas(1);
estErrorSigma_py = estErrorSigmas(2);
estErrorSigma_pz = estErrorSigmas(3);
estErrorSigma_vx = estErrorSigmas(4);
estErrorSigma_vy = estErrorSigmas(5);
estErrorSigma_vz = estErrorSigmas(6);

% Printing results

fprintf('\n===== Impact Estimation Performance Metrics =====\n');

fprintf('ECI Impact Error Norm      : %10.3f km\n', norm(missDist)/1000);
fprintf('Crossrange Error           : %10.3f km\n', crossrangeError/1000);
fprintf('Latitude Error             : %10.6f deg\n', latError * 180/pi);
fprintf('Longitude Error            : %10.6f deg\n', lonError * 180/pi);
fprintf('Circular Error Probable    : %10.3f km\n', CEP/1000);

fprintf('\n===== Estimation Uncertainty =====\n');

fprintf('Stddev Position X          : %10.3f km\n', estErrorSigma_px/1000);
fprintf('Stddev Position Y          : %10.3f km\n', estErrorSigma_py/1000);
fprintf('Stddev Position Z          : %10.3f km\n', estErrorSigma_pz/1000);
fprintf('Stddev Velocity X          : %10.3f km/s\n', estErrorSigma_vx/1000);
fprintf('Stddev Velocity Y          : %10.3f km/s\n', estErrorSigma_vy/1000);
fprintf('Stddev Velocity Z          : %10.3f km/s\n', estErrorSigma_vz/1000);

fprintf('\n====================================\n\n');

%Plot of the ground plane with predicted missile location and actual
%missile location, with CEP 
[east_error, north_error]= getCEPPlotVals(xs_missile(1:3, end), xs_estimated(1:3,i));
figure(5); clf
plot(0,0, 'go', 'MarkerSize', 10, 'LineWidth', 2)
hold on
plot(east_error/1000, north_error/1000, 'rx', 'MarkerSize', 12, 'LineWidth', 2)

theta = linspace(0, 2*pi, 250);
x_circ = (CEP/1000)*cos(theta);
y_circ = (CEP/1000)*sin(theta);
plot(x_circ, y_circ, 'k-', 'LineWidth', 2)
xlabel('East Error (km)')
ylabel('North Error (km)')
legend('Actual Impact', 'Estimated Impact', 'CEP (50%)', 'Location', 'best')
title('Impact Prediction Miss with Circular Error Probable')
grid on
box on
axis equal

%% Functions
function [fs]= spreadSatellitesEvenly(N_constellation)
% Function to spread satellites in some orbital path around that orbital
% path evenly by adjusting their desired true anomalies

    max_f = 360;
    min_f = 0;

    range = 360 - (max_f - min_f)/N_constellation;

    fs = linspace(0, range, N_constellation);

end

function [xdot]= assumedMissileDynamics(x_missile_ECI)
    % Function showing what the satellite constellation assumes the missile
    % dynamics are - this is intentionally mismodeled to reflect a
    % real-world scenario, where one would not know the dynamics of their
    % adversary's vehicle

    mu = 3.986e14; % gravitational parameter of Earth, m^3/s^2

    xdot = zeros(length(x_missile_ECI),1);

    r = x_missile_ECI(1:3);
    v = x_missile_ECI(4:6);

    xdot(1:3) = v;
    xdot(4:6) = -mu/(norm(r)^3)*r;

end

function [A]= computeA(x_missile_ECI)
    % Function to compute the state partials matrix, A = df/dx, for
    % later formation of the state transition matrix
    mu = 3.986e14; %m^3/s^2
    n = length(x_missile_ECI);
    A = zeros(n,n);

    A(1:3, 1:3) = 0;
    A(1:3, 4:6) = eye(3);
    A(4:6, 4:6) = 0;

    G_sub = [x_missile_ECI(1)^2, x_missile_ECI(1)*x_missile_ECI(2), x_missile_ECI(1)*x_missile_ECI(3);
             x_missile_ECI(1)*x_missile_ECI(2), x_missile_ECI(2)^2, x_missile_ECI(2)*x_missile_ECI(3);
             x_missile_ECI(1)*x_missile_ECI(3), x_missile_ECI(2)*x_missile_ECI(3), x_missile_ECI(3)^2];

    r = x_missile_ECI(1:3);
    G = -mu*( 1/(norm(r)^3)*eye(3) - (3/(norm(r)^5))*G_sub  ); % gravity gradient matrix

    A(4:6, 1:3) = G;

end

function [Phi]= computeSTM(A, dt, order)
    % Function to compute numerical approximation of state transition
    % matrix to order "order"
    n = length(A);
    Phi = eye(n);
    for i = 1:order
        Phi = Phi + (A^i)*(dt^i)/factorial(i);
    end

end

function [x_next]= propagateMissile(x, dt)
    % Function to use nonlinear missile dynamics model to propagate missile
    % with RK4 (more lightweight than ode45; more realistic for a satellite
    % OBC)

    k1 = assumedMissileDynamics(x);
    k2 = assumedMissileDynamics(x + (dt/2)*k1);
    k3 = assumedMissileDynamics(x + (dt/2)*k2);
    k4 = assumedMissileDynamics(x + dt*k3);

    x_next = x + (dt/6)*(k1 + 2*k2 + 2*k3 + k4);
    
end

function [yk]= formSimulatedMeasurement(SatelliteArray, measCovs, x_missile)
    % Function to package the simulated measurements from the satellites in
    % the constellation
    R_tilde = diag([measCovs.az, measCovs.el, measCovs.azrate, measCovs.elrate]);

    yk = [];

    for i = 1:length(SatelliteArray)
        if SatelliteArray(i).canObserve
            yk_tilde = SatelliteArray(i).getMeasurement(x_missile, R_tilde);
            yk = [yk; yk_tilde]; % Combine all individual satellite meas into one big meas vector
        end
    end

end

function [ykhat]= formPredictedMeasurement(SatelliteArray, x_prior)
    % Function to package the predicted measurement based on the estimated
    % state of the missile at the current iteration
    ykhat = [];
    R_tilde = zeros(4,4);
    for i = 1:length(SatelliteArray)
        if SatelliteArray(i).canObserve
            ykhat_tilde = SatelliteArray(i).getMeasurement(x_prior, R_tilde);
            ykhat = [ykhat; ykhat_tilde]; % Combine all individual predicted satellite meas into one big meas vector
        end
    end

end

function [R]= formR(measCovs, NObserving, m)
    % Function to form an overall R matrix corresponding to the whole
    % system
    R_tilde = diag([measCovs.az, measCovs.el, measCovs.azrate, measCovs.elrate]); % R for one satellite
    R = zeros(NObserving*m, NObserving*m); % overall R matrix

    for i = 1:NObserving
        if i == 1
            R = R_tilde;
        else
            R = blkdiag(R, R_tilde);
        end
    end

end

function [H]= formH(SatelliteArray, x_prior_missile)
    % Function to form an overall H matrix corresponding to all the
    % satellites in the constellation(s)
    H = [];

    for i = 1:length(SatelliteArray)
        if SatelliteArray(i).canObserve
            x_satellite = SatelliteArray(i).x;
            H_tilde = SatelliteArray(i).evaluateH(x_prior_missile - x_satellite); % uses relative missile pos for H
            H = [H; H_tilde];
        end
    end

end

function [x0, P0]= initializeEKF(yk, SatelliteArray)
    % Function to initialize the state and covariance of the filter after
    % an initial measurement from 1 or more satellites

    observingIdx = [];
    for i = 1:length(SatelliteArray)
        if SatelliteArray(i).canObserve
            observingIdx(end+1) = i;
        end
    end
    NObs = length(observingIdx);

    LOS_vecs = zeros(3,NObs); % stores the line of sight from observing sat to missile
    A_sum = zeros(3,3); % for least squares
    b_sum = zeros(3,1); % for least squares

    for k = 1:NObs
        az = yk(1 + (k-1)*4);
        el = yk(2 + (k-1)*4);

        LOS = [cos(az)*cos(el);
               sin(az)*cos(el);
               -sin(el)]; % LOS vector for current sat

        LOS_vecs(:,k) = LOS;

        A_k = eye(3) - LOS*LOS';
        i = observingIdx(k);
        r_Sat = SatelliteArray(i).x(1:3);

        A_sum = A_sum + A_k;
        b_sum = b_sum + A_k*r_Sat;

    end

    % Estimate pos
    if NObs >= 2 % triangulate
        r0 = A_sum\b_sum;
    else % N = 1, can't triangulate
        r_Sat = SatelliteArray(i).x(1:3);
        az = yk(1);
        el = yk(2);
        LOS = [cos(az)*cos(el);
               sin(az)*cos(el);
               -sin(el)];

        alt_Sat = norm(r_Sat) - 6371000;
        guess_minObservableMissileHeight = 100*1000;
        guess_range_nominal = alt_Sat - guess_minObservableMissileHeight;
        r0 = r_Sat + guess_range_nominal*LOS;

    end

    % Estimate velocity (may make more robust in future)
    v0 = 5000*(r0/norm(r0)); % 5000 m/s in the direction from Earth to guessed missile location

    % Covariance
    sigma_perpendicular = 25*1000; % variance in position perpendicular to missile LOS, m
    sigma_parallel = 250*1000; % variance in position parallel to missile LOS, m

    LOS_avg = mean(LOS_vecs, 2);
    LOS_avg = LOS_avg/norm(LOS_avg);

    P_pos = (sigma_perpendicular^2)*(eye(3) - LOS_avg*LOS_avg') + (sigma_parallel^2)*(LOS_avg*LOS_avg');
    
    sigma_vel = 2*1000; % velocity covariance in any direction, m/s
    P_vel = (sigma_vel^2)*eye(3);

    P0 = blkdiag(P_pos, P_vel);

    x0 = [r0; v0];

end

function [threeSigmas]= processCovarianceToSigmas(Ps_estimated, ts_missile)

    threeSigmas = zeros(6, length(ts_missile));
    for i = 1:length(ts_missile)
        threeSigmas(1,i) = 3*sqrt(Ps_estimated(i,1,1));
        threeSigmas(2,i) = 3*sqrt(Ps_estimated(i,2,2));
        threeSigmas(3,i) = 3*sqrt(Ps_estimated(i,3,3));
        threeSigmas(4,i) = 3*sqrt(Ps_estimated(i,4,4));
        threeSigmas(5,i) = 3*sqrt(Ps_estimated(i,5,5));
        threeSigmas(6,i) = 3*sqrt(Ps_estimated(i,6,6));
    end

end

function [missVector, crossrangeError, CEP]= getMissStatistics(xs_estimated, Ps_estimated, xs_missile, ts_missile)
    % Function to compute some desired miss statistics
    predictedTimeImpact = nan;
    for i = 1:length(ts_missile)
        if norm(xs_estimated(1:3, i)) - 6371000 < 0
            predictedImpactLocation = xs_estimated(1:3, i);
            predictedTimeImpact = ts_missile(i);
        end
    end
    if isnan(predictedTimeImpact) % state estimate results in estimated impact happening after actual impact
        predictedImpactLocation = (xs_estimated(1:3,end)/norm(xs_estimated(1:3,end)))*6371000;
        hFinal = norm(xs_estimated(1:3,i)) - 6371000; % height above ground at last propagation pt
        tFall = sqrt(2*hFinal/9.81); % how much time to fall that height
        predictedTimeImpact = ts_missile(end) + tFall;
    end

    actualImpactLocation = xs_missile(1:3, end);
    missVector = actualImpactLocation - predictedImpactLocation; % error vector

    uhat = xs_missile(1:3, end)/norm(xs_missile(1:3, end)); % direction from Earth to actual impact location
    e_radial = uhat'*missVector; % radial error 
    crossrangeError = norm(missVector - e_radial*uhat); % crossrange error
    try
        impactP = squeeze(Ps_estimated(predictedTimeImpact, :, :));
    catch
        impactP = squeeze(Ps_estimated(end, :, :)) + tFall*25*eye(6);
    end
    impactP_pos = impactP(1:3, 1:3); % covariance of position at predicted impact location

    projMatrix = eye(3) - uhat*uhat'; % projection matrix into ground plane

    impactP_groundplane = projMatrix*impactP_pos*projMatrix'; % covariance in ground plane

    eigs_impactP = sort(eig(impactP_groundplane));

    CEP = 1.177*sqrt((eigs_impactP(1) + eigs_impactP(2))/2); % circular error probability, the radius of a circle
    % inside which there is a 50% probability of the true impact being
    % there

end

function [LLA]= ECI_2_LLA(r_ECI, t)
    % Function to compute LLA coordinates given some ECI coordinates and a
    % time since ECI-ECEF frame alignment

    wE = 7.2921e-5; % rad/s
    RE = 6371000; % m

    R_ECI_2_ECEF = [cos(wE*t), sin(wE*t), 0;
                    -sin(wE*t), cos(wE*t), 0;
                    0, 0, 1];

    r_ECEF = R_ECI_2_ECEF*r_ECI;

    lon = atan2(r_ECEF(2), r_ECEF(1));
    rho = sqrt(r_ECEF(1)^2 + r_ECEF(2)^2);
    lat = atan2(r_ECEF(3), rho);

    alt = norm(r_ECEF) - RE;

    LLA = [lat; lon; alt];

end

function [east_error, north_error]= getCEPPlotVals(rImpactTrue, rImpactEst)
    % Function to calculate the east and north errors of the impact
    % location state estimate

    uhat = rImpactTrue/norm(rImpactTrue); % up direction

    zhat = [0;0;1];
    e_dir = cross(zhat, uhat);
    ehat = e_dir/norm(e_dir); % east unit vector
    nhat = cross(uhat, ehat); % north unit vector

    R_ENU_2_ECI = [ehat nhat uhat]; % rotation from ENU to ECI given unit vectors
    dr = rImpactEst - rImpactTrue;
    ENU_dr = R_ENU_2_ECI'*dr;
    east_error = ENU_dr(1);
    north_error = ENU_dr(2);


end

function [xdot]= assumedMissileDynamicsReentry(x_missile_ECI)
    % Function showing what the satellite constellation assumes the missile
    % dynamics are - this is intentionally mismodeled to reflect a
    % real-world scenario, where one would not know the dynamics of their
    % adversary's vehicle. This function applies when the satellites can no
    % longer "see" the missile due to it reentering the lower atmosphere

    mu = 3.986e14; % gravitational parameter of Earth, m^3/s^2

    xdot = zeros(length(x_missile_ECI),1);

    r = x_missile_ECI(1:3);
    v = x_missile_ECI(4:6);

    xdot(1:3) = v;

    a_grav = -mu/(norm(r)^3)*r;

    RE = 6371000;
    h = norm(r) - RE; % m
    rho0 = 1.225; % kg/m^3
    H = 7000; % scale height, m (doesn't exactly match missile dyn scale height)
    rho = rho0*exp(-h/H);
    wE = 7.2921e-5; % rad/s
    v_rel = v + cross([0;0;wE], r); % relative velocity to Earth's atmosphere
    CDA_div_M_guess = 1/1000; % general approximation of CD*A/m

    a_drag = (-0.5)*CDA_div_M_guess*rho*norm(v_rel)*v_rel;
    %a_drag = 0;
    xdot(4:6) = a_grav + a_drag;

end

function [x_next]= propagateMissileReentry(x, dt)
    % Function to use nonlinear missile dynamics model to propagate missile
    % with RK4 (more lightweight than ode45; more realistic for a satellite
    % OBC) - this one is for reentry specifically

    k1 = assumedMissileDynamicsReentry(x);
    k2 = assumedMissileDynamicsReentry(x + (dt/2)*k1);
    k3 = assumedMissileDynamicsReentry(x + (dt/2)*k2);
    k4 = assumedMissileDynamicsReentry(x + dt*k3);

    x_next = x + (dt/6)*(k1 + 2*k2 + 2*k3 + k4);
    
end

function [A]= computeAReentry(x_missile_ECI)
    % Function to compute A for assumed reentry dynamics of missile

    A = zeros(length(x_missile_ECI), length(x_missile_ECI));
    h = 1e-6; % step size base

    for i = 1:length(x_missile_ECI)

        dx = zeros(length(x_missile_ECI), 1);
        eps = h*max(abs(x_missile_ECI(i)), 1); % perturbation of state
        dx(i) = eps; % perturb only the current state coordinate to look at ith column of A

        f_plus = assumedMissileDynamicsReentry(x_missile_ECI + dx);
        f_minus = assumedMissileDynamicsReentry(x_missile_ECI - dx);

        A(:, i) = (f_plus - f_minus)/(2*eps); % Set ith column of A

    end

end



