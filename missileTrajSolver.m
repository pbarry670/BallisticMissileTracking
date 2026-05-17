%% Ballistic Missile Trajectory Generator
% AE 6505 Final Project - Patrick Barry

%% Constants

params.boostTime = 5*60; % 5 minute boost phase
params.ballisticTime = 30*60; % 30 minute midcourse + reentry phase
params.trajTime = params.boostTime + params.ballisticTime;

params.RE = 6371000; % Earth radius, m
params.wE = 7.2921e-5; % Earth angular rate, rad/s
params.mu = 3.986e14; % Earth gravitational parameter, m^3/s^2

params.rho0 = 1.225; % sea level rho, kg/m^3
params.H = 8500; % scale height, m (for atmospheric model)

params.CD = 2.2; % coefficient of drag
params.A = 1; % effective area, m^2
params.mDry = 2000; % dry mass, kg
params.mWet = 10000; % wet mass, kg
params.mFuel = params.mWet - params.mDry; % fuel mass, kg
params.Isp = 290; % specific impulse of given solid fuel, s
params.g0 = 9.81; % sea level gravity, m/s^2
params.Tmax = 150000; % max thrust, N

params.dt = 1;


%% Coords of Launch Site, Target

launchSite_LLA = [38.890262*(pi/180); % lat, rads
                  -77.036271*(pi/180); % lon, rads
                  12]; % alt, m

target_LLA = [55.752268*(pi/180);  % lat, rads
              37.623276*(pi/180); % lon, rads
              512];  % alt, m

r_launchSite_ECEF = [(params.RE+launchSite_LLA(3))*cos(launchSite_LLA(1))*cos(launchSite_LLA(2));
                     (params.RE+launchSite_LLA(3))*cos(launchSite_LLA(1))*sin(launchSite_LLA(2));
                     (params.RE+launchSite_LLA(3))*sin(launchSite_LLA(1))]; % launch site location in ECEF frame

r_target_ECEF = [(params.RE+target_LLA(3))*cos(target_LLA(1))*cos(target_LLA(2));
                     (params.RE+target_LLA(3))*cos(target_LLA(1))*sin(target_LLA(2));
                     (params.RE+target_LLA(3))*sin(target_LLA(1))]; % target location in ECEF frame

figure(1); clf
scatter3(r_target_ECEF(1), r_target_ECEF(2), r_target_ECEF(3), 80, 'filled')
hold on
scatter3(r_launchSite_ECEF(1), r_launchSite_ECEF(2), r_launchSite_ECEF(3), 80, 'filled')
[xe, ye, ze] = sphere(50); % arg is resolution
xe = params.RE*xe;
ye = params.RE*ye;
ze = params.RE*ze;
surf(xe, ye, ze, ...
    'FaceColor', [0.2 0.6 1], ...
    'EdgeColor', 'none', ...
    'FaceAlpha', 0.5); % semi-transparent
axis equal
xlabel('ECEF x (m)')
ylabel('ECEF y (m)')
zlabel('ECEF z (m)')

R_ECEF_2_ECI_t0 = eye(3);
R_ECEF_2_ECI_tf = [cos(params.wE*params.trajTime), -sin(params.wE*params.trajTime), 0;
                   sin(params.wE*params.trajTime), cos(params.wE*params.trajTime), 0;
                    0, 0, 1]; % rotation from ECEF to ECI at final time

r_launchSite_ECI = R_ECEF_2_ECI_t0*r_launchSite_ECEF;
r_target_ECI = R_ECEF_2_ECI_tf*r_target_ECEF;

%% Part I: Trajectory solver for midcourse/coast phase from some point in the air to target impact with no burn

r0_guess = r_launchSite_ECI + [0; 0; 200000]; % 200 km altitude
v0_guess = cross([0;0;params.wE], r0_guess) + [0; 0; 7500]; % rough orbital-scale guess
tspan = [params.boostTime params.trajTime];

options = optimoptions('fsolve', ...
    'Display','iter', ...
    'TolFun',1e-6, ...
    'TolX',1e-6);

v0_solved = fsolve(@(v) computeResidual(r0_guess, v, tspan, r_target_ECI, params), ...
                v0_guess, options);

x0_ballistic = [r0_guess; v0_solved];

tspan = params.boostTime:params.dt:params.trajTime;
[ts, xs_ballistic] = ode45(@(t,x) ballisticDynamics(t,x,params), tspan, x0_ballistic);
rs_ballistic = xs_ballistic(:,1:3);

figure(2); clf
scatter3(r_launchSite_ECI(1), r_launchSite_ECI(2), r_launchSite_ECI(3), 80, 'filled')
hold on
scatter3(r_target_ECI(1), r_target_ECI(2), r_target_ECI(3), 80, 'filled')
plot3(rs_ballistic(:,1), rs_ballistic(:,2), rs_ballistic(:,3), 'LineWidth',2)
surf(xe, ye, ze, ...
    'FaceColor', [0.2 0.6 1], ...
    'EdgeColor', 'none', ...
    'FaceAlpha', 0.5); % semi-transparent
axis equal
xlabel('ECEF x (m)')
ylabel('ECEF y (m)')
zlabel('ECEF z (m)')
axis equal
grid on
xlabel('ECI x', 'FontSize', 20);
ylabel('ECI y', 'FontSize', 20);
zlabel('ECI z', 'FontSize', 20);
title('Ballistic Trajectory (ECI)', 'FontSize', 20)

disp("Final Missile Position:")
disp(xs_ballistic(end,1:3))

disp("Final Missile Velocity:")
disp(xs_ballistic(end, 4:6))

disp("End of Boost Phase Position:")
disp(xs_ballistic(1,1:3))

disp("End of Boost Phase Velocity:")
disp(xs_ballistic(1,4:6))

%% Part II. Trajectory Approximation for Boost Phase

r0_boost = r_launchSite_ECI;
v0_boost = [0;0;0] + cross([0;0;params.wE], r_launchSite_ECI);
x0_boost = [r0_boost; v0_boost];

xf_boost = xs_ballistic(1,:)';

% Boundary states
r0 = x0_boost(1:3);
v0 = x0_boost(4:6);
rf = xf_boost(1:3);
vf = xf_boost(4:6);
tf = params.boostTime;

dt_boost = params.dt;
ts = 0:dt_boost:tf-1;  % time vector

N = length(ts);

rs_boost = zeros(3,N);
vs_boost = zeros(3,N);

for i = 1:3
    % Quadratic coefficients - match initial vel
    a = (rf(i) - r0(i) - v0(i)*tf) / tf^2;
    b = v0(i);
    c = r0(i);

    % Quadratic - match final vel
    %a = (r0(i) - rf(i) + vf(i)*tf)/(tf^2);
    %b = vf(i) - 2*a*tf;
    %c = r0(i);

    % Cubic coefficients
    %a = (2*(r0(i)-rf(i)) + (v0(i)+vf(i))*tf)/(tf^3);
    %b = (3*(rf(i)-r0(i)) - (2*v0(i)+vf(i))*tf)/(tf^2);
    %c = v0(i);
    %d = r0(i);

    rs_boost(i,:) = a*ts.^2 + b*ts + c;
    s = ts/params.boostTime;
    vs_boost(i,:) = v0(i)*(1-s) + vf(i)*s;
end

xs_boost = [rs_boost; vs_boost];

figure(3); clf
plot3(xs_boost(1,:)/1000, xs_boost(2,:)/1000, xs_boost(3,:)/1000,'LineWidth',2)
hold on
scatter3(x0_boost(1)/1000,  x0_boost(2)/1000, x0_boost(3)/1000, 80, 'filled')
scatter3(xf_boost(1)/1000,xf_boost(2)/1000, xf_boost(3)/1000, 80, 'filled')
axis equal
grid on
xlabel('ECI x (km)', 'FontSize', 20)
ylabel('ECI y (km)', 'FontSize', 20)
zlabel('ECI z (km)', 'FontSize', 20)
title('Boost Phase Trajectory', 'FontSize', 20)

%% Total Trajectory

xs_ballistic = xs_ballistic';

xs_missile = [xs_boost, xs_ballistic];
ts_missile = 0:params.dt:params.trajTime;

% Add phase of flight labels and missile apparent brightness coefficients
phaseOfFlight = cell(length(ts_missile), 1);
betas = zeros(length(ts_missile), 1);
for i = 1:length(ts_missile)
    if i < params.boostTime
        phaseOfFlight(i) = {'BOOST'};
        betas(i) = 500;
    elseif i < params.boostTime + 20*60
        phaseOfFlight(i) = {'MIDCOURSE'};
        betas(i) = 100;
    else
        phaseOfFlight(i) = {'REENTRY'};
        betas(i) = 250;
    end
end

figure(4); clf
scatter3(r_launchSite_ECI(1), r_launchSite_ECI(2), r_launchSite_ECI(3), 80, 'filled')
hold on
scatter3(r_target_ECI(1), r_target_ECI(2), r_target_ECI(3), 80, 'filled')
plot3(xs_missile(1,:), xs_missile(2,:), xs_missile(3,:), 'LineWidth',2)
surf(xe, ye, ze, ...
    'FaceColor', [0.2 0.6 1], ...
    'EdgeColor', 'none', ...
    'FaceAlpha', 0.5); % semi-transparent
axis equal
grid on
xlabel('ECI x (m)', 'FontSize', 20);
ylabel('ECI y (m)', 'FontSize', 20);
zlabel('ECI z (m)', 'FontSize', 20);
title('Missile Trajectory', 'FontSize', 20)


alts = vecnorm(xs_missile(1:3,:)) - params.RE;
figure(5); clf
plot(ts_missile, alts)
xlabel('Time (s)')
ylabel('Missile Altitude (m)')
title('Altitude Plot')

disp("Max Altitude Achieved (m):")
disp(max(alts))

disp("Miss Distance (m):")
disp(xs_missile(1:3,end)' - r_target_ECI')

%% Save variables

save('missileTraj.mat', 'xs_missile', 'ts_missile', 'phaseOfFlight', 'betas', 'r_target_ECI', 'r_launchSite_ECI');


%% Functions

function [xdot]= ballisticDynamics(t,x,params) 
    % Dynamics during midcourse/coast phase

    r = x(1:3);
    v = x(4:6);

    rdot = v;
    
    a_grav = (-params.mu/(norm(r)^3))*r;
    
    h = norm(r) - params.RE;
    rho = params.rho0*exp(-h/params.H);
    v_rel = v - cross([0;0;params.wE], r);
    a_drag = (-0.5)*params.CD*(params.A/params.mDry)*rho*norm(v_rel)*v_rel;

    vdot = a_grav + a_drag;

    xdot = [rdot; vdot];

end

function [miss_vec] = computeResidual(r0, v0, tspan, r_target, params)

    x0 = [r0; v0];

    opts = odeset('RelTol',1e-9,'AbsTol',1e-9);
    [~, xs] = ode45(@(t,x) ballisticDynamics(t,x,params), tspan, x0, opts);

    r_final = xs(end,1:3)';

    miss_vec = r_final - r_target;

end

function [xdot]= boostDynamics(t,x,u,params)
    % Dynamics during boost phase
    r = x(1:3);
    v = x(4:6);
    m = x(7);

    rdot = v;
    
    a_grav = (-params.mu/(norm(r)^3))*r;
    
    h = norm(r) - params.RE;
    rho = params.rho0*exp(-h/params.H);
    v_rel = v - cross([0;0;params.wE], r);
    a_drag = (-0.5)*params.CD*(params.A/m)*rho*norm(v_rel)*v_rel;

    T = u(4)*params.Tmax;
    a_thrust = (T/m)*u(1:3);

    vdot = a_grav + a_drag + a_thrust;

    mdot = -T/(params.Isp*params.g0);

    xdot = [rdot; vdot; mdot];

end

function [cost]= boostCost(z, x0, xf_des, params)

    [~, xs_hist] = propagateBoost(z, x0, params);
    xf = xs_hist(end,:)';

    r_err = xf(1:3) - xf_des(1:3);
    v_err = xf(4:6) - xf_des(4:6);
    %m_err = xf(7) - xf_des(7); % maybe omit

    cost = [r_err/1e6;
            v_err/1e3];
            %m_err/1e3]; % maybe omit

end

function [t_hist, x_hist] = propagateBoost(z, x0, params)

    dt = params.boostTime/params.N;

    x = x0; % initial state
    t0 = 0;

    t_hist = []; % allocate for storage of traj history
    x_hist = [];

    for i = 1:params.N

        % Extract decision variables
        theta = z((i-1)*3 + 1);
        phi   = z((i-1)*3 + 2);
        throttle = z((i-1)*3 + 3);

        u = [cos(theta)*cos(phi);
             cos(theta)*sin(phi);
             sin(theta)]; % thrust vector direction in ECI
        u = [u; throttle]; % full control vector

        tspan = [t0, t0 + dt];
        [t, xs] = ode45(@(t,x) boostDynamics(t,x,u,params), tspan, x);

        t_hist = [t_hist; t];
        x_hist = [x_hist; xs];

        x = xs(end,:)'; % make column vector for ode45()
        t0 = t0 + dt; % update time

    end

end

