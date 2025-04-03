close all; clear; clc;

rng(1,"twister");

% Use same random seed

set(groot,'DefaultLineLineWidth',1)

simulationParameters;

close all;

%% Sensor Model Aspects

% Simulation time
startTime = 0;
finishTime = 20;
tspan = [startTime finishTime]; % [s]

% dt = 1/imu_rate;  % [s]
dt = 0.0001;
t = (tspan(1):dt:tspan(2))';
t_count = length(t);
indeces = @(t) floor(t/dt)+1;

defineSignals
% defineSignalsNoNoise

%% 3D motion equations

a.real_pos_xI = @(t) 0.4/(beta^2)*sin(beta*t/2);
a.real_pos_yI = @(t) 1.6/(beta^2)*sin(beta*t/4);
a.real_pos_zI = @(t) alpha*sin(beta*t) + hdeck;

a.real_vel_xI = @(t) 0.2/beta*cos(beta/2*t);
a.real_vel_yI = @(t) 0.4/beta*cos(beta/4*t);
a.real_vel_zI = @(t) beta*alpha*cos(beta*t);

a.real_accel_xI = @(t) -0.1*sin(beta/2*t); % [m*s^-2]
a.real_accel_yI = @(t) -0.1*sin(beta/4*t); % [m*s^-2]
a.real_accel_zI = @(t) -beta^2*alpha*sin(beta*t) - 9.81; % [m*s^-2]

% a.theta = @(t) -atan(beta*alpha*cos(beta*t)); % [rad]
% a.phi = @(t) atan(0.2/beta*cos(beta/2*t)); % [rad]
% a.psi = @(t) atan(0.4/beta*cos(beta/4*t)); % [rad]


a.psi = @(t) -atan(beta*alpha*cos(beta*t)); % [rad]
a.theta = @(t) atan(0.2/beta*cos(beta/2*t)); % [rad]
a.phi = @(t) atan(0.4/beta*cos(beta/4*t)); % [rad]

% a.theta_dot = @(t) ((alpha*beta^2*sin(beta*t))./(alpha^2*beta^2*(cos(beta*t).^2)+1)); % [rad/s]
% a.phi_dot = @(t) ((-0.1*sin(beta/2*t))/(((0.04*(cos(beta*t/2).^2))/(beta^2))+1)); % [rad/s]
% a.psi_dot = @(t) ((-0.1*sin(beta/4*t))/(((0.16*(cos(beta*t/4).^2))/(beta^2))+1)); % [rad/s]

% Angular rate needs to be taken manually since the derivative equation is +/-
theta_vals = a.theta(t);
phi_vals = a.phi(t);
psi_vals = a.psi(t);

theta_dot = zeros(1, length(t));
phi_dot = zeros(1, length(t));
psi_dot = zeros(1, length(t));

theta_dot(2:end) = diff(theta_vals)/dt;
phi_dot(2:end) = diff(phi_vals)/dt;
psi_dot(2:end) = diff(psi_vals)/dt;

a.theta_dot_eq = @(t) theta_dot(floor(t./dt)+1);
a.phi_dot_eq = @(t) phi_dot(floor(t./dt)+1);
a.psi_dot_eq = @(t) psi_dot(floor(t./dt)+1);

s = @(x) sin(x);
c = @(x) cos(x);

%% Simplifying to just z-direction equations for now
a.real_pos = @(t) alpha*sin(beta*t) + hdeck;
a.real_vel = @(t) beta*alpha*cos(beta*t);
a.real_accel = @(t, y) -beta^2*alpha*sin(beta*t) - 9.81; % [m*s^-2]
a.real_ang = @(t) atan(0.4/beta*cos(beta/4*t)); % [rad]
a.real_ang_rate = @(t, y) theta_dot(floor(t./dt)+1);

%% Run Control Dynamics Integration

a.dt = dt;

fprintf('\nStarting Integration with NO Sensor Error.')
fprintf("\nTime: ")

% % Initial States
% p0 =  a.d(tspan(1))+a.pr_d;   % Platform position [m]
% p_dot0 = 0;   % Platform velocity [m/s]
% pr_err_accum0 = 0;            % Integral of relative position error [m*s]
% pm0 = p0;                     % Platform inetegrated position [m]
% pm_dot = p_dot0;              % Platform integrated velocity [m/s]
% pm_ddot = a.d_ddot(tspan(1)); % Platform measured acceleration [m*s^-2]
% p_theta0 = a.real_ang(tspan(1)); % Platform inertial angle [deg]

% Initial States
p0 =  a.d(tspan(1))+a.pr_d;             % Platform position [m]
p_dot0 = 0;             % Platform velocity [m/s]
pr_err_accum0 = 0;                      % Integral of relative position error [m*s]
pm0 = p0;                               % Platform integrated position [m]
pm_dot = p_dot0;                        % Platform integrated velocity [m/s]
pm_ddot = 0;           % Platform measured acceleration [m*s^-2]
p_theta0 = a.real_ang(tspan(1));   % Platform inertial angle [deg]
p_theta_err_accum0 = 0;

s0 = [p0 p_dot0 pr_err_accum0 pm0 pm_dot pm_ddot];

control_dynamics = @(t, state) NoError_FixedInt(t, a, state);

[t_control, sol_control]= rk4_solver(control_dynamics, tspan, s0, dt);


fprintf('\nFinished Integration with NO Sensor Error.\n')

%% Running Control Law Simulation WITH Sensor Error

fprintf('\nStarting Integration WITH Sensor Error.')
fprintf("\nTime: ")

% Running simulation with sensor error

p0_x = 0;
p0_y = 0;
p0_z = 0;
p0 = [p0_x p0_y p0_z];

p_dot0_x = 0;
p_dot0_y = 0;
p_dot0_z = 0;
p_dot0 = [p_dot0_x p_dot0_y p_dot0_z];

pr_err_accum0_x = 0;
pr_err_accum0_y = 0;
pr_err_accum0_z = 0;
pr_err_accum0 = [pr_err_accum0_x pr_err_accum0_y pr_err_accum0_z];

pm0_x = 0;
pm0_y = 0;
pm0_z = 0;
pm0 = [pm0_x pm0_y pm0_z];

pm_dot_x = 0;
pm_dot_y = 0;
pm_dot_z = 0;
pm_dot = [pm_dot_x pm_dot_y pm_dot_z];

pm_ddot_x = 0;
pm_ddot_y = 0;
pm_ddot_z = 0;
pm_ddot = [pm_ddot_x pm_ddot_y pm_ddot_z];

p_theta0_x = 0;
p_theta0_y = 0;
p_theta0_z = 0;
p_theta0 = [p_theta0_x p_theta0_y p_theta0_z];

p_theta_err_accum0_x = 0;
p_theta_err_accum0_y = 0;
p_theta_err_accum0_z = 0;
p_theta_err_accum0 = [p_theta_err_accum0_x p_theta_err_accum0_y p_theta_err_accum0_z];

s0 = [p0 p_dot0 pr_err_accum0 pm0 pm_dot pm_ddot p_theta0 p_theta_err_accum0];
control_dynamics_err = @(t, state) rigidArmControl_3D(t, a, state);

[t_error, sol_error]= rk4_solver(control_dynamics_err, tspan, s0, dt);

% Get platform position in inertial frame, with deck as reference zero
plat_pos = sol_error(:,1);

fprintf('\nFinished Integration WITH Sensor Error.\n')

% %% Plotting other states
% 
% figure
% plot(t_error, sol_error(:, 1))
% xlabel('Time (sec)')
% ylabel('Velocity (m/s)')
% title('Platform Inertial Velocity')
% 
% figure
% plot(t_error, sol_error(:, 2))
% xlabel('Time (sec)')
% ylabel('Angle (rad)')
% title('Platform Angle')

%% Plotting Platform Position

figure;
plot(t_control, sol_control(:,3))
hold on
plot(t_error, plat_pos)
hold on
plot(t,a.d(t))
hold on
plot(t, a.d(t)+1)
hold on
plot(t, a.d(t)+0.09, '--')
hold on
plot(t, a.d(t)+0.5+0.41, '--')
title('Platform Inertial Position vs Time')
xlabel('Time (s)')
ylabel('Position (m)')
title('Platform Inertial Position over Time')
legend('Fixed-Step (without sensor error) Integration', 'Fixed-Step (with sensor error) Integration')


%% Plotting Platform Angle

figure;
plot(t_error, 180/pi*a.theta(t))
hold on
plot(t_error, 180/pi*sol_error(:,19))
title('Platform Angle vs Time')
xlabel('Time (s)')
ylabel('Angle (deg)')
legend('Fixed-Step (without sensor error) Integration', 'Fixed-Step (with sensor error) Integration')

%% Plotting Platform Acceleration

figure;
plot(t, a.real_accel_zI(t))
hold on
plot(t, sol_error(:,18))
% ylim([-0.01 0.01])
title('Platform Inertial Acceleration vs Time')
xlabel('Time (s)')
ylabel('Acceleration (m/s^2)')
title('Platform Inertial Acceleration over Time')
legend('Deck Disturbance', 'Corrected Acceleration')

%% Getting and Plotting error

pos_err = plat_pos - sol_control(:, 3);

sz = 2;
plot_scale = 0.001;

figure
scatter(t_error, pos_err*100, sz, 'filled', displayName="Positional Error")
% hold on
% plot(t,a.d(t)/200, displayName="Deck Disturbance")
% hold on
% plot(t,plot_scale*a.d(t))
% hold on
% plot(t, plot_scale*(a.d(t)+1))
% hold on
% plot(t, plot_scale*(a.d(t)+0.09), '--')
% hold on
% plot(t, plot_scale*(a.d(t)+0.5+0.41), '--')
title('Worst Case Relative Position Error vs Time')
xlabel('Time (s)')
ylabel('Error (cm)')
legend

plot_scale = 0.00001;
growth = diff(pos_err);

% figure
% scatter(t_error(2:end), growth*100, sz, 'filled', displayName="Positional Error")
% hold on
% plot(t,plot_scale*a.d(t))
% hold on
% plot(t, plot_scale*(a.d(t)+1))
% hold on
% plot(t, plot_scale*(a.d(t)+0.09), '--')
% hold on
% plot(t, plot_scale*(a.d(t)+0.5+0.41), '--')
% % hold on
% % plot(t,a.d(t)/200, displayName="Deck Disturbance")
% title('Error Growth over Time')
% xlabel('Time (s)')
% ylabel('Error (cm)')
% legend
