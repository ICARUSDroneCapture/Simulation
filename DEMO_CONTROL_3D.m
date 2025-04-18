close all; clear; clc;

rng(1,"twister");

% Use same random seed

set(groot,'DefaultLineLineWidth',1)

simulationParameters

%% Run Control Dynamics Integration

fprintf('\nStarting Integration with NO Sensor Error.')
fprintf("\nTime: ")

% Initial States

% Initial States
p0 =  hdeck + a.pr_d;             % Platform position [m]
p_dot0 = 0;             % Platform velocity [m/s]
pr_err_accum0 = 0;                      % Integral of relative position error [m*s]
pm0 = p0;                               % Platform integrated position [m]
pm_dot = p_dot0;                        % Platform integrated velocity [m/s]
pm_ddot = 0;           % Platform measured acceleration [m*s^-2]
p_theta0 = a.psi(tspan(1));   % Platform inertial angle [deg]
p_theta_err_accum0 = 0;

s0 = [p0 p_dot0 pr_err_accum0 pm0 pm_dot pm_ddot];

control_dynamics = @(t, state) NoError_FixedInt(t, a, state);

[~, sol_control]= rk4_solver(control_dynamics, tspan, s0, dt);


fprintf('\nFinished Integration with NO Sensor Error.\n')

%% Running Control Law Simulation WITH Sensor Error

fprintf('\nStarting Integration WITH Sensor Error.')
fprintf("\nTime: ")

% Running simulation with sensor error

p0_x = 0;
p0_y = 0;
p0_z = hdeck + a.pr_d;
p0 = [p0_x p0_y p0_z];

p_dot0_x = 0;
p_dot0_y = 0;
p_dot0_z = 0;
p_dot0 = [p_dot0_x p_dot0_y p_dot0_z];

pr_err_accum0_x = 0;
pr_err_accum0_y = 0;
pr_err_accum0_z = 0;
pr_err_accum0 = [pr_err_accum0_x pr_err_accum0_y pr_err_accum0_z];

pm0 = p0;

pm_dot = p_dot0;

pm_ddot_x = 0;
pm_ddot_y = 0;
pm_ddot_z = 0;
pm_ddot = [pm_ddot_x pm_ddot_y pm_ddot_z];

p_theta0_x = a.theta(tspan(1));
p_theta0_y = a.phi(tspan(1));
p_theta0_z = a.psi(tspan(1));
p_theta0 = [p_theta0_x p_theta0_y p_theta0_z];

p_theta_err_accum0_x = 0;
p_theta_err_accum0_y = 0;
p_theta_err_accum0_z = 0;
p_theta_err_accum0 = [p_theta_err_accum0_x p_theta_err_accum0_y p_theta_err_accum0_z];

s0 = [p0 p_dot0 pr_err_accum0 pm0 pm_dot pm_ddot p_theta0 p_theta_err_accum0];
control_dynamics_err = @(t, state) rigidArmControl_3D(t, a, state);

[~, sol_error]= rk4_solver(control_dynamics_err, tspan, s0, dt);

% Get platform position in inertial frame, with deck as reference zero
plat_pos = sol_error(:,3);

fprintf('\nFinished Integration WITH Sensor Error.\n')

% %% Plotting other states
% 
% figure
% plot(t, sol_error(:, 1))
% xlabel('Time (sec)')
% ylabel('Velocity (m/s)')
% title('Platform Inertial Velocity')
% 
% figure
% plot(t, sol_error(:, 2))
% xlabel('Time (sec)')
% ylabel('Angle (rad)')
% title('Platform Angle')

%% Plotting Platform Position

figure;
plot(t, sol_control(:,1))
hold on
plot(t, plat_pos)
hold on
plot(t,a.real_pos_zI(t))
hold on
plot(t, a.real_pos_zI(t)+1)
hold on
plot(t, a.real_pos_zI(t)+0.09, '--')
hold on
plot(t, a.real_pos_zI(t)+0.5+0.41, '--')
ylim([0 3])
title('Platform Inertial Position vs Time')
xlabel('Time (s)')
ylabel('Position (m)')
title('Platform Inertial Position over Time')
legend('Fixed-Step (without sensor error) Integration', 'Fixed-Step (with sensor error) Integration')


%% Plotting Platform Angle

figure;
plot(t, 180/pi*a.phi(t))
hold on
plot(t, 180/pi*sol_error(:,20))
title('Platform Angle vs Time')
xlabel('Time (s)')
ylabel('Angle (deg)')
legend('Fixed-Step (without sensor error) Integration', 'Fixed-Step (with sensor error) Integration')

%% Plotting Platform Acceleration

figure;
plot(t, sol_control(:, 6))
hold on
plot(t, sol_error(:,18))
% ylim([-0.01 0.01])
title('Platform Inertial Acceleration vs Time')
xlabel('Time (s)')
ylabel('Acceleration (m/s^2)')
title('Platform Inertial Acceleration over Time')
legend('Deck Disturbance', 'Corrected Acceleration')

%% Getting and Plotting error

pos_err = plat_pos - sol_control(:, 1);

sz = 2;
plot_scale = 0.001;

figure
scatter(t, pos_err*100, sz, 'filled', displayName="Positional Error")
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
% scatter(t(2:end), growth*100, sz, 'filled', displayName="Positional Error")
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
