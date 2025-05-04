close all; clear; clc;

rng(1,"twister");

% Use same random seed

set(groot,'DefaultLineLineWidth',1)

simulationParameters

dir = "NoAccelControl3D";
% dir = "AccelControl3D";

%% Run Control Dynamics Integration

fprintf('\nStarting Integration with NO Sensor Error.')
fprintf("\nTime: ")

% Initial States
p0_x = a.pr_d(1);
p0_y = a.pr_d(2);
p0_z = hdeck + a.pr_d(3);
p0 = [p0_x p0_y p0_z];

p_dot0_x = a.real_vel_xI(tspan(1));
p_dot0_y = a.real_vel_yI(tspan(1));
p_dot0_z = a.real_vel_zI(tspan(1));
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

s0 = [p0 p_dot0 pr_err_accum0 pm0 pm_dot pm_ddot];

control_dynamics = @(t, state) NoError_FixedInt_3D(t, a, state);

[~, sol_control]= rk4_solver(control_dynamics, tspan, s0, dt);

fprintf('\nFinished Integration with NO Sensor Error.\n')

fig = figure;

subplot(3,1,1)
plot(t, sol_control(:,1))
hold on
plot(t,a.real_pos_xI(t))
yline(0, 'g--')
% ylim([-0.5 0.5])
title('Platform Inertial X Position vs Time')
xlabel('Time (s)')
ylabel('Position (m)')
title('Platform Inertial X Position over Time')
legend('Fixed-Step (without sensor error) Integration', 'Deck Motion')

subplot(3,1,2)
plot(t, sol_control(:,2))
hold on
plot(t,a.real_pos_yI(t))
hold on
yline(0, 'g--')
% ylim([-0.5 0.5])
title('Platform Inertial Y Position vs Time')
xlabel('Time (s)')
ylabel('Position (m)')
title('Platform Inertial Y Position over Time')
legend('Fixed-Step (without sensor error) Integration', 'Deck Motion')

subplot(3,1,3)
plot(t, sol_control(:,3))
hold on
plot(t,a.real_pos_zI(t))
hold on
plot(t, a.real_pos_zI(t)+1)
hold on
plot(t, a.real_pos_zI(t)+0.09, '--')
hold on
plot(t, a.real_pos_zI(t)+0.5+0.41, '--')
hold on
yline(1.5, 'g--')
ylim([0.5 2.5])
title('Platform Inertial Z Position vs Time')
xlabel('Time (s)')
ylabel('Position (m)')
title('Platform Inertial Z Position over Time')
legend('Fixed-Step (without sensor error) Integration', 'Deck Motion')

saveas(fig, "figures/" + dir + "/inertial_pos_no_error.png")

%% Running Control Law Simulation WITH Sensor Error

fprintf('\nStarting Integration WITH Sensor Error.')
fprintf("\nTime: ")

% Running simulation with sensor error

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
% fig = figure;
% plot(t, sol_error(:, 1))
% xlabel('Time (sec)')
% ylabel('Velocity (m/s)')
% title('Platform Inertial Velocity')
% 
% fig = figure;
% plot(t, sol_error(:, 2))
% xlabel('Time (sec)')
% ylabel('Angle (rad)')
% title('Platform Angle')

%% Plotting Platform Position

fig = figure;

subplot(3,1,1)
plot(t, sol_control(:,1))
hold on
plot(t, sol_error(:,1))
hold on
plot(t,a.real_pos_xI(t))
title('Platform Inertial X Position vs Time')
xlabel('Time (s)')
ylabel('Position (m)')
title('Platform Inertial X Position over Time')
legend('Fixed-Step (without sensor error) Integration', 'Fixed-Step (with sensor error) Integration')

subplot(3,1,2)
plot(t, sol_control(:,2))
hold on
plot(t, sol_error(:,2))
hold on
plot(t,a.real_pos_yI(t))
title('Platform Inertial Y Position vs Time')
xlabel('Time (s)')
ylabel('Position (m)')
title('Platform Inertial Y Position over Time')
legend('Fixed-Step (without sensor error) Integration', 'Fixed-Step (with sensor error) Integration')

subplot(3,1,3)
plot(t, sol_control(:,3))
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
ylim([0.5 2.5])
title('Platform Inertial Z Position vs Time')
xlabel('Time (s)')
ylabel('Position (m)')
title('Platform Inertial Z Position over Time')
legend('Fixed-Step (without sensor error) Integration', 'Fixed-Step (with sensor error) Integration')

saveas(fig, "figures/" + dir + "/inertial_pos_with_error.png")

%% Plotting Platform Angle

fig = figure;

subplot(3,1,1)
plot(t, 180/pi*a.theta(t))
hold on
plot(t, 180/pi*sol_error(:,19))
title('Platform Angle (Theta) vs Time')
xlabel('Time (s)')
ylabel('Angle (deg)')
legend('Fixed-Step (without sensor error) Integration', 'Fixed-Step (with sensor error) Integration')

subplot(3,1,2)
plot(t, 180/pi*a.phi(t))
hold on
plot(t, 180/pi*sol_error(:,20))
title('Platform Angle (Phi) vs Time')
xlabel('Time (s)')
ylabel('Angle (deg)')
legend('Fixed-Step (without sensor error) Integration', 'Fixed-Step (with sensor error) Integration')

subplot(3,1,3)
plot(t, 180/pi*a.psi(t))
hold on
plot(t, 180/pi*sol_error(:,21))
title('Platform Angle (Psi) vs Time')
xlabel('Time (s)')
ylabel('Angle (deg)')
legend('Fixed-Step (without sensor error) Integration', 'Fixed-Step (with sensor error) Integration')

saveas(fig, "figures/" + dir + "/platform_angle.png")

%% Plotting Platform Acceleration

fig = figure;

subplot(3,1,1)
plot(t, sol_control(:, 16))
hold on
plot(t, sol_error(:,16))
% ylim([-0.01 0.01])
title('Platform Inertial Acceleration vs Time')
xlabel('Time (s)')
ylabel('Acceleration (m/s^2)')
title('Platform Inertial (X) Acceleration over Time')
legend('Deck Disturbance', 'Corrected Acceleration')

subplot(3,1,2)
plot(t, sol_control(:, 17))
hold on
plot(t, sol_error(:,17))
% ylim([-0.01 0.01])
title('Platform Inertial Acceleration vs Time')
xlabel('Time (s)')
ylabel('Acceleration (m/s^2)')
title('Platform Inertial (Y) Acceleration over Time')
legend('Deck Disturbance', 'Corrected Acceleration')

subplot(3,1,3)
plot(t, sol_control(:, 18))
hold on
plot(t, sol_error(:,18))
% ylim([-0.01 0.01])
title('Platform Inertial Acceleration vs Time')
xlabel('Time (s)')
ylabel('Acceleration (m/s^2)')
title('Platform Inertial (Z) Acceleration over Time')
legend('Deck Disturbance', 'Corrected Acceleration')

saveas(fig, "figures/" + dir + "/platform_accel.png")


%% Getting and Plotting error

fig = figure;

sz = 2;

plat_pos_x = sol_error(:,1);
plat_pos_y = sol_error(:,2);

pos_err_x = plat_pos_x - sol_control(:, 1);
pos_err_y = plat_pos_y - sol_control(:, 2);
pos_err_z = plat_pos - sol_control(:, 3);

subplot(3,1,1)
scatter(t, pos_err_x*100, sz, 'filled', displayName="Positional Error")
title('Worst Case Relative Position Error (X) vs Time')
xlabel('Time (s)')
ylabel('Error (cm)')
legend

subplot(3,1,2)
scatter(t, pos_err_y*100, sz, 'filled', displayName="Positional Error")
title('Worst Case Relative Position Error (Y) vs Time')
xlabel('Time (s)')
ylabel('Error (cm)')
legend

subplot(3,1,3)
scatter(t, pos_err_z*100, sz, 'filled', displayName="Positional Error")
title('Worst Case Relative Position Error (Z) vs Time')
xlabel('Time (s)')
ylabel('Error (cm)')
legend

saveas(fig, "figures/" + dir + "/position_error.png")
