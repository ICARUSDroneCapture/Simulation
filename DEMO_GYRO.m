close all; clear; clc;

rng(1,"twister");

% Use same random seed

set(groot,'DefaultLineLineWidth',1)

simulationParameters

%% Get Error Signal

fprintf('\nIntegrating Angular Rate (no Sensor Error)')
fprintf("\nTime: ")

% Initial Angle
phi_base = a.phi(t);
p_phi0 = phi_base(1); % Platform inertial angle [deg]

% Running simulation with sensor error

phi_over_time = a.phi(t);
phi_dot_over_time = a.phi_dot(t);
gyro_with_error_signal = a.measured_gyro_3D(a, t, phi_dot_over_time);

% Basic gyro signal (real angular rate with error)
figure
plot(t, phi_dot_over_time)
hold on
plot(t, gyro_with_error_signal)
% hold on
xlabel('Time (sec)')
ylabel('Angular Rate (rad/s)')
legend('Real Value', 'Gyro Measured Signal')
title('Angular Rate Real vs Measured')

control_dynamics = @(t_i, state) a.measured_gyro_3D(a, t_i, a.phi_dot(t_i));

[~, phi_error]= rk4_solver(control_dynamics, tspan, p_phi0, dt);

fprintf('\nFinished Integration WITH Sensor Error.\n')

%% Plotting Platform Angle

figure
plot(t, 180/pi*phi_base)
hold on
plot(t, 180/pi*phi_error)
xlabel('Time (sec)')
ylabel('Angle (deg)')
legend('NO Sensor Error', 'WITH Sensor Error')
title('Basic Integrated Angle Plot')
ylim([-30 30])

%% Getting and Plotting Angle Error

ang_err = phi_error(:, 1) - phi_base(:, 1);

sz = 2;
plot_scale = 10;

figure
scatter(t, ang_err*180/pi, sz, 'filled', displayName="Angle Error")
title('Integrated Angle Error vs Time')
xlabel('Time (s)')
ylabel('Error (deg)')

%% Control Law Drift Compensation

fprintf('\nStarting Integration WITH Sensor Compensation Control Law.')
fprintf("\nTime: ")

control_dynamics = @(t_i, state) GyroDriftCorrection(t_i, a, state);

angle_0 = 0;
angle_dot_0 = 0;

s0 = [angle_0 angle_0 angle_0 angle_dot_0 angle_dot_0 angle_dot_0];

[~, angle_corr]= rk4_solver(control_dynamics, tspan, s0, dt);

fprintf('\nFinished Integration WITH Sensor Compensation Control Law.\n')

angle_der = diff(angle_corr(:, 2)) / dt;
figure
plot(t(2:end), angle_der)
hold on
plot(t, a.phi_dot(t))
xlabel('Time (sec)')
ylabel('Angular Rate (rad/s)')
title('Gyroscope Angular Velocity')
legend('Controlled Signal', 'Real Deck Motion')

figure
plot(t, 180/pi*phi_base)
hold on
plot(t, 180/pi*phi_error)
hold on
plot(t, 180/pi*angle_corr(:, 2))
xlabel('Time (sec)')
ylabel('Angle (deg)')
legend('NO Sensor Error', 'WITH Sensor Error', 'Corrected Sensor Error')
title('Basic Integrated Angle Plot')
ylim([-30 30])


theta_dot_m_controlled_x = diff(angle_corr(:, 1))/dt;
phi_dot_m_controlled_y = diff(angle_corr(:, 2))/dt;
psi_dot_m_controlled_z = diff(angle_corr(:, 3))/dt;

figure
subplot(3,1,1)
plot(t(2:end), theta_dot_m_controlled_x/pi*180)
hold on
plot(t, a.theta_dot(t)/pi*180)
xlabel('Time (sec)')
ylabel('Angular Velocity (deg/s)')
legend('NO Sensor Error', 'Controlled Error Integrated')
title('Controlled Measured Angular Velocity Theta')
% ylim([-5 5])

subplot(3,1,2)
plot(t(2:end), phi_dot_m_controlled_y/pi*180)
hold on
plot(t, a.phi_dot(t)/pi*180)
xlabel('Time (sec)')
ylabel('Angular Velocity (deg/s)')
legend('NO Sensor Error', 'Controlled Error Integrated')
title('Controlled Measured Angular Velocity Phi')
% ylim([-5 5])

subplot(3,1,3)
plot(t(2:end), psi_dot_m_controlled_z/pi*180)
hold on
plot(t, a.psi_dot(t)/pi*180)
xlabel('Time (sec)')
ylabel('Angular Velocity (deg/s)')
legend('NO Sensor Error', 'Controlled Error Integrated')
title('Controlled Measured Angular Velocity Psi')
% ylim([-5 5])


function s_dot = GyroDriftCorrection(time_i, a, prev_state)

    % Current states
    theta = prev_state(1);
    phi = prev_state(2);
    psi = prev_state(3);
    theta_err_accum = prev_state(4);
    phi_err_accum = prev_state(5);
    psi_err_accum = prev_state(6);

    specs = a.specs;

    angle = [theta; phi; psi];
    angle_err_accum = [theta_err_accum; phi_err_accum; psi_err_accum];

    kw = a.kw(4:6);
    kt = a.kt(4:6);

    angle_0 = 0;
    angle_dot_0 = 0;

    % For this script, we are focusing on the only non-zero angle, phi
    ang_rate_i = [a.theta_dot(time_i); a.phi_dot(time_i); a.psi_dot(time_i)];
    measured_g = a.measured_gyro_3D(a, time_i, ang_rate_i);
   
    measuredState = [0 0 0 measured_g'];
    corrected_state = compensateError(measuredState, specs, time_i);

    angle_dot_m = [corrected_state(4); corrected_state(5); corrected_state(6)];

    angle_control = angle - angle_0;

    angle_dot_comp = kt .* angle_err_accum + kw .* angle_control;

    angle_dot = angle_dot_m + angle_dot_0 - angle_dot_comp;

    s_dot = zeros(6,1);

    s_dot(1:3) = angle_dot;
    s_dot(4:6) = angle_control;

end
