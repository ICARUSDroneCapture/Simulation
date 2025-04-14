close all; clear; clc;

rng(1,"twister");

% Use same random seed

set(groot,'DefaultLineLineWidth',1)

simulationParameters

%% Get Error Signal

fprintf('\nIntegrating Acceleration (with Sensor Error)')
fprintf("\nTime: ")

% Initial Angle
vel_base = a.real_vel_zI(t);
p_vel0 = vel_base(1); % Platform inertial angle [deg]

% Running simulation with sensor error

vel_z_over_time = a.real_vel_zI(t);
accel_z_over_time = a.real_accel_zI(t);
accel_with_error_signal = a.measured_accel_3D(a, t, accel_z_over_time);

% Basic gyro signal (real angular rate with error)
figure
plot(t, a.real_accel_zI(t))
hold on
plot(t, accel_with_error_signal)
xlabel('Time (sec)')
ylabel('Acceleration (m/s^2)')
legend('Real Value (with g)', 'Accel Measured Signal (inertial z)')
title('Acceleration Real vs Measured')

control_dynamics = @(t_i, state) a.measured_accel_3D(a, t_i, a.real_accel_zI(t_i)+a.g);

[~, vel_error]= rk4_solver(control_dynamics, tspan, p_vel0, dt);

fprintf('\nFinished Integration WITH Sensor Error.\n')

%% Plotting Platform Velocity

figure
plot(t, a.real_accel_zI(t))
hold on
plot(t, vel_error)
xlabel('Time (sec)')
ylabel('Velocity (m/s)')
legend('NO Sensor Error', 'WITH Sensor Error')
title('Basic Integrated Velocity Plot')

%% Getting and Plotting Velocity Error

vel_diff = vel_error(:, 1) - vel_base(:, 1);

sz = 2;
plot_scale = 10;

figure
scatter(t, vel_diff, sz, 'filled', displayName="Velocity Error")
title('Integrated (vertical) Velocity Error vs Time')
xlabel('Time (s)')
ylabel('Error (m/s)')

%% Control Law Drift Compensation

fprintf('\nStarting Integration WITH Sensor Compensation Control Law.')
fprintf("\nTime: ")

control_dynamics = @(t_i, state) AccelDriftCorrection(t_i, a, state);

vel_0 = 0;
accel_0 = 0;

s0 = [vel_0 vel_0 vel_0 accel_0 accel_0 accel_0];

[~, vel_corr]= rk4_solver(control_dynamics, tspan, s0, dt);

a_der = diff(vel_corr(:, 3)) / dt;
figure
plot(t(2:end), a_der)
hold on
plot(t, a.real_accel_zI(t) + a.g)
xlabel('Time (sec)')
ylabel('Acceleration (m/s^2)')
legend('Controlled Acceleration', 'Real Deck Acceleration')

fprintf('\nFinished Integration WITH Sensor Compensation Control Law.\n')

figure
plot(t, vel_base)
hold on
plot(t, vel_error(:, 1))
hold on
plot(t, vel_corr(:, 3))
xlabel('Time (sec)')
ylabel('Velocity (m/s)')
legend('NO Sensor Error', 'WITH Sensor Error', 'Corrected (z axis) Accel Signal')
title('Basic Integrated Velocity Plot')
% ylim([-5 5])


accel_m_controlled_x = diff(vel_corr(:, 1))/dt;
accel_m_controlled_y = diff(vel_corr(:, 2))/dt;
accel_m_controlled_z = diff(vel_corr(:, 3))/dt;

figure
subplot(3,1,1)
plot(t(2:end), accel_m_controlled_x)
hold on
plot(t, a.real_accel_xI(t))
xlabel('Time (sec)')
ylabel('Acceleration (m/s^2)')
legend('NO Sensor Error', 'Controlled Error Integrated')
title('Controlled Measured Acceleration X')
% ylim([-5 5])

subplot(3,1,2)
plot(t(2:end), accel_m_controlled_y)
hold on
plot(t, a.real_accel_yI(t))
xlabel('Time (sec)')
ylabel('Acceleration (m/s^2)')
legend('NO Sensor Error', 'Controlled Error Integrated')
title('Controlled Measured Acceleration Y')
% ylim([-5 5])

subplot(3,1,3)
plot(t(2:end), accel_m_controlled_z)
hold on
plot(t, a.real_accel_zI(t) + 9.81)
xlabel('Time (sec)')
ylabel('Acceleration (m/s^2)')
legend('NO Sensor Error', 'Controlled Error Integrated')
title('Controlled Measured Acceleration Z')
% ylim([-5 5])


function s_dot = AccelDriftCorrection(time_i, a, prev_state)

    % Current states
    vel_x = prev_state(1);
    vel_y = prev_state(2);
    vel_z = prev_state(3);
    vel_err_accum_x = prev_state(4);
    vel_err_accum_y = prev_state(5);
    vel_err_accum_z = prev_state(6);

    specs = a.specs;

    vel = [vel_x; vel_y; vel_z];
    vel_err_accum = [vel_err_accum_x; vel_err_accum_y; vel_err_accum_z];

    kw = a.kw(1:3);
    kt = a.kt(1:3);

    vel_0 = 0;
    vel_dot_0 = 0;

    a_I = [a.real_accel_xI(time_i); a.real_accel_yI(time_i); a.real_accel_zI(time_i)];
    
    theta_real = a.theta(time_i);
    phi_real = a.phi(time_i);
    psi_real = a.psi(time_i);

    a_S = Rotate_I_S(a_I, theta_real, phi_real, psi_real);
    
    accel_S = a.measured_accel_3D(a, time_i, a_S);

    accel_I = Rotate_S_I(accel_S, theta_real, phi_real, psi_real);

    accel_state = [accel_I(1) accel_I(2) accel_I(3)+a.g];

    measuredState = [accel_state 0 0 0];
    corrected_state = compensateError(measuredState, specs, time_i);
    
    vel_dot_m = [corrected_state(1); corrected_state(2); corrected_state(3)];

    vel_control = vel - vel_0;

    vel_dot_comp = kt .* vel_err_accum + kw .* vel_control;

    vel_dot = vel_dot_m + vel_dot_0 - vel_dot_comp;

    s_dot = zeros(6,1);

    s_dot(1:3) = vel_dot;
    s_dot(4:6) = vel_control;

end
