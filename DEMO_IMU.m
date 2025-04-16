close all; clear; clc;

rng(1,"twister");

% Use same random seed

set(groot,'DefaultLineLineWidth',1)

simulationParameters

%% Sensor Frame Accelerations

a_I_over_time = [a.real_accel_xI(t), a.real_accel_yI(t), a.real_accel_zI(t)];

theta_over_time = a.theta(t);
phi_over_time = a.phi(t);
psi_over_time = a.psi(t);

a_S_true = zeros(length(t), 3);
a_S_measured = zeros(length(t), 3);
a_I_ref = zeros(length(t), 3);

for i = 1:length(t)

    t_i = t(i);

    theta_i = theta_over_time(i);
    phi_i = phi_over_time(i);
    psi_i = psi_over_time(i);

    a_I_i = a_I_over_time(i, :)';
    
    % Convert inertial frame accelerations to sensor frame
    a_S_i = Rotate_I_S(a_I_i, theta_i, phi_i, psi_i);
    
    % Get measured sensor frame accelerations
    a_S_m = a.measured_accel_3D(a, t_i, a_S_i);
    
    % Convert Measured Accelerations back to Inertial without any Correction
    a_I_i = Rotate_S_I(a_S_m, theta_i, phi_i, psi_i);

    a_S_true(i, :) = a_S_i';
    a_S_measured(i, :) = a_S_m';
    a_I_ref(i, :) = a_I_i';

end

figure
subplot(3,1,1)
plot(t, a.real_accel_xI(t))
hold on
plot(t, a_S_measured(:, 1))
xlabel('Time (sec)')
ylabel('Acceleration (m/s^2)')
title('Measured Sensor Acceleration X (m/s^2)')
legend('Inertial Frame Acceleration', 'Sensor Frame Acceleration')

subplot(3,1,2)
plot(t, a.real_accel_yI(t))
hold on
plot(t, a_S_measured(:, 2))
xlabel('Time (sec)')
ylabel('Acceleration (m/s^2)')
title('Measured Sensor Acceleration Y (m/s^2)')
legend('Inertial Frame Acceleration', 'Sensor Frame Acceleration')

subplot(3,1,3)
plot(t, a.real_accel_zI(t))
hold on
plot(t, a_S_measured(:, 3))
xlabel('Time (sec)')
ylabel('Acceleration (m/s^2)')
title('Measured Sensor Acceleration Z (m/s^2)')
legend('Inertial Frame Acceleration', 'Sensor Frame Acceleration')


figure
subplot(3,1,1)
plot(t, a.real_accel_xI(t))
hold on
plot(t, a_I_ref(:, 1))
xlabel('Time (sec)')
ylabel('Acceleration (m/s^2)')
title('Inertial Acceleration X Converted From Sensor (m/s^2)')
legend('Inertial Frame Acceleration', 'Sensor Frame Acceleration')

subplot(3,1,2)
plot(t, a.real_accel_yI(t))
hold on
plot(t, a_I_ref(:, 3))
xlabel('Time (sec)')
ylabel('Acceleration (m/s^2)')
title('Inertial Acceleration Y Converted From Sensor (m/s^2)')
legend('Inertial Frame Acceleration', 'Sensor Frame Acceleration')

subplot(3,1,3)
plot(t, a.real_accel_zI(t))
hold on
plot(t, a_I_ref(:, 3))
xlabel('Time (sec)')
ylabel('Acceleration (m/s^2)')
title('Inertial Acceleration Z Converted From Sensor (m/s^2)')
legend('Inertial Frame Acceleration', 'Sensor Frame Acceleration')

% %% Get Error Signal
% 
% fprintf('\nStarting Integration WITH Sensor Error.')
% fprintf("\nTime: ")
% 
% % Running simulation with sensor error
% 
% error_signals = zeros(t_count, 6);
% 
% for i = 1:length(t)
%     time_i = t(i);
% 
%     ang_rate_i = a.real_ang_rate(time_i);
%     accel_i = a.real_accel(time_i);
% 
%     measured_gyro = a.measured_gyro(time_i, a.o_d_n_g_c, a.biasStabDistGyro, a.biasTempDistGyro, a.gyro_drift, a.noiseDistGyro, ang_rate_i);
%     measured_accel_v = a.measured_accel_vert(time_i, a.o_d_n_a_c_v, a.biasStabDistAccel, a.biasTempDistAccel, a.accel_drift_vert, a.noiseDistAccel, accel_i, a.theta_err);
%     measured_accel_h = a.measured_accel_horz(time_i, a.o_d_n_a_c_h, a.biasStabDistAccel, a.biasTempDistAccel, a.accel_drift_horz, a.noiseDistAccel, accel_i, a.theta_err);
% 
%     measuredState = [measured_accel_h measured_accel_h measured_accel_v measured_gyro measured_gyro measured_gyro];
% 
%     error_signals(i, :) = measuredState;
% end
% 
% theta0_x = 0;
% theta0_y = 0;
% theta0_z = 0;
% vel0_x = 0;
% vel0_y = 0;
% vel0_z = 0;
% 
% s0 = [theta0_x theta0_y theta0_z vel0_x vel0_y vel0_z];
% 
% control_dynamics = @(t_i, state) error_signals(floor(t_i./dt)+1)';
% 
% [t_error, int_signals_error]= rk4_solver(control_dynamics, tspan, s0, dt);
% 
% fprintf('\nFinished Integration WITH Sensor Error.\n')
% 
% figure
% plot(t_error, error_signals(:, 1))
% hold on
% plot(t_error, error_signals(:, 3))
% hold on
% plot(t, a.real_accel(t))
% xlabel('Time (sec)')
% ylabel('Acceleration (m/s^2)')
% legend('WITH Sensor Error (horizontal)', 'WITH Sensor Error (vertical)', 'NO Sensor Error')
% title('Measured Accelerometer Signal')
% % ylim([-5 5])
% 
% figure
% plot(t_error, error_signals(:, 4))
% hold on
% plot(t, a.real_ang_rate(t))
% xlabel('Time (sec)')
% ylabel('Angular Rate (deg/s)')
% legend('WITH Sensor Error', 'NO Sensor Error')
% title('Measured Gyroscope Signal')
% % ylim([-5 5])
% 
% %% Error Compensation
% 
% fprintf('\nGetting Corrected Signal')
% 
% corrected_signals = zeros(length(t_count), 6);
% 
% for i = 1:length(t)
%     time_i = t(i);
% 
%     accel_i = a.real_accel(time_i);
%     ang_rate_i = a.real_ang_rate(time_i);
%     curr_angle = a.real_ang(time_i);
% 
%     measured_gyro = a.measured_gyro(time_i, a.o_d_n_g_c, a.biasStabDistGyro, a.biasTempDistGyro, a.gyro_drift, a.noiseDistGyro, ang_rate_i);
% 
%     measured_accel_v = a.measured_accel_vert(time_i, a.o_d_n_a_c_v, a.biasStabDistAccel, a.biasTempDistAccel, a.accel_drift_vert, a.noiseDistAccel, accel_i, a.theta_err);
%     measured_accel_h = a.measured_accel_horz(time_i, a.o_d_n_a_c_h, a.biasStabDistAccel, a.biasTempDistAccel, a.accel_drift_horz, a.noiseDistAccel, accel_i, a.theta_err);
%     measured_accel_h = measured_accel_h / cos(curr_angle);
%     measured_accel_v = -measured_accel_v / cos(curr_angle) - a.g;
% 
%     measuredState = [measured_accel_h measured_accel_h measured_accel_v measured_gyro measured_gyro measured_gyro];
%     corrected_state = compensateError(measuredState, specs, time_i);
% 
%     corrected_signals(i, :) = corrected_state;
% end
% 
% fprintf('\nFinished Corrected Signal\n')
% 
% % Initial States
% p0 =  a.d(tspan(1))+a.pr_d;   % Platform position [m]
% p_dot0 = a.d_dot(tspan(1));   % Platform velocity [m/s]
% pr_err_accum0 = 0;            % Integral of relative position error [m*s]
% pm0 = p0;                     % Platform inetegrated position [m]
% pm_dot = p_dot0;              % Platform integrated velocity [m/s]
% pm_ddot = a.d_ddot(tspan(1)); % Platform measured acceleration [m*s^-2]
% p_theta0 = a.real_ang(tspan(1)); % Platform inertial angle [deg]
% 
% s0 = [p0; p_dot0; pr_err_accum0; pm0; pm_dot; pm_ddot];
% 
% fprintf('\nStarting Integration WITH Sensor Error.')
% fprintf("\nTime: ")
% 
% % Integrate corrected gyro signal
% 
% control_dynamics = @(t_i, state) corrected_signals((floor(t_i./dt)+1), :)';
% 
% [t_corr, int_sig_corr]= rk4_solver(control_dynamics, tspan, s0, dt);
% 
% fprintf('\nFinished Integration WITH Sensor Error.\n')
% 
% figure
% plot(t_error, corrected_signals(:, 1))
% hold on
% plot(t_error, corrected_signals(:, 3))
% hold on
% plot(t, a.real_accel(t))
% xlabel('Time (sec)')
% ylabel('Acceleration (m/s^2)')
% legend('Corrected Sensor Error (horizontal)', 'WITH Sensor Error (vertical)', 'NO Sensor Error')
% title('Measured Accelerometer Signal')
% % ylim([-5 5])
% 
% figure
% plot(t_error, corrected_signals(:, 4))
% hold on
% plot(t, a.real_ang_rate(t))
% xlabel('Time (sec)')
% ylabel('Angular Rate (deg/s)')
% legend('Corrected Sensor Error', 'NO Sensor Error')
% title('Measured Gyroscope Signal')
% % ylim([-5 5])
% 
% figure
% plot(t_corr, int_sig_corr(:, 1))
% hold on
% plot(t_corr, int_sig_corr(:, 3))
% hold on
% plot(t, a.real_vel(t))
% xlabel('Time (sec)')
% ylabel('Velocity (m/s)')
% legend('Corrected Sensor Error (horizontal)', 'Corrected Sensor Error (vertical)', 'NO Sensor Error')
% title('Integrated Velocity Signal')
% % ylim([-5 5])
% 
% figure
% plot(t_corr, int_sig_corr(:, 4))
% hold on
% plot(t, a.real_ang(t))
% xlabel('Time (sec)')
% ylabel('Angle (deg)')
% legend('Corrected Sensor Error', 'NO Sensor Error')
% title('Integrated Angle Signal')
% % ylim([-5 5])

%% Control Law Drift Compensation

finishCalibrationTime = 60; % seconds

fprintf('\nStarting Integration WITH Sensor Compensation Control Law.')
fprintf("\nTime: ")

control_dynamics = @(t_i, state) IMUDriftCorrection(t_i, a, state, finishCalibrationTime);

vel0_x = 0;
vel0_y = 0;
vel0_z = 0;
theta0 = 0;
phi0 = 0;
psi0 = 0;
p0_x = 0;
p0_y = 0;
p0_z = 0;
a_err0_x = 0;
a_err0_y = 0;
a_err0_z = 0;

s0 = [vel0_x vel0_y vel0_z theta0 phi0 psi0 p0_x p0_y p0_z a_err0_x a_err0_y a_err0_z];

[t_corr, int_state_corr]= rk4_solver(control_dynamics, tspan, s0, dt);

fprintf('\nFinished Integration WITH Sensor Compensation Control Law.\n')


%% Plotting

figure
subplot(3,1,1)
plot(t, a.theta(t))
hold on
plot(t_corr, int_state_corr(:, 4))
xlabel('Time (sec)')
ylabel('Angle (deg)')
legend('NO Sensor Error', 'Controlled Error')
title('Controlled Integrated Angle')
ylim([-0.3 0.3])

subplot(3,1,2)
plot(t, a.phi(t))
hold on
plot(t_corr, int_state_corr(:, 5))
xlabel('Time (sec)')
ylabel('Angle (deg)')
legend('NO Sensor Error', 'Controlled Error')
title('Controlled Integrated Angle')
ylim([-0.3 0.3])

subplot(3,1,3)
plot(t, a.psi(t))
hold on
plot(t_corr, int_state_corr(:, 6))
xlabel('Time (sec)')
ylabel('Angle (deg)')
legend('NO Sensor Error', 'Controlled Error')
title('Controlled Integrated Angle')
ylim([-0.3 0.3])


figure
subplot(3,1,1)
plot(t, a.real_vel_xI(t))
hold on
plot(t_corr, int_state_corr(:, 1))
xlabel('Time (sec)')
ylabel('Velocity (m/s)')
legend('NO Sensor Error', 'Controlled Error Integrated')
title('Controlled Integrated (vertical) Velocity')
% ylim([-5 5])

subplot(3,1,2)
plot(t, a.real_vel_yI(t))
hold on
plot(t_corr, int_state_corr(:, 2))
xlabel('Time (sec)')
ylabel('Velocity (m/s)')
legend('NO Sensor Error', 'Controlled Error Integrated')
title('Controlled Integrated (vertical) Velocity')
% ylim([-5 5])

subplot(3,1,3)
plot(t, a.real_vel_zI(t))
hold on
plot(t_corr, int_state_corr(:, 3))
xlabel('Time (sec)')
ylabel('Velocity (m/s)')
legend('NO Sensor Error', 'Controlled Error Integrated')
title('Controlled Integrated (vertical) Velocity')
% ylim([-5 5])

%% Plotting Acceleration


accel_m_controlled_x = diff(int_state_corr(:, 1))/dt;
accel_m_controlled_y = diff(int_state_corr(:, 2))/dt;
accel_m_controlled_z = diff(int_state_corr(:, 3))/dt;

figure
subplot(3,1,1)
plot(t, a.real_accel_xI(t))
hold on
plot(t(2:end), accel_m_controlled_x)
hold on
xline(finishCalibrationTime, 'b--')
xlabel('Time (sec)')
ylabel('Acceleration (m/s^2)')
legend('NO Sensor Error', 'Controlled Error Integrated')
title('Controlled Measured Acceleration X')
% ylim([-0.4 0.4])

subplot(3,1,2)
plot(t, a.real_accel_yI(t))
hold on
plot(t(2:end), accel_m_controlled_y)
hold on
xline(finishCalibrationTime, 'b--')
xlabel('Time (sec)')
ylabel('Acceleration (m/s^2)')
legend('NO Sensor Error', 'Controlled Error Integrated')
title('Controlled Measured Acceleration Y')
% ylim([-0.4 0.4])

subplot(3,1,3)
plot(t, a.real_accel_zI(t) + a.g)
hold on
plot(t(2:end), accel_m_controlled_z)
hold on
xline(finishCalibrationTime, 'b--')
xlabel('Time (sec)')
ylabel('Acceleration (m/s^2)')
legend('NO Sensor Error', 'Controlled Error Integrated')
title('Controlled Measured Acceleration Z')
% ylim([-0.4 0.4]) 

%% Plotting Angular Velocity


gyro_m_controlled_theta = diff(int_state_corr(:, 4))/dt;
gyro_m_controlled_phi = diff(int_state_corr(:, 5))/dt;
gyro_m_controlled_psi = diff(int_state_corr(:, 6))/dt;

figure
subplot(3,1,1)
plot(t, a.theta_dot(t))
hold on
plot(t(2:end), gyro_m_controlled_theta)
hold on
xline(finishCalibrationTime, 'b--')
xlabel('Time (sec)')
ylabel('Angular Velocity (rad/s)')
legend('NO Sensor Error', 'Controlled Error Integrated')
title('Controlled Measured Angular Velocity Theta')
% ylim([-5 5])

subplot(3,1,2)
plot(t, a.phi_dot(t))
hold on
plot(t(2:end), gyro_m_controlled_phi)
hold on
xline(finishCalibrationTime, 'b--')
xlabel('Time (sec)')
ylabel('Angular Velocity (rad/s)')
legend('NO Sensor Error', 'Controlled Error Integrated')
title('Controlled Measured Angular Velocity Phi')
% ylim([-5 5])

subplot(3,1,3)
plot(t, a.psi_dot(t))
hold on
plot(t(2:end), gyro_m_controlled_psi)
hold on
xline(finishCalibrationTime, 'b--')
xlabel('Time (sec)')
ylabel('Angular Velocity (rad/s)')
legend('NO Sensor Error', 'Controlled Error Integrated')
title('Controlled Measured Angular Velocity Psi')
% ylim([-5 5])


%% Functions

function s_dot = IMUDriftCorrection(time_i, a, prev_state, finishCalibrationTime)

    % Current states
    vel_x = prev_state(1);
    vel_y = prev_state(2);
    vel_z = prev_state(3);

    angle_theta = prev_state(4);
    angle_phi = prev_state(5);
    angle_psi = prev_state(6);

    p_x = prev_state(7);
    p_y = prev_state(8);
    p_z = prev_state(9);

    theta_err_accum = prev_state(10);
    phi_err_accum = prev_state(11);
    psi_err_accum = prev_state(12);

    specs = a.specs;

    state = [vel_x; vel_y; vel_z; angle_theta; angle_phi; angle_psi];
    state_err_accum = [p_x; p_y; p_z; theta_err_accum; phi_err_accum; psi_err_accum];

    kw = a.kw;
    kt = a.kt;

    state_0 = 0;
    state_dot_0 = 0;
    
    % Get real inertial accelerations
    a_I = [a.real_accel_xI(time_i); a.real_accel_yI(time_i); a.real_accel_zI(time_i)];
    ang_rate_i = [a.theta_dot(time_i); a.phi_dot(time_i); a.psi_dot(time_i)];
    
    % Get real angles
    theta_real = a.theta(time_i);
    phi_real = a.phi(time_i);
    psi_real = a.psi(time_i);

    % Get real sensor frame accelerations
    a_S = Rotate_I_S(a_I, phi_real, theta_real, psi_real);
    
    accel_m = a.measured_accel_3D(a, time_i, a_S);

    gyro_m = a.measured_gyro_3D(a, time_i, ang_rate_i);

    theta_use = theta_real;
    phi_use = phi_real;
    psi_use = psi_real;

    if time_i > finishCalibrationTime
        theta_use = angle_theta;
        phi_use = angle_phi;
        psi_use = angle_psi;
    end
    
    accel_I = Rotate_S_I(accel_m, phi_use, psi_use, theta_use);

    accel_state = [accel_I(1) accel_I(2) accel_I(3)+a.g];

    measuredState = [accel_state, gyro_m'];
    corrected_state = compensateError(measuredState, specs, time_i);
    
    state_dot_m = corrected_state';

    state_control = state - state_0;

    state_dot_comp = kt .* state_err_accum + kw .* state_control;

    state_dot = state_dot_m + state_dot_0 - state_dot_comp;

    s_dot = zeros(12,1);

    s_dot(1:6) = state_dot;
    s_dot(7:12) = state_control;

end
