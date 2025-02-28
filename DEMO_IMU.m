close all; clear; clc;

rng(1,"twister");

% Use same random seed

set(groot,'DefaultLineLineWidth',1)

simulationParameters;

close all;

% Redefining acceleration/gyro curves for clarity

a.real_pos = @(t) alpha*sin(beta*t) + hdeck;
a.real_vel = @(t) beta*alpha*cos(beta*t);
a.real_accel = @(t, y) -beta^2*alpha*sin(beta*t); % [m*s^-2]
a.real_ang = @(t) atan(beta*alpha*cos(beta*t)); % [rad]
a.real_ang_rate = @(t, y) (-(alpha*beta^2*sin(beta*t))./(alpha^2*beta^2*(cos(beta*t).^2)+1)); % [rad/s]

%% Sensor Model Aspects

% Simulation time
startTime = 0;
finishTime = 60;
tspan = [startTime finishTime]; % [s]

% dt = 1/imu_rate;  % [s]
dt = 0.0001;
t = (tspan(1):dt:tspan(2))';
t_count = length(t);
indeces = @(t) floor(t/dt)+1;

defineSignals

scale_w = 1;
scale_t = 1;

a.kw = scale_w*a.beta_min; % 
a.kt = scale_t*a.beta_max; % 

%% Get Error Signal

fprintf('\nStarting Integration WITH Sensor Error.')
fprintf("\nTime: ")

% Running simulation with sensor error

error_signals = zeros(t_count, 6);

for i = 1:length(t)
    time_i = t(i);

    ang_rate_i = a.real_ang_rate(time_i);
    accel_i = a.real_accel(time_i);
    
    measured_gyro = a.measured_gyro(time_i, a.o_d_n_g_c, a.biasStabDistGyro, a.biasTempDistGyro, a.gyro_drift, a.noiseDistGyro, ang_rate_i);
    measured_accel_v = a.measured_accel_vert(time_i, a.o_d_n_a_c_v, a.biasStabDistAccel, a.biasTempDistAccel, a.accel_drift_vert, a.noiseDistAccel, accel_i, a.theta_err);
    measured_accel_h = a.measured_accel_horz(time_i, a.o_d_n_a_c_h, a.biasStabDistAccel, a.biasTempDistAccel, a.accel_drift_horz, a.noiseDistAccel, accel_i, a.theta_err);

    measuredState = [measured_accel_h measured_accel_h measured_accel_v measured_gyro measured_gyro measured_gyro];

    error_signals(i, :) = measuredState;
end

theta0_x = 0;
theta0_y = 0;
theta0_z = 0;
vel0_x = 0;
vel0_y = 0;
vel0_z = 0;

s0 = [theta0_x theta0_y theta0_z vel0_x vel0_y vel0_z];

control_dynamics = @(t_i, state) error_signals(floor(t_i./dt)+1)';

[t_error, int_signals_error]= rk4_solver(control_dynamics, tspan, s0, dt);

fprintf('\nFinished Integration WITH Sensor Error.\n')

figure
plot(t_error, error_signals(:, 1))
hold on
plot(t_error, error_signals(:, 3))
hold on
plot(t, a.real_accel(t))
xlabel('Time (sec)')
ylabel('Acceleration (m/s^2)')
legend('WITH Sensor Error (horizontal)', 'WITH Sensor Error (vertical)', 'NO Sensor Error')
title('Measured Accelerometer Signal')
% ylim([-5 5])

figure
plot(t_error, error_signals(:, 4))
hold on
plot(t, a.real_ang_rate(t))
xlabel('Time (sec)')
ylabel('Angular Rate (deg/s)')
legend('WITH Sensor Error', 'NO Sensor Error')
title('Measured Gyroscope Signal')
% ylim([-5 5])

%% Error Compensation

fprintf('\nGetting Corrected Signal')

corrected_signals = zeros(length(t_count), 6);

for i = 1:length(t)
    time_i = t(i);

    accel_i = a.real_accel(time_i);
    ang_rate_i = a.real_ang_rate(time_i);
    curr_angle = a.real_ang(time_i);
    
    measured_gyro = a.measured_gyro(time_i, a.o_d_n_g_c, a.biasStabDistGyro, a.biasTempDistGyro, a.gyro_drift, a.noiseDistGyro, ang_rate_i);
    
    measured_accel_v = a.measured_accel_vert(time_i, a.o_d_n_a_c_v, a.biasStabDistAccel, a.biasTempDistAccel, a.accel_drift_vert, a.noiseDistAccel, accel_i, a.theta_err);
    measured_accel_h = a.measured_accel_horz(time_i, a.o_d_n_a_c_h, a.biasStabDistAccel, a.biasTempDistAccel, a.accel_drift_horz, a.noiseDistAccel, accel_i, a.theta_err);
    measured_accel_h = measured_accel_h / cos(curr_angle);
    measured_accel_v = -measured_accel_v / cos(curr_angle) - a.g;

    measuredState = [measured_accel_h measured_accel_h measured_accel_v measured_gyro measured_gyro measured_gyro];
    corrected_state = compensateError(measuredState, specs, time_i);

    corrected_signals(i, :) = corrected_state;
end

fprintf('\nFinished Corrected Signal\n')

% Initial States
p0 =  a.d(tspan(1))+a.pr_d;   % Platform position [m]
p_dot0 = a.d_dot(tspan(1));   % Platform velocity [m/s]
pr_err_accum0 = 0;            % Integral of relative position error [m*s]
pm0 = p0;                     % Platform inetegrated position [m]
pm_dot = p_dot0;              % Platform integrated velocity [m/s]
pm_ddot = a.d_ddot(tspan(1)); % Platform measured acceleration [m*s^-2]
p_theta0 = a.real_ang(tspan(1)); % Platform inertial angle [deg]

s0 = [p0; p_dot0; pr_err_accum0; pm0; pm_dot; pm_ddot];

fprintf('\nStarting Integration WITH Sensor Error.')
fprintf("\nTime: ")

% Integrate corrected gyro signal

control_dynamics = @(t_i, state) corrected_signals((floor(t_i./dt)+1), :)';

[t_corr, int_sig_corr]= rk4_solver(control_dynamics, tspan, s0, dt);

fprintf('\nFinished Integration WITH Sensor Error.\n')

figure
plot(t_error, corrected_signals(:, 1))
hold on
plot(t_error, corrected_signals(:, 3))
hold on
plot(t, a.real_accel(t))
xlabel('Time (sec)')
ylabel('Acceleration (m/s^2)')
legend('Corrected Sensor Error (horizontal)', 'WITH Sensor Error (vertical)', 'NO Sensor Error')
title('Measured Accelerometer Signal')
% ylim([-5 5])

figure
plot(t_error, corrected_signals(:, 4))
hold on
plot(t, a.real_ang_rate(t))
xlabel('Time (sec)')
ylabel('Angular Rate (deg/s)')
legend('Corrected Sensor Error', 'NO Sensor Error')
title('Measured Gyroscope Signal')
% ylim([-5 5])

figure
plot(t_corr, int_sig_corr(:, 1))
hold on
plot(t_corr, int_sig_corr(:, 3))
hold on
plot(t, a.real_vel(t))
xlabel('Time (sec)')
ylabel('Velocity (m/s)')
legend('Corrected Sensor Error (horizontal)', 'Corrected Sensor Error (vertical)', 'NO Sensor Error')
title('Integrated Velocity Signal')
% ylim([-5 5])

figure
plot(t_corr, int_sig_corr(:, 4))
hold on
plot(t, a.real_ang(t))
xlabel('Time (sec)')
ylabel('Angle (deg)')
legend('Corrected Sensor Error', 'NO Sensor Error')
title('Integrated Angle Signal')
% ylim([-5 5])

%% Control Law Drift Compensation

fprintf('\nStarting Integration WITH Sensor Compensation Control Law.')
fprintf("\nTime: ")

control_dynamics = @(t_i, state) IMUDriftCorrection(t_i, a, state);

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

figure
plot(t, a.real_vel(t))
hold on
plot(t_corr, int_state_corr(:, 1))
hold on
plot(t_corr, int_state_corr(:, 3))
xlabel('Time (sec)')
ylabel('Velocity (m/s)')
legend('NO Sensor Error', 'Controlled Error (horizontal)', 'Controlled Error (vertical)')
title('Controlled Integrated Velocity')
% ylim([-5 5])

figure
plot(t, a.real_ang(t))
hold on
plot(t_corr, int_state_corr(:, 4))
xlabel('Time (sec)')
ylabel('Angle (deg)')
legend('NO Sensor Error', 'Controlled Error')
title('Controlled Integrated Angle')
% ylim([-5 5])

figure
plot(t, a.real_vel(t))
hold on
plot(t, int_signals_error(:, 3))
hold on
plot(t_corr, int_state_corr(:, 3))
xlabel('Time (sec)')
ylabel('Velocity (m/s)')
legend('NO Sensor Error', 'Integrated Raw Signal', 'Controlled Error Integrated')
title('Controlled Integrated (vertical) Velocity')
% ylim([-5 5])

figure
plot(t, a.real_ang(t))
hold on
plot(t, int_signals_error(:, 4))
hold on
plot(t_corr, int_state_corr(:, 4))
xlabel('Time (sec)')
ylabel('Angle (deg)')
legend('NO Sensor Error', 'Integrated Raw Signal', 'Controlled Error Integrated')
title('Controlled Integrated Angle')
% ylim([-5 5])

function s_dot = IMUDriftCorrection(time_i, a, prev_state)

    % Current states
    vel_x = prev_state(1);
    vel_y = prev_state(2);
    vel_z = prev_state(3);
    theta = prev_state(4);
    phi = prev_state(5);
    psi = prev_state(6);
    p_x = prev_state(7);
    p_y = prev_state(8);
    p_z = prev_state(9);
    theta_err_accum = prev_state(10);
    phi_err_accum = prev_state(11);
    psi_err_accum = prev_state(12);

    specs = a.specs;

    state = [vel_x; vel_y; vel_z; theta; phi; psi];
    state_err_accum = [p_x; p_y; p_z; theta_err_accum; phi_err_accum; psi_err_accum];

    kw = a.kw;
    kt = a.kt;

    state_0 = 0;
    state_dot_0 = 0;

    accel_i = a.real_accel(time_i);
    ang_rate_i = a.real_ang_rate(time_i);
    curr_angle = a.real_ang(time_i);

    a_h = a.measured_accel_vert(time_i, a.o_d_n_a_c_v, a.biasStabDistAccel, a.biasTempDistAccel, a.accel_drift_vert, a.noiseDistAccel, accel_i, a.real_ang);
    a_v = a.measured_accel_horz(time_i, a.o_d_n_a_c_h, a.biasStabDistAccel, a.biasTempDistAccel, a.accel_drift_horz, a.noiseDistAccel, accel_i, a.real_ang);
    
    accel_state = [a_h a_h a_v];

    corr_a = AccelRemoveGrav(accel_state, curr_angle, a);

    measured_a_h = corr_a(1);
    measured_a_v = corr_a(3);

    measured_g = a.measured_gyro(time_i, a.o_d_n_g_c, a.biasStabDistGyro, a.biasTempDistGyro, a.gyro_drift, a.noiseDistGyro, ang_rate_i);
   
    measuredState = [measured_a_h measured_a_h measured_a_v measured_g measured_g measured_g];
    corrected_state = compensateError(measuredState, specs, time_i);
    
    state_dot_m = corrected_state';

    state_control = state - state_0;

    state_dot_comp = kt * state_err_accum + kw * state_control;

    state_dot = state_dot_m + state_dot_0 - state_dot_comp;

    s_dot = zeros(12,1);

    s_dot(1:6) = state_dot;
    s_dot(7:12) = state_control;

end