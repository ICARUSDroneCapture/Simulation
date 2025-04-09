close all; clear; clc;

rng(1,"twister");

% Use same random seed

set(groot,'DefaultLineLineWidth',1)

simulationParameters;

% Sensor Drift Control
a.kw = 1200; % 
a.kt = 0.6; % 

close all;

% Redefining acceleration/gyro curves for clarity

a.real_pos = @(t) alpha*sin(beta*t) + hdeck;
a.real_vel = @(t) beta*alpha*cos(beta*t);
a.real_accel = @(t, y) -beta^2*alpha*sin(beta*t); % [m*s^-2]
% a.real_accel = @(t, y) 0*t; % [m*s^-2]
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

fprintf('\nIntegrating Acceleration (no Sensor Error)')
fprintf("\nTime: ")

% Initial Angle
p_vel0 = a.real_vel(tspan(1)); % Platform inertial angle [deg]

[t, vel_base]= rk4_solver(a.real_accel, tspan, p_vel0, dt);

% Integration Verification Plot
figure
plot(t, vel_base)
hold on
plot(t, a.real_vel(t))
xlabel('Time (sec)')
ylabel('Velocity (m/s)')
legend('Integrated Velocity', 'Velocity Equation')
title('Basic Integrated Velocity Plot')

fprintf('\nFinished Integration with NO Sensor Error.\n')
fprintf('\nStarting Integration WITH Sensor Error.')
fprintf("\nTime: ")

% Running simulation with sensor error

accel_error_signal = zeros(length(t), 2);

for i = 1:length(t)
    time_i = t(i);

    accel_i = a.real_accel(time_i);
    curr_angle = a.real_ang(time_i);

    measured_a_h = a.measured_accel_vert(time_i, a.o_d_n_a_c_v, a.biasStabDistAccel, a.biasTempDistAccel, a.accel_drift_vert, a.noiseDistAccel, accel_i, a.real_ang);
    measured_a_v = a.measured_accel_horz(time_i, a.o_d_n_a_c_h, a.biasStabDistAccel, a.biasTempDistAccel, a.accel_drift_horz, a.noiseDistAccel, accel_i, a.real_ang);
    
    accel_error_signal(i, 1) = measured_a_h / cos(curr_angle);
    accel_error_signal(i, 2) = -measured_a_v / cos(curr_angle) - a.g;
end

% Basic gyro signal (real angular rate with error)
figure
plot(t, a.real_accel(t))
hold on
plot(t, accel_error_signal(:, 1))
hold on
plot(t, accel_error_signal(:, 2))
xlabel('Time (sec)')
ylabel('Acceleration (m/s^2)')
legend('Real Value', 'Accel Measured Signal (vertical)', 'Accel Measured Signal (horizontal)')
title('Acceleration Real vs Measured')

control_dynamics = @(t_i, state) accel_error_signal((floor(t_i./dt)+1), :)';

[t_error, vel_error]= rk4_solver(control_dynamics, tspan, [p_vel0 p_vel0], dt);

fprintf('\nFinished Integration WITH Sensor Error.\n')

%% Plotting Platform Velocity

figure
plot(t, a.real_vel(t))
hold on
plot(t_error, vel_error(:, 1))
hold on
plot(t_error, vel_error(:, 2))
xlabel('Time (sec)')
ylabel('Velocity (m/s)')
legend('NO Sensor Error', 'WITH Sensor Error (vertical)', 'WITH Sensor Error (horizontal)')
title('Basic Integrated Velocity Plot')

%% Getting and Plotting Velocity Error

vel_diff = vel_error(:, 1) - vel_base(:, 1);

sz = 2;
plot_scale = 10;

figure
scatter(t_error, vel_diff, sz, 'filled', displayName="Velocity Error")
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

[t_corr, vel_corr]= rk4_solver(control_dynamics, tspan, s0, dt);

a_der = diff(vel_corr(:, 3)) / dt;
figure
plot(t(2:end), a_der)
hold on
plot(t, a.real_accel(t))
xlabel('Time (sec)')
ylabel('Acceleration (m/s^2)')
legend('Controlled Acceleration', 'Real Deck Acceleration')

fprintf('\nFinished Integration WITH Sensor Compensation Control Law.\n')

figure
plot(t, vel_base)
hold on
plot(t_error, vel_error(:, 1))
hold on
plot(t_corr, vel_corr(:, 3))
xlabel('Time (sec)')
ylabel('Velocity (m/s)')
legend('NO Sensor Error', 'WITH Sensor Error', 'Corrected (vertical) Accel Signal')
title('Basic Integrated Velocity Plot')
% ylim([-5 5])

% %% Gravity Check
% 
% fprintf('\nStarting Integration WITH Sensor Compensation Control Law.')
% fprintf("\nTime: ")
% 
% control_dynamics = @(t_i, state) AccelGravAddRemove(t_i, a, state);
% 
% vel_0 = 0;
% accel_0 = 0;
% 
% s0 = [vel_0 vel_0 vel_0 accel_0 accel_0 accel_0];
% 
% [t_corr, vel_corr]= rk4_solver(control_dynamics, tspan, s0, dt);
% 
% a_der = diff(vel_corr(:, 3)) / dt;
% figure
% plot(t(2:end), a_der)
% hold on
% plot(t, a.real_accel(t))
% xlabel('Time (sec)')
% ylabel('Acceleration (m/s^2)')
% legend('Controlled Acceleration', 'Real Deck Acceleration')
% 
% fprintf('\nFinished Integration WITH Sensor Compensation Control Law.\n')
% 
% figure
% plot(t, vel_base)
% hold on
% plot(t_error, vel_error(:, 1))
% hold on
% plot(t_corr, vel_corr(:, 3))
% xlabel('Time (sec)')
% ylabel('Velocity (m/s)')
% legend('NO Sensor Error', 'WITH Sensor Error', 'Corrected (vertical) Accel Signal')
% title('Basic Integrated Velocity Plot')
% % ylim([-5 5])

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

    kw = a.kw;
    kt = a.kt;

    vel_0 = 0;
    vel_dot_0 = 0;

    accel_i = a.real_accel(time_i);
    curr_angle = a.real_ang(time_i);

    a_h = a.measured_accel_vert(time_i, a.o_d_n_a_c_v, a.biasStabDistAccel, a.biasTempDistAccel, a.accel_drift_vert, a.noiseDistAccel, accel_i, a.real_ang);
    a_v = a.measured_accel_horz(time_i, a.o_d_n_a_c_h, a.biasStabDistAccel, a.biasTempDistAccel, a.accel_drift_horz, a.noiseDistAccel, accel_i, a.real_ang);
    
    accel_state = [a_h a_h a_v];

    corr_a = AccelRemoveGrav(accel_state, curr_angle, a);

    measured_a_h = corr_a(1);
    measured_a_v = corr_a(3);

    measuredState = [measured_a_h measured_a_h measured_a_v 0 0 0];
    corrected_state = compensateError(measuredState, specs, time_i);
    
    vel_dot_m = [corrected_state(1); corrected_state(2); corrected_state(3)];

    vel_control = vel - vel_0;

    vel_dot_comp = kt * vel_err_accum + kw * vel_control;

    vel_dot = vel_dot_m + vel_dot_0 - vel_dot_comp;

    s_dot = zeros(6,1);

    s_dot(1:3) = vel_dot;
    s_dot(4:6) = vel_control;

end


function s_dot = AccelGravAddRemove(time_i, a, prev_state)

    % Current states
    vel_x = prev_state(1);
    vel_y = prev_state(2);
    vel_z = prev_state(3);

    specs = a.specs;

    vel = [vel_x; vel_y; vel_z];

    accel_i = a.real_accel(time_i);
    curr_angle = a.real_ang(time_i);

    a_h = a.measured_accel_vert(time_i, a.o_d_n_a_c_v, a.biasStabDistAccel, a.biasTempDistAccel, a.accel_drift_vert, a.noiseDistAccel, accel_i, a.real_ang);
    a_v = a.measured_accel_horz(time_i, a.o_d_n_a_c_h, a.biasStabDistAccel, a.biasTempDistAccel, a.accel_drift_horz, a.noiseDistAccel, accel_i, a.real_ang);
    
    accel_state = [a_h a_h a_v];

    corr_a = AccelRemoveGrav(accel_state, curr_angle, a);

    measured_a_h = corr_a(1);
    measured_a_v = corr_a(3);

    measuredState = [measured_a_h measured_a_h measured_a_v 0 0 0];
    corrected_state = compensateError(measuredState, specs, time_i);
    
    vel_dot_m = [corrected_state(1); corrected_state(2); corrected_state(3)];

    s_dot = zeros(3,1);

    s_dot(1:3) = vel_dot_m;

end
