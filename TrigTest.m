close all; clear; clc;

rng(1,"twister");

% Use same random seed

set(groot,'DefaultLineLineWidth',1)

simulationParameters;

close all;

% Redefining acceleration/gyro curves for clarity

a.real_pos = @(t) alpha*sin(beta*t) + hdeck;
a.real_vel = @(t) beta*alpha*cos(beta*t);
a.real_accel = @(t) -beta^2*alpha*sin(beta*t); % [m*s^-2]
a.real_ang = @(t) atan(beta*alpha*cos(beta*t)); % [rad]
a.real_ang_rate = @(t) (-(alpha*beta^2*sin(beta*t))./(alpha^2*beta^2*(cos(beta*t).^2)+1)); % [rad/s]

% Testing simple equations to verify integration

real_vel = @(t, y) 1/3*t.^3;
real_ang = @(t, y) 1/6*t.^4; % [deg]

real_accel = @(t, y) t.^2; % [m*s^-2]
real_ang_rate = @(t, y) 2/3*t.^3; % [deg/s]

dynamics = @(t, y) [ t.^2; 2/3*t.^3 ]; % [ real_accel real_ang_rate]

%% Sensor Model Aspects

% Simulation time
startTime = 0;
finishTime = 30;
tspan = [startTime finishTime]; % [s]

% dt = 1/imu_rate;  % [s]
dt = 0.0001;
t = (tspan(1):dt:tspan(2))';
t_count = length(t);
indeces = @(t) floor(t/dt)+1;

defineSignals

%% Plot prelim environment disturbances

x_y_z_wave = zeros(length(t), 3);

for i = 1:length(t)
    t_i = t(i);
    x_y_z_wave(i, :) = x_y_z_measurements(t_i, a);
end

figure
plot(t,x_y_z_wave(:, 1))
hold on
plot(t,x_y_z_wave(:, 2))
hold on
plot(t,x_y_z_wave(:, 3))
xlabel('Time (sec)')
ylabel('Acceleration (m/s^2)')
title('Wave Measured Signals (without sensor error)')
legend('x-axis', 'y-axis', 'z-axis')

%% Remove gravity (using true angle)

true_waves = zeros(length(t), 3);

for i = 1:length(t)

    t_i = t(i);
    state = x_y_z_wave(i, :);
    curr_angle = a.real_ang(t_i);

    true_waves(i, :) = removeGrav(a, state, curr_angle);
end


% All plots on the following graph should be the same as each other, as
% well as the same as the a.real_accel signal
figure
plot(t,true_waves(:, 1))
hold on
plot(t,true_waves(:, 2))
hold on
plot(t,true_waves(:, 3))
hold on
plot(t, a.real_accel(t))
xlabel('Time (sec)')
ylabel('Acceleration (m/s^2)')
title('True Signals (without sensor error)')
legend('x-axis', 'y-axis', 'z-axis', 'True Disturbance')

%% Get Integrated with Error Angle

fprintf('\nStarting Integration WITH Sensor Error.')
fprintf("\nTime: ")

% Running simulation with sensor error

gyro_error_signal = zeros(length(t), 1);

for i = 1:length(t)
    time_i = t(i);
    ang_rate_i = a.real_ang_rate(time_i);
    gyro_error_signal(i) = a.measured_gyro(time_i, a.o_d_n_g_c, a.biasStabDistGyro, a.biasTempDistGyro, a.gyro_drift, a.noiseDistGyro, ang_rate_i);
end

% Basic gyro signal (real angular rate with error)
figure
plot(t, a.real_ang_rate(t))
hold on
plot(t, gyro_error_signal)
% hold on
xlabel('Time (sec)')
ylabel('Angular Rate (rad/s)')
legend('Real Value', 'Gyro Measured Signal')
title('Angular Rate Real vs Measured')

p_theta0 = 0;

control_dynamics = @(t_i, state) a.measured_gyro(t_i, a.o_d_n_g_c, a.biasStabDistGyro, a.biasTempDistGyro, a.gyro_drift, a.noiseDistGyro, a.real_ang_rate(t_i));

[t_error, theta_error]= rk4_solver(control_dynamics, tspan, p_theta0, dt);

fprintf('\nFinished Integration WITH Sensor Error.\n')

figure
plot(t, 180/pi*a.real_ang(t))
hold on
plot(t, 180/pi*theta_error)
xlabel('Time (sec)')
ylabel('Angle (deg)')
legend('NO Sensor Error', 'WITH Sensor Error')
title('Basic Integrated Angle Plot')
% ylim([-30 30])

%% Remove gravity (using error angle)

error_waves = zeros(length(t), 3);

for i = 1:length(t)

    t_i = t(i);
    state = x_y_z_wave(i, :);
    error_angle = theta_error(i);

    error_waves(i, :) = removeGrav(a, state, error_angle);
end


% All plots on the following graph should be the same as each other, as
% well as the same as the a.real_accel signal
figure
plot(t,error_waves(:, 1))
hold on
plot(t,error_waves(:, 2))
hold on
plot(t,error_waves(:, 3))
hold on
plot(t, a.real_accel(t))
xlabel('Time (sec)')
ylabel('Acceleration (m/s^2)')
title('Error Signals (WITH sensor error)')
legend('x-axis', 'y-axis', 'z-axis', 'True Disturbance')


%% Correct Signal with Error

fprintf('\nStarting Integration WITH Sensor Compensation Control Law.')
fprintf("\nTime: ")

control_dynamics = @(t_i, state) GyroDriftCorrection(t_i, a, state);

theta_0 = 0;
theta_dot_0 = 0;

s0 = [theta_0 theta_dot_0];

[t_corr, theta_corr]= rk4_solver(control_dynamics, tspan, s0, dt);

fprintf('\nFinished Integration WITH Sensor Compensation Control Law.\n')

figure
plot(t, 180/pi*a.real_ang(t))
hold on
plot(t, 180/pi*theta_error)
hold on
plot(t, 180/pi*theta_corr(:, 1))
xlabel('Time (sec)')
ylabel('Angle (deg)')
legend('NO Sensor Error', 'WITH Sensor Error', 'Corrected Gyro Signal')
title('Basic Integrated Angle Plot')
ylim([-30 30])

%% Remove gravity (using corrected angle)

corr_waves = zeros(length(t), 3);

for i = 1:length(t)

    t_i = t(i);
    state = x_y_z_wave(i, :);
    int_angle = theta_corr(i, 1);

    corr_waves(i, :) = removeGrav(a, state, int_angle);
end


% All plots on the following graph should be the same as each other, as
% well as the same as the a.real_accel signal
figure
plot(t,corr_waves(:, 1))
hold on
plot(t,corr_waves(:, 2))
hold on
plot(t,corr_waves(:, 3))
hold on
plot(t, a.real_accel(t))
xlabel('Time (sec)')
ylabel('Acceleration (m/s^2)')
title('Fixed Signals (with CORRECTED sensor error)')
legend('x-axis', 'y-axis', 'z-axis', 'True Disturbance')


%% Functions


function signals = x_y_z_measurements(t, a)
    
    accel_true = a.real_accel(t);
    angle_true = a.real_ang(t);

    a_signal_x = accel_true * cos(angle_true);
    a_signal_y = accel_true * cos(angle_true);
    a_signal_z = -(accel_true + a.g) * cos(angle_true);

    signals = [a_signal_x a_signal_y a_signal_z];

end


function true_state = removeGrav(a, state, angle)

    a_x = state(1) / cos(angle);
    a_y = state(2) / cos(angle);
    a_z = -state(3) / cos(angle) - a.g;

    true_state = [a_x a_y a_z];

end

function s_dot = GyroDriftCorrection(time_i, a, prev_state)


    % Current states
    theta = prev_state(1);
    theta_err_accum = prev_state(2);

    specs = a.specs;

    % gyro_noise_std = specs.gyro_noiseDensity * sqrt(specs.gyro_bandwidth);
    % theta_0 = normrnd(0, gyro_noise_std);

    kw = a.kw; % 
    kt = a.kt; % 

    theta_0 = 0;
    theta_dot_0 = 0;

    ang_rate_i = a.real_ang_rate(time_i);
    measured_g = a.measured_gyro(time_i, a.o_d_n_g_c, a.biasStabDistGyro, a.biasTempDistGyro, a.gyro_drift, a.noiseDistGyro, ang_rate_i);
   
    measuredState = [0 0 0 measured_g measured_g measured_g];
    corrected_state = compensateError(measuredState, specs, time_i);

    theta_dot_m = corrected_state(4);

    theta_control = theta - theta_0;

    theta_dot_comp = kt * theta_err_accum + kw * theta_control;

    theta_dot = theta_dot_m + theta_dot_0 - theta_dot_comp;

    s_dot = zeros(2,1);

    s_dot(1) = theta_dot;
    s_dot(2) = theta_control;

end
