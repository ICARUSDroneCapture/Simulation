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
a.real_ang_rate = @(t, y) (-(alpha*beta^2*sin(beta*t))./(alpha^2*beta^2*(cos(beta*t).^2)+1)); % [rad/s]

%% Sensor Model Aspects

imx_5_specs
% test3

% ----------- Set gyro parameters here for testing if desired -------------

% specs.accel_resolution = 0.122 / 1000 * a.g; % m/s
% specs.accel_resolution = 1*10^-30; % rad/s
% specs.accel_noiseDensity = 0 * (10^-3) / 180*pi; % rad/s/sqrt(Hz)
% specs.accel_temp_bias = 0 / 180*pi; % rad/s RMS

% Gyro Specs
% specs.b_g = 0 / 3600 / 180*pi; % Time Varying Bias (rad/s)
% specs.ARW = 0 / 60 / 180*pi; % Angle Random Walk (rad/sqrt(s))

% -------------------------------------------------------------------------

specs.g = a.g;
a.specs = specs;

k = specs.k;
nonlinearity = specs.dk;

V_err_0 = specs.V_err_0;
P_err_0 = specs.P_err_0;

accel_resolution = specs.accel_resolution;
accel_samplingRate = specs.accel_samplingRate;
accel_noiseDensity = specs.accel_noiseDensity;
accel_bandwidth = specs.accel_bandwidth;
accel_temp_bias = specs.accel_temp_bias;

gyro_resolution = specs.gyro_resolution;
gyro_samplingRate = specs.gyro_samplingRate;
gyro_noiseDensity = specs.gyro_noiseDensity;
gyro_bandwidth = specs.gyro_bandwidth;
gyro_temp_bias = specs.gyro_temp_bias;

b_a = specs.b_a;
VRW = specs.VRW;
b_g = specs.b_g; 
ARW = specs.ARW;

g = a.g;

sz = 3;

% Simulation time
startTime = 0;
finishTime = 3*60;
% finishTime = 30;
tspan = [startTime finishTime]; % [s]

%solving the system
freq = 160; % Hz
dt = 1 / freq; % Timestep (s)

imu_accel_sampling_rate = 3*specs.accel_bandwidth;
imu_gyro_sampling_rate = 3*specs.gyro_bandwidth;

imu_rate = max(imu_accel_sampling_rate, imu_gyro_sampling_rate);

% dt = 1/imu_rate;  % [s]
% dt = 0.0001;
t = (tspan(1):dt:tspan(2))';
t_count = length(t);
indeces = @(t) floor(t/dt)+1;

%% Standard Noise Visualization

accel_noise_std = specs.accel_noiseDensity * sqrt(specs.accel_bandwidth); % Noise standard deviation (m/s^2)
gyro_noise_std = specs.gyro_noiseDensity * sqrt(specs.gyro_bandwidth); % Noise standard deviation (rad/s)

% Getting noise values based on standard distribution
a.noiseDistAccel = @(t) accel_noise_std*randn(length(t), 1);
a.noiseDistGyro = @(t, y) gyro_noise_std*randn(length(t), 1);

% Noisy deck signal
a.noisyAccelCurve = @(t, real_accel) real_accel(t) + accel_noise_std*randn(length(t),1);
a.noisyGyroCurve = @(t, real_ang_rate) real_ang_rate(t) + gyro_noise_std*randn(length(t),1);

figure
plot(t, a.noisyAccelCurve(t, a.real_accel));
hold on
plot(t, a.real_accel(t))
xlabel('Time (sec)')
ylabel('Acceleration (m/s^2)')
title('Accelerometer Noise Effect')
legend('Measured Accelerometer Signal', 'Real Acceleration')

% Test integration with just this source of noise

% Integrate gyro noise
noise_handle = @(t, y) a.noisyAccelCurve(t, a.real_accel);
[t_vrw, VRW_curve]= rk4_solver(noise_handle, tspan, a.real_vel(tspan(1)), dt);

% Integration Verification Plot
figure
plot(t_vrw, VRW_curve)
hold on
plot(t_vrw, a.real_accel(t_vrw))
xlabel('Time (sec)')
ylabel('Velocity (m/s)')
title('Integrated Noisy Accel Curve - Velocity Random Walk')
legend('Integrated with Noise', 'Integrated without Noise')

%% Adding bias instability

% But first, turn-on bias
accel_turn_on_bias_offset = 0.01;
gyro_turn_on_bias_offset = 0.01;

sigma_bias_accel = b_a * sqrt(imu_rate);
sigma_bias_gyro = b_g * sqrt(imu_rate);

bias_noise_accel = sigma_bias_accel * randn(length(t), 1);
bias_noise_gyro = sigma_bias_gyro * randn(length(t), 1);

scale = 100000;

freq_pass_accel = specs.accel_bandwidth / scale;
b_accel = (1 - exp(-freq_pass_accel*dt)) / freq_pass_accel;

freq_pass_gyro = specs.gyro_bandwidth / scale;
b_gyro = (1 - exp(-freq_pass_gyro*dt)) / freq_pass_gyro;

B = [0 0 b_accel];
A = [1 -2*(exp(-freq_pass_accel*dt)) (exp(-freq_pass_accel*dt))^2];
y_a_filter = filter(B,A,bias_noise_accel);

B = [0 0 b_gyro];
A = [1 -2*(exp(-freq_pass_gyro*dt)) (exp(-freq_pass_gyro*dt))^2];
y_g_filter = filter(B,A,bias_noise_gyro);

accel_bias_norm = normalize(y_a_filter, "range")*b_a;
gyro_bias_norm = normalize(y_g_filter, "range")*b_g;

% Bias Error over Temp ----------------------------------------------------

sensor_temp_range = 85 - (-40); % deg C

% Setting environmental temperature change to full range for worst case
env_temp_range = 6;

bias_offset_scale = env_temp_range / sensor_temp_range;

accel_temp_std = bias_offset_scale * accel_temp_bias;
gyro_temp_std = bias_offset_scale * gyro_temp_bias;

sigma_bias_temp_accel = accel_temp_std * sqrt(imu_rate);
sigma_bias_temp_gyro = gyro_temp_std * sqrt(imu_rate);

bias_temp_noise_accel = sigma_bias_temp_accel * randn(length(t), 1);
bias_temp__gyro = sigma_bias_temp_gyro * randn(length(t), 1);

temp_freq_pass = 1/200;

freq_pass_accel = temp_freq_pass;
b_accel = (1 - exp(-freq_pass_accel*dt)) / freq_pass_accel;

freq_pass_gyro = temp_freq_pass;
b_gyro = (1 - exp(-freq_pass_gyro*dt)) / freq_pass_gyro;

B = [0 0 b_accel];
A = [1 -2*(exp(-freq_pass_accel*dt)) (exp(-freq_pass_accel*dt))^2];
y_a_filter_temp = filter(B,A,bias_temp_noise_accel);

B = [0 0 b_gyro];
A = [1 -2*(exp(-freq_pass_gyro*dt)) (exp(-freq_pass_gyro*dt))^2];
y_g_filter_temp = filter(B,A,bias_temp__gyro);

accel_bias_temp_norm = normalize(y_a_filter_temp, "range")*accel_temp_std;
gyro_bias_temp_norm = normalize(y_g_filter_temp, "range")*gyro_temp_std;

a.biasStabDistAccel = @(t) accel_bias_norm(floor(t./dt)+1);
a.biasStabDistGyro = @(t) gyro_bias_norm(floor(t./dt)+1);

a.biasTempDistAccel = @(t) accel_bias_temp_norm(floor(t./dt)+1);
a.biasTempDistGyro = @(t) gyro_bias_temp_norm(floor(t./dt)+1);

a.biasedAccelCurve = @(t, real_accel, biasStabDistAccel, biasTempDistAccel) real_accel(t) + a.biasStabDistAccel(t) + a.biasTempDistAccel(t) + accel_turn_on_bias_offset;
a.biasedGyroCurve = @(t, real_ang_rate, biasStabDistGyro, biasTempDistGyro) real_ang_rate(t) + a.biasStabDistGyro(t) + a.biasTempDistGyro(t) + gyro_turn_on_bias_offset;

figure
plot(t, a.biasStabDistAccel(t) + a.biasTempDistAccel(t));
xlabel('Time (sec)')
ylabel('Acceleration (m/s^2)')
title('Accelerometer Drift (Bias instability and Temperature Drift)')

a.theta_err = @(t) a.biasStabDistGyro(t).*t + ARW.*sqrt(t);

a.accel_drift_vert = @(t, accel, theta_err, biasStabDistAccel, biasTempDistAccel) (1 + k)*accel + biasStabDistAccel(t)  + biasTempDistAccel(t) + g*(1-cos(theta_err(t))) + accel_turn_on_bias_offset;
a.accel_drift_horz = @(t, accel, theta_err, biasStabDistAccel, biasTempDistAccel) (1 + k)*accel + biasStabDistAccel(t)  + biasTempDistAccel(t)  + g*sin(theta_err(t)) + accel_turn_on_bias_offset;

a.gyro_drift = @(t, ang_rate, biasStabDistGyro, biasTempDistGyro) (1 + k)*ang_rate + a.biasStabDistGyro(t) + biasTempDistGyro(t) + gyro_turn_on_bias_offset;

figure
plot(t, a.accel_drift_vert(t, a.real_accel(t), a.theta_err, a.biasStabDistAccel, a.biasTempDistAccel));
hold on
plot(t, a.accel_drift_horz(t, a.real_accel(t), a.theta_err, a.biasStabDistAccel, a.biasTempDistAccel));
hold on
plot(t, a.real_accel(t))
xlabel('Time (sec)')
ylabel('Acceleration (m/s^2)')
title('Accelerometer Drifting Signal')
legend('Measured (vertical) Accelerometer Signal', 'Measured (horizontal) Accelerometer Signal', 'Real Acceleration')

% Test integration with just this source of noise

% Integrate gyro noise
bias_handle = @(t, y) a.biasedAccelCurve(t, a.real_accel, a.biasStabDistAccel, a.biasTempDistAccel);
[t_bias, biased_curve]= rk4_solver(bias_handle, tspan, a.real_vel(tspan(1)), dt);

% Integration Verification Plot
figure
plot(t_bias, biased_curve)
hold on
plot(t_bias, a.real_vel(t_bias))
xlabel('Time (sec)')
ylabel('Velocity (m/s)')
title('Integrated Biased Accel Curve')
legend('Integrated with Bias', 'Integrated without Bias')

%% Combining bias drift, temp bias, turn-on bias, noise, and quantization

% o_d_n_a_c_v: offset drifting noisy accel_curve vert
% o_d_n_a_c_h: offset drifting noisy accel curve horz
% o_d_n_g_c: offset drifting noisy gyro curve

a.o_d_n_a_c_v = @(t, biasStabDistAccel, biasTempDistAccel, accel_drift_vert, noiseDistAccel, real_accel, theta_err) accel_drift_vert(t, real_accel, theta_err, biasStabDistAccel, biasTempDistAccel) + noiseDistAccel(t);
a.o_d_n_a_c_h = @(t, biasStabDistAccel, biasTempDistAccel, accel_drift_horz, noiseDistAccel, real_accel, theta_err) accel_drift_horz(t, real_accel, theta_err, biasStabDistAccel, biasTempDistAccel) + noiseDistAccel(t);
a.o_d_n_g_c = @(t, biasStabDistGyro, biasTempDistGyro, gyro_drift, noiseDistGyro, real_ang_rate) gyro_drift(t, real_ang_rate, biasStabDistGyro, biasTempDistGyro) + noiseDistGyro(t);

a.measured_accel_vert = @(t, o_d_n_a_c_v, biasStabDistAccel, biasTempDistAccel, accel_drift_vert, noiseDistAccel, real_accel, theta_err) accel_resolution*floor(o_d_n_a_c_v(t, biasStabDistAccel, biasTempDistAccel, accel_drift_vert, noiseDistAccel, real_accel, theta_err)/accel_resolution);
a.measured_accel_horz = @(t, o_d_n_a_c_h, biasStabDistAccel, biasTempDistAccel, accel_drift_horz, noiseDistAccel, real_accel, theta_err) accel_resolution*floor(o_d_n_a_c_h(t, biasStabDistAccel, biasTempDistAccel, accel_drift_horz, noiseDistAccel, real_accel, theta_err)/accel_resolution);
a.measured_gyro = @(t, o_d_n_g_c, biasStabDistGyro, biasTempDistGyro, gyro_drift, noiseDistGyro, real_ang_rate) gyro_resolution*floor(o_d_n_g_c(t, biasStabDistGyro, biasTempDistGyro, gyro_drift, noiseDistGyro, real_ang_rate)/gyro_resolution);

% measured_accel_vert = a.measured_accel_vert(t, a.o_d_n_a_c_v, a.biasStabDistAccel, a.biasTempDistAccel, a.accel_drift_vert, a.noiseDistAccel, a.real_accel, a.theta_err);
% measured_accel_horz = a.measured_accel_horz(t, a.o_d_n_a_c_h, a.biasStabDistAccel, a.biasTempDistAccel, a.accel_drift_horz, a.noiseDistAccel, a.real_accel, a.theta_err);
% measured_gyro = a.measured_gyro(t, a.o_d_n_g_c, a.biasStabDistGyro, a.biasTempDistGyro, a.gyro_drift, a.noiseDistGyro, a.real_ang_rate);

%% Get Error Signal

fprintf('\nIntegrating Acceleration (no Sensor Error)')
fprintf("\nTime: ")

% Initial Angle
p_vel0 = a.real_vel(tspan(1)); % Platform inertial angle [deg]

[t_base, vel_base]= rk4_solver(a.real_accel, tspan, p_vel0, dt);

% Integration Verification Plot
figure
plot(t_base, vel_base)
hold on
plot(t_base, a.real_vel(t_base))
xlabel('Time (sec)')
ylabel('Velocity (m/s)')
legend('Integrated Velocity', 'Velocity Equation')
title('Basic Integrated Velocity Plot')

fprintf('\nFinished Integration with NO Sensor Error.\n')
fprintf('\nStarting Integration WITH Sensor Error.')
fprintf("\nTime: ")

% Running simulation with sensor error

accel_error_signal = zeros(length(t_base), 2);

for i = 1:length(t_base)
    time_i = t_base(i);
    accel_i = a.real_accel(time_i);

    accel_error_signal(i, 1) = a.measured_accel_vert(time_i, a.o_d_n_a_c_v, a.biasStabDistAccel, a.biasTempDistAccel, a.accel_drift_vert, a.noiseDistAccel, accel_i, a.theta_err);
    accel_error_signal(i, 2) = a.measured_accel_horz(time_i, a.o_d_n_a_c_h, a.biasStabDistAccel, a.biasTempDistAccel, a.accel_drift_horz, a.noiseDistAccel, accel_i, a.theta_err);
end

% Basic gyro signal (real angular rate with error)
figure
plot(t_base, a.real_accel(t_base))
hold on
plot(t_base, accel_error_signal(:, 1))
hold on
plot(t_base, accel_error_signal(:, 2))
xlabel('Time (sec)')
ylabel('Acceleration (m/s^2)')
legend('Real Value', 'Accel Measured Signal (vertical)', 'Accel Measured Signal (horizontal)')
title('Acceleration Real vs Measured')

control_dynamics = @(t_i, state) [a.measured_accel_vert(time_i, a.o_d_n_a_c_v, a.biasStabDistAccel, a.biasTempDistAccel, a.accel_drift_vert, a.noiseDistAccel, a.real_accel(t_i), a.theta_err); 
    a.measured_accel_horz(time_i, a.o_d_n_a_c_h, a.biasStabDistAccel, a.biasTempDistAccel, a.accel_drift_horz, a.noiseDistAccel, a.real_accel(t_i), a.theta_err)];

[t_error, vel_error]= rk4_solver(control_dynamics, tspan, [p_vel0 p_vel0], dt);

fprintf('\nFinished Integration WITH Sensor Error.\n')

%% Plotting Platform Position

figure
plot(t_base, a.real_vel(t_base))
hold on
plot(t_error, vel_error(:, 1))
hold on
plot(t_error, vel_error(:, 2))
xlabel('Time (sec)')
ylabel('Velocity (m/s)')
legend('NO Sensor Error', 'WITH Sensor Error (vertical)', 'WITH Sensor Error (horizontal)')
title('Basic Integrated Velocity Plot')


%% Getting and Plotting error

pos_err = vel_error(:, 1) - vel_base(:, 1);

sz = 2;
plot_scale = 10;

figure
scatter(t_error, pos_err, sz, 'filled', displayName="Velocity Error")
title('Integrated (vertical) Velocity Error vs Time')
xlabel('Time (s)')
ylabel('Error (m/s)')


%% Error Compensation

corrected_accel_signal = zeros(length(t_base), 2);

for i = 1:length(t_base)
    time_i = t_base(i);
    measured_accel = accel_error_signal(i, :);
    measuredState = [measured_accel(2) measured_accel(2) measured_accel(1) 0 0 0];
    corrected_state = compensateError(measuredState, specs, time_i);

    corrected_accel_signal(i, :) = [corrected_state(1) corrected_state(3)];
end


% Basic gyro signal (real angular rate with error)
figure
plot(t_base, accel_error_signal(:, 1))
hold on
plot(t_base, corrected_accel_signal(:, 1))
hold on
plot(t_base, corrected_accel_signal(:, 2))
hold on
plot(t_base, a.real_accel(t_base))
xlabel('Time (sec)')
ylabel('Acceleration (m/s^2)')
legend('Accel Measured Signal', 'Corrected Accel Signal (vertical)', 'Corrected Accel Signal (horizontal)', 'Real Value')
title('Acceleration Real vs Measured')


% Integrate corrected gyro signal

control_dynamics = @(t_i, state) corrected_accel_signal((floor(t_i./dt)+1), :)';

[t_corr, vel_corr]= rk4_solver(control_dynamics, tspan, [p_vel0 p_vel0], dt);

figure
plot(t_base, vel_base)
hold on
plot(t_error, vel_error(:, 1))
hold on
plot(t_corr, vel_corr(:, 1))
xlabel('Time (sec)')
ylabel('Velocity (m/s)')
legend('NO Sensor Error', 'WITH Sensor Error', 'Corrected Accel Signal')
title('Basic Integrated (Vertical) Velocity Plot')
% ylim([-5 5])

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
plot(t_base, vel_base)
hold on
plot(t_error, vel_error(:, 1))
hold on
plot(t_corr, vel_corr(:, 3))
xlabel('Time (sec)')
ylabel('Velocity (m/s)')
legend('NO Sensor Error', 'WITH Sensor Error', 'Corrected (vertical) Accel Signal')
title('Basic Integrated Velocity Plot')
% ylim([-5 5])

function s_dot = AccelDriftCorrection(time_i, a, prev_state)

    % Current states
    vel_x = prev_state(1);
    vel_y = prev_state(2);
    vel_z = prev_state(3);
    vel_err_accum_x = prev_state(4);
    vel_err_accum_y = prev_state(5);
    vel_err_accum_z = prev_state(6);

    vel = [vel_x; vel_y; vel_z];
    vel_err_accum = [vel_err_accum_x; vel_err_accum_y; vel_err_accum_z];

    specs = a.specs;

    % gyro_noise_std = specs.gyro_noiseDensity * sqrt(specs.gyro_bandwidth);
    % vel_0 = normrnd(0, gyro_noise_std);

    kw = a.kw; % 
    kt = a.kt; % 

    % kw = 0.5; % 
    % kt = 0.5; % 

    vel_0 = 0;
    vel_dot_0 = 0;

    accel_i = a.real_accel(time_i);
    vel_dot_m_v = a.measured_accel_vert(time_i, a.o_d_n_a_c_v, a.biasStabDistAccel, a.biasTempDistAccel, a.accel_drift_vert, a.noiseDistAccel, accel_i, a.theta_err);
    vel_dot_m_h = a.measured_accel_horz(time_i, a.o_d_n_a_c_h, a.biasStabDistAccel, a.biasTempDistAccel, a.accel_drift_horz, a.noiseDistAccel, accel_i, a.theta_err);
    vel_dot_m = [vel_dot_m_h; vel_dot_m_h; vel_dot_m_v];

    vel_control = vel - vel_0;

    vel_dot_comp = kt * vel_err_accum + kw * vel_control;

    vel_dot = vel_dot_m + vel_dot_0 - vel_dot_comp;

    s_dot = zeros(2,1);

    s_dot(1:3) = vel_dot;
    s_dot(4:6) = vel_control;

end