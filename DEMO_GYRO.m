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

% Testing simple equations to verify integration

real_vel = @(t, y) 1/3*t.^3;
real_ang = @(t, y) 1/6*t.^4; % [rad]

real_accel = @(t, y) t.^2; % [m*s^-2]
real_ang_rate = @(t, y) 2/3*t.^3; % [rad/s]

dynamics = @(t, y) [ t.^2; 2/3*t.^3 ]; % [ real_accel real_ang_rate]

%% Sensor Model Aspects

imx_5_specs
% test3

% specs.k = 0;

% ----------- Set gyro parameters here for testing if desired -------------

% specs.gyro_resolution = 0.0076 / 180*pi; % rad/s
% specs.gyro_resolution = 1*10^-30; % rad/s
% specs.gyro_noiseDensity = 0 * (10^-3) / 180*pi; % rad/s/sqrt(Hz)
% specs.gyro_temp_bias = 0 / 180*pi; % rad/s RMS

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
plot(t, a.noisyGyroCurve(t, a.real_ang_rate));
hold on
plot(t, a.real_ang_rate(t))
xlabel('Time (sec)')
ylabel('Angular Rate (rad/s)')
title('Gyroscope Noise Effect')
legend('Measured Gyroscope Signal', 'Real Angular Rate')

% Test integration with just this source of noise

% Integrate gyro noise
noise_handle = @(t, y) a.noisyGyroCurve(t, a.real_ang_rate);
[t_arw, ARW_curve]= rk4_solver(noise_handle, tspan, a.real_ang(tspan(1)), dt);

% Integration Verification Plot
figure
plot(t_arw, ARW_curve)
hold on
plot(t_arw, a.real_ang(t_arw))
xlabel('Time (sec)')
ylabel('Angle (rad)')
title('Integrated Noisy Gyro Curve - Angle Random Walk')
legend('Integrated with Noise', 'Truth Angle')

%% Adding bias instability

% But first, turn-on bias
accel_turn_on_bias_offset = 0.01; % m/s^2
gyro_turn_on_bias_offset = 0.01; % rad/s

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
plot(t, a.biasStabDistGyro(t) + a.biasTempDistGyro(t));
xlabel('Time (sec)')
ylabel('Angular Rate (rad/s)')
title('Gyroscope Drift (Bias instability and Temperature Drift)')

a.theta_err = @(t) a.biasStabDistGyro(t).*t + ARW.*sqrt(t);

a.accel_drift_vert = @(t, accel, theta_err, biasStabDistAccel, biasTempDistAccel) (1 + k)*accel + biasStabDistAccel(t)  + biasTempDistAccel(t) + g*(1-cos(theta_err(t))) + accel_turn_on_bias_offset;
a.accel_drift_horz = @(t, accel, theta_err, biasStabDistAccel, biasTempDistAccel) (1 + k)*accel + biasStabDistAccel(t)  + biasTempDistAccel(t)  + g*sin(theta_err(t)) + accel_turn_on_bias_offset;

a.gyro_drift = @(t, ang_rate, biasStabDistGyro, biasTempDistGyro) (1 + k)*ang_rate + a.biasStabDistGyro(t) + biasTempDistGyro(t) + gyro_turn_on_bias_offset;

figure
plot(t, a.gyro_drift(t, a.real_ang_rate(t), a.biasStabDistGyro, a.biasTempDistGyro));
hold on
plot(t, a.real_ang_rate(t))
xlabel('Time (sec)')
ylabel('Angular Rate (rad/s)')
title('Gyroscope Drifting Signal')
legend('Real Angular Rate', 'Measured Gyroscope Signal')

% Test integration with just this source of noise

% Integrate gyro noise
bias_handle = @(t, y) a.biasedGyroCurve(t, a.real_ang_rate, a.biasStabDistGyro, a.biasTempDistGyro);
[t_bias, biased_curve]= rk4_solver(bias_handle, tspan, a.real_ang(tspan(1)), dt);

% Integration Verification Plot
figure
plot(t_bias, 180/pi*biased_curve)
hold on
plot(t_bias, 180/pi*a.real_ang(t_bias))
ylim([-30 30])
xlabel('Time (sec)')
ylabel('Angle (deg)')
title('Integrated Biased Gyro Curve')
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

fprintf('\nIntegrating Angular Rate (no Sensor Error)')
fprintf("\nTime: ")

% Initial Angle
p_theta0 = a.real_ang(tspan(1)); % Platform inertial angle [deg]

[t_base, theta_base]= rk4_solver(a.real_ang_rate, tspan, p_theta0, dt);

% Integration Verification Plot
figure
plot(t_base, 180/pi*theta_base)
hold on
plot(t_base, 180/pi*a.real_ang(t_base))
xlabel('Time (sec)')
ylabel('Angle (deg)')
legend('Integrated Angle', 'Angle Equation')
title('Basic Integrated Angle Plot')

fprintf('\nFinished Integration with NO Sensor Error.\n')
fprintf('\nStarting Integration WITH Sensor Error.')
fprintf("\nTime: ")

% Running simulation with sensor error

gyro_error_signal = zeros(length(t_base), 1);

for i = 1:length(t_base)
    time_i = t_base(i);
    ang_rate_i = a.real_ang_rate(time_i);
    gyro_error_signal(i) = a.measured_gyro(time_i, a.o_d_n_g_c, a.biasStabDistGyro, a.biasTempDistGyro, a.gyro_drift, a.noiseDistGyro, ang_rate_i);
end

% Basic gyro signal (real angular rate with error)
figure
plot(t_base, a.real_ang_rate(t_base))
hold on
plot(t_base, gyro_error_signal)
% hold on
xlabel('Time (sec)')
ylabel('Angular Rate (rad/s)')
legend('Real Value', 'Gyro Measured Signal')
title('Angular Rate Real vs Measured')

control_dynamics = @(t_i, state) a.measured_gyro(t_i, a.o_d_n_g_c, a.biasStabDistGyro, a.biasTempDistGyro, a.gyro_drift, a.noiseDistGyro, a.real_ang_rate(t_i));

[t_error, theta_error]= rk4_solver(control_dynamics, tspan, p_theta0, dt);

fprintf('\nFinished Integration WITH Sensor Error.\n')

%% Plotting Platform Position

figure
plot(t_base, 180/pi*theta_base)
hold on
plot(t_error, 180/pi*theta_error)
xlabel('Time (sec)')
ylabel('Angle (deg)')
legend('NO Sensor Error', 'WITH Sensor Error')
title('Basic Integrated Angle Plot')
ylim([-30 30])

%% Getting and Plotting error

pos_err = theta_error(:, 1) - theta_base(:, 1);

sz = 2;
plot_scale = 10;

figure
scatter(t_error, pos_err*180/pi, sz, 'filled', displayName="Angle Error")
title('Integrated Angle Error vs Time')
xlabel('Time (s)')
ylabel('Error (deg)')


%% Error Compensation

corrected_gyro_signal = zeros(size(t_base));

for i = 1:length(t_base)
    time_i = t_base(i);
    measured_gyro = gyro_error_signal(i);
    measuredState = [0 0 0 measured_gyro measured_gyro measured_gyro];
    corrected_state = compensateError(measuredState, specs, time_i);

    corrected_gyro_signal(i) = corrected_state(4);
end


% Basic gyro signal (real angular rate with error)
figure
plot(t_base, gyro_error_signal)
hold on
plot(t_base, corrected_gyro_signal)
hold on
plot(t_base, a.real_ang_rate(t_base))
xlabel('Time (sec)')
ylabel('Angular Rate (rad/s)')
legend('Gyro Measured Signal', 'Corrected Gyro Signal', 'Real Value')
title('Angular Rate Real vs Measured')


% Integrate corrected gyro signal

control_dynamics = @(t_i, state) corrected_gyro_signal(floor(t_i./dt)+1);

[t_corr, theta_corr]= rk4_solver(control_dynamics, tspan, p_theta0, dt);

figure
plot(t_base, 180/pi*theta_base)
hold on
plot(t_error, 180/pi*theta_error)
hold on
plot(t_corr, 180/pi*theta_corr)
xlabel('Time (sec)')
ylabel('Angle (deg)')
legend('NO Sensor Error', 'WITH Sensor Error', 'Corrected Gyro Signal')
title('Basic Integrated Angle Plot')
ylim([-30 30])


%% Control Law Drift Compensation


fprintf('\nStarting Integration WITH Sensor Compensation Control Law.')
fprintf("\nTime: ")

control_dynamics = @(t_i, state) GyroDriftCorrection(t_i, a, state);

theta_0 = 0;
theta_dot_0 = 0;

s0 = [theta_0 theta_dot_0];

[t_corr, theta_corr]= rk4_solver(control_dynamics, tspan, s0, dt);

fprintf('\nFinished Integration WITH Sensor Compensation Control Law.\n')

figure
plot(t_base, 180/pi*theta_base)
hold on
plot(t_error, 180/pi*theta_error)
hold on
plot(t_corr, 180/pi*theta_corr(:, 1))
xlabel('Time (sec)')
ylabel('Angle (deg)')
legend('NO Sensor Error', 'WITH Sensor Error', 'Corrected Gyro Signal')
title('Basic Integrated Angle Plot')
ylim([-30 30])

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
    theta_dot_m = a.measured_gyro(time_i, a.o_d_n_g_c, a.biasStabDistGyro, a.biasTempDistGyro, a.gyro_drift, a.noiseDistGyro, ang_rate_i);

    theta_control = theta - theta_0;

    theta_dot_comp = kt * theta_err_accum + kw * theta_control;

    theta_dot = theta_dot_m + theta_dot_0 - theta_dot_comp;

    s_dot = zeros(2,1);

    s_dot(1) = theta_dot;
    s_dot(2) = theta_control;

end