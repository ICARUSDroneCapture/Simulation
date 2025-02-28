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

% Testing simple equations to verify integration

real_vel = @(t, y) 1/3*t.^3;
real_ang = @(t, y) 1/6*t.^4; % [rad]

real_accel = @(t, y) t.^2; % [m*s^-2]
real_ang_rate = @(t, y) 2/3*t.^3; % [rad/s]

dynamics = @(t, y) [ t.^2; 2/3*t.^3 ]; % [ real_accel real_ang_rate]

%% Sensor Model Aspects

% imx_5_specs
imx_5_specs_deg

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

% Test integration with just this source of noise

% Integrate gyro noise
noise_handle = @(t, y) a.noisyAccelCurve(t, a.real_accel);
[t_vrw, VRW_curve]= rk4_solver(noise_handle, tspan, a.real_vel(tspan(1)), dt);

%% Adding bias instability

% But first, turn-on bias
accel_turn_on_bias_offset = 0.1; % m/s^2
gyro_turn_on_bias_offset = 0.1; % rad

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

corrected_signals = zeros(length(t_count), 6);

for i = 1:length(t)
    time_i = t(i);

    ang_rate_i = a.real_ang_rate(time_i);
    accel_i = a.real_accel(time_i);
    

    measured_gyro = a.measured_gyro(time_i, a.o_d_n_g_c, a.biasStabDistGyro, a.biasTempDistGyro, a.gyro_drift, a.noiseDistGyro, ang_rate_i);
    measured_accel_v = a.measured_accel_vert(time_i, a.o_d_n_a_c_v, a.biasStabDistAccel, a.biasTempDistAccel, a.accel_drift_vert, a.noiseDistAccel, accel_i, a.theta_err);
    measured_accel_h = a.measured_accel_horz(time_i, a.o_d_n_a_c_h, a.biasStabDistAccel, a.biasTempDistAccel, a.accel_drift_horz, a.noiseDistAccel, accel_i, a.theta_err);

    measuredState = [measured_accel_h measured_accel_h measured_accel_v measured_gyro measured_gyro measured_gyro];
    corrected_state = compensateError(measuredState, specs, time_i);

    corrected_signals(i, :) = corrected_state;
end

% Initial States
p0 =  a.d(tspan(1))+a.pr_d;   % Platform position [m]
p_dot0 = a.d_dot(tspan(1));   % Platform velocity [m/s]
pr_err_accum0 = 0;            % Integral of relative position error [m*s]
pm0 = p0;                     % Platform inetegrated position [m]
pm_dot = p_dot0;              % Platform integrated velocity [m/s]
pm_ddot = a.d_ddot(tspan(1)); % Platform measured acceleration [m*s^-2]
p_theta0 = a.real_ang(tspan(1)); % Platform inertial angle [deg]

s0 = [p0; p_dot0; pr_err_accum0; pm0; pm_dot; pm_ddot];

% Integrate corrected gyro signal

control_dynamics = @(t_i, state) corrected_signals((floor(t_i./dt)+1), :)';

[t_corr, int_sig_corr]= rk4_solver(control_dynamics, tspan, s0, dt);

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

control_dynamics = @(t_i, state) DriftCorrection(t_i, a, state);

vel0_x = 0;
vel0_y = 0;
vel0_z = 0;
p0_x = 0;
p0_y = 0;
p0_z = 0;
ang0_x = 0;
ang0_y = 0;
ang0_z = 0;
a_err0_x = 0;
a_err0_y = 0;
a_err0_z = 0;

s0 = [vel0_x vel0_y vel0_z p0_x p0_y p0_z ang0_x ang0_y ang0_z a_err0_x a_err0_y a_err0_z];

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


%% Control Law Drift Compensation - 1D

fprintf('\nStarting Integration WITH Sensor Compensation Control Law.')
fprintf("\nTime: ")

control_dynamics = @(t_i, state) DriftCorrection1D(t_i, a, state);

vel0 = 0;                    % Platform velocity [m/s]
ang0 = 0; % Platform inertial angle [deg]
p0 = 0;   % Platform position [m]
a_err0 = 0;

s0 = [vel0 ang0 p0 a_err0];

[t_corr, int_state_corr]= rk4_solver(control_dynamics, tspan, s0, dt);

fprintf('\nFinished Integration WITH Sensor Compensation Control Law.\n')

figure
plot(t, a.real_vel(t))
hold on
plot(t_corr, int_state_corr(:, 1))
xlabel('Time (sec)')
ylabel('Velocity (m/s)')
legend('NO Sensor Error', 'Controlled Error (vertical)')
title('Controlled Integrated Velocity')
% ylim([-5 5])


figure
plot(t, a.real_ang(t))
hold on
plot(t_corr, int_state_corr(:, 2))
xlabel('Time (sec)')
ylabel('Angle (deg)')
legend('NO Sensor Error', 'Controlled Error')
title('Controlled Integrated Angle')
% ylim([-5 5])

figure
plot(t, int_state_corr(:, 3))
xlabel('Time (sec)')
ylabel('Position (m)')
legend('Controlled Error')
title('Platform Position (from integration)')

figure
plot(t, a.real_vel(t))
hold on
plot(t, int_signals_error(:, 3))
hold on
plot(t_corr, int_state_corr(:, 1))
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
plot(t_corr, int_state_corr(:, 2))
xlabel('Time (sec)')
ylabel('Angle (deg)')
legend('NO Sensor Error', 'Integrated Raw Signal', 'Controlled Error Integrated')
title('Controlled Integrated Angle')
% ylim([-5 5])

