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
a.real_ang = @(t) atan(beta*alpha*cos(beta*t)); % [deg]
a.real_ang_rate = @(t) 180/pi*(-(alpha*beta^2*sin(beta*t))./(alpha^2*beta^2*(cos(beta*t).^2)+1)); % [deg/s]

% Testing simple equations to verify integration

% real_vel = @(t) 1/3*t.^3;
% real_ang = @(t) 1/6*t.^4; % [deg]
% 
% real_accel = @(t) t.^2; % [m*s^-2]
% real_ang_rate = @(t) 2/3*t.^3; % [deg/s]

%% Sensor Specs and Error

% ARW/VRW is our noise measurement
% Bias is our drift measurement, note:
%   measurement of how the bias will drift during operation over time at a constant temperature
% Resolution doesn't matter, as per this:
%   https://www.vectornav.com/resources/inertial-navigation-primer/specifications--and--error-budgets/specs-imuspecs

% Sensor
%   IMX-5: https://docs.inertialsense.com/datasheets/IMX-5_IMU_AHRS_GNSS-INS_Datasheet.pdf

% Simulation time
startTime = 0;
finishTime = 180; %
tspan = [startTime finishTime]; % [s]
dt = 0.01;  % [s]
t = (tspan(1):dt:tspan(2))';

imx_5_specs
% gx5_specs

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

%% Standard Noise Visualization

accel_noise_std = specs.accel_noiseDensity * sqrt(specs.accel_bandwidth); % Noise standard deviation (m/s^2)
gyro_noise_std = specs.gyro_noiseDensity * sqrt(specs.gyro_bandwidth); % Noise standard deviation (dps)

% Plotting noise distributions

a.noiseDistAccel = @(t) accel_noise_std*randn(length(t),1);
a.noiseDistGyro = @(t) gyro_noise_std*randn(length(t),1);

figure
subplot(1,2,1)
plot(t, a.noiseDistAccel(t)/a.g*1000)
xlabel('Time (s)')
ylabel('Noise Offset (mg)')
title('Accelerometer Noise Offset Distribution')
subplot(1,2,2)
xlabel('Time (s)')
plot(t, a.noiseDistGyro(t))
ylabel('Noise Offset (dps)')
title('Gyroscope Noise Offset Distribution')

% Noisy signal
a.noisyAccelCurve = @(t, real_accel) real_accel(t) + accel_noise_std*randn(length(t),1);
a.noisyGyroCurve = @(t, real_ang_rate) real_ang_rate(t) + gyro_noise_std*randn(length(t),1);

figure
subplot(1,2,1)
plot(t, a.noisyAccelCurve(t, a.real_accel)/a.g)
xlabel('Time (s)')
ylabel('Noise Offset (g)')
title('Accelerometer Noise Offset Distribution')
subplot(1,2,2)
xlabel('Time (s)')
plot(t, a.noisyGyroCurve(t, a.real_ang_rate))
ylabel('Noise Offset (dps)')
title('Gyroscope Noise Offset Distribution')

% Columns represent each variable, so [accel gyro]
% Rows represents timesteps, max limit of 5 timesteps kept track of
integrator_type = 2;
state_record = zeros(5, 2);
deriv_record = zeros(5, 2);
state_int = zeros(length(t), 2);

for i = 1:length(t)

    time = t(i);

    val_dot = [a.real_accel(time) a.real_ang_rate(time)];
    deriv_record = insertVector(deriv_record, val_dot);

    state_vec = fdm_integrator(state_record, deriv_record, dt, integrator_type);

    state_record = insertVector(state_record, state_vec);

    state_int(i, :) = state_vec;
end

vel = state_int(:, 1);
angle = state_int(:, 2);

figure
subplot(2,1,1)
plot(t, a.real_accel(t));
xlabel('Time (sec)')
ylabel('Linear Acceleration (m/s^s)')
title('Wave Linear Acceleration over Time')

subplot(2,1,2)
plot(t, vel);
xlabel('Time (sec)')
ylabel('Linear Velocity (m/s^2)')
title('Wave Integrated Velocity (from acceleration)')

figure
subplot(2,1,1)
plot(t, a.real_ang_rate(t));
xlabel('Time (sec)')
ylabel('Angular Velocity (rad/s)')
title('Wave Angular Velocity over Time')

subplot(2,1,2)
plot(t, angle);
xlabel('Time (sec)')
ylabel('Angle (deg)')
title('Wave Integrated Angle (from angular velocity)')

% Integrate Noisy Signal to show Random Walk

% Columns represent each variable, so [accel gyro]
% Rows represents timesteps, max limit of 5 timesteps kept track of
integrator_type = 2;
state_record_noise = zeros(5, 2);
deriv_record_noise = zeros(5, 2);
state_int_noise = zeros(length(t), 2);

for i = 1:length(t)

    time = t(i);

    val_dot_noise = [a.noisyAccelCurve(time, a.real_accel) a.noisyGyroCurve(time, a.real_ang_rate)];
    deriv_record_noise = insertVector(deriv_record_noise, val_dot_noise);

    state_vec_noise = fdm_integrator(state_record_noise, deriv_record_noise, dt, integrator_type);

    state_record_noise = insertVector(state_record_noise, state_vec_noise);

    state_int_noise(i, :) = state_vec_noise;
end

vel_noise = state_int_noise(:, 1);
angle_noise = state_int_noise(:, 2);

figure
subplot(2,1,1)
plot(t, a.noisyAccelCurve(t, a.real_accel))
hold on 
plot(t, a.real_accel(t));
xlabel('Time (s)')
ylabel('Linear Acceleration with Noise (m/s^2)')
title('Noisy Linear Acceleration Curve')

subplot(2,1,2)
plot(t, vel_noise);
hold on
plot(t, vel);
xlabel('Time (sec)')
ylabel('Linear Velocity (m/s)')
title('Integrated Velocity with Noise (from Linear Acceleration)')

figure
subplot(2,1,1)
xlabel('Time (s)')
plot(t, a.noisyGyroCurve(t, a.real_ang_rate))
hold on
plot(t, a.real_ang_rate(t));
ylabel('Angular Velocity with Noise (dps)')
title('Noisy Angular Velocity Curve')

subplot(2,1,2)
plot(t, angle_noise);
hold on
plot(t, angle);
xlabel('Time (sec)')
ylabel('Angle (deg)')
title('Integrated Noisy Angle with Noise (from Angular Velocity)')

% Random Walk Difference

% n_a = @(t) 0.5*VRW*t.^(-0.5);
% n_g = @(t) 0.5*ARW*t.^(-0.5);

n_a = @(t) VRW*t.^(0.5);
n_g = @(t) ARW*t.^(0.5);

vel_noise_walk_diff = abs(vel - vel_noise);
angle_noise_walk_diff = abs(angle - angle_noise);

figure
subplot(2,1,1)
scatter(t, vel_noise_walk_diff, sz, 'filled')
hold on
plot(t, n_a(t))
xlabel('Time (sec)')
ylabel('Velocity Error (m/s)')
title('Velocity Random Walk')
legend('Noise Walk', 'Accelerometer Gaussian Noise Walk')

subplot(2,1,2)
scatter(t, angle_noise_walk_diff, sz, 'filled')
hold on
plot(t, n_g(t))
xlabel('Time (sec)')
ylabel('Angle Error (deg)')
title('Angle Random Walk')
legend('Noise Walk', 'Gyroscope Gaussian Noise Walk')

%% Showing drift for regular signals

accel_curve = a.real_accel(t);
gyro_curve = a.real_ang_rate(t);

gyro_tau = 150; % [s]
accel_tau = 20; % [s]

% driftPeriod = 5 * 60;  % drift changes every 5 minutes
driftPeriod = 2;  % drift changes every 5 seconds
driftAlterations = floor(tspan(2)/driftPeriod)+1; % drift changes every 5 minutes

% b_a_drift_vals = b_a.*ones(driftAlterations,1);
% b_g_drift_vals = b_g.*ones(driftAlterations,1);

% a.biasStabDistAccel = @(t) b_a*(1-exp(-t/accel_tau));
% a.biasStabDistGyro = @(t) b_g*(1-exp(-t/gyro_tau));

b_a_drift_vals = b_a.*rand(driftAlterations,1);
b_g_drift_vals = b_g.*rand(driftAlterations,1);

bias_indeces = floor(t./(length(t)/driftAlterations)*100)+1;

a.biasStabDistAccel = b_a_drift_vals(bias_indeces);
a.biasStabDistGyro = b_g_drift_vals(bias_indeces);

figure
subplot(2,1,1)
plot(t, a.biasStabDistAccel);
hold on
yline(b_a)
hold on
yline(-b_a)
xlabel('Time (sec)')
ylabel('Acceleration (m/s^2)')
title('Accelerometer Bias Offset Model')
subplot(2,1,2)
plot(t, a.biasStabDistGyro);
hold on
yline(b_g)
hold on
yline(-b_g)
xlabel('Time (sec)')
ylabel('Angular Rate (m/s^2)')
title('Gyroscope Bias Biass Instability Offset Model')

a.biasedAccelCurve = @(t, real_accel) real_accel(t) + b_a_drift_vals(bias_indeces);
a.biasedGyroCurve = @(t, real_ang_rate) real_ang_rate(t) + b_g_drift_vals(bias_indeces);

figure
subplot(2,1,1)
plot(t, a.biasedAccelCurve(t, a.real_accel));
xlabel('Time (sec)')
ylabel('Acceleration (m/s^2)')
title('Accelerometer Bias Offset Model')
subplot(2,1,2)
plot(t, a.biasedGyroCurve(t, a.real_ang_rate));
xlabel('Time (sec)')
ylabel('Angular Rate (m/s^2)')
title('Gyroscope Bias Offset Model')

% Bias Error over Temp

% a.biasTempDistAccel = @(t) accel_temp_bias*randn(length(t),1);
% a.biasTempDistGyro = @(t) gyro_temp_bias*randn(length(t),1);

a.biasTempDistAccel = @(t) 0;
a.biasTempDistGyro = @(t) 0;

figure
subplot(1,2,1)
plot(t, a.biasTempDistAccel(t)/a.g*1000)
xlabel('Time (s)')
ylabel('Noise Offset (mg)')
title('Accelerometer Temp Bias Offset Distribution')
subplot(1,2,2)
xlabel('Time (s)')
plot(t, a.biasTempDistGyro(t))
ylabel('Noise Offset (dps)')
title('Gyroscope Temp Bias Offset Distribution')

figure
subplot(1,2,1)
plot(t, accel_curve + a.biasTempDistAccel(t)); hold on
plot(t, accel_curve);
xlabel('Time (s)')
ylabel('Noise Offset (mg)')
title('Accelerometer Temp Bias Offset Distribution')
legend('Temp Noise Curve', 'Real Curve')
subplot(1,2,2)
plot(t, gyro_curve + a.biasTempDistGyro(t)); hold on
plot(t, gyro_curve);
xlabel('Time (s)')
ylabel('Noise Offset (dps)')
title('Gyroscope Temp Bias Offset Distribution')
legend('Temp Noise Curve', 'Real Curve')

a.theta_err = @(t) a.biasStabDistGyro.*t + ARW.*sqrt(t);

figure
plot(t, a.theta_err(t))
xlabel('Time (sec)')
ylabel('Angle (def)')
title('Angle Error Growth over Time')

a.accel_drift_vert = @(t, real_accel, theta_err, biasTempDistAccel) (1 + k)*real_accel(t) + a.biasStabDistAccel  + biasTempDistAccel(t) + g*(1-cos(theta_err(t)));
a.accel_drift_horz = @(t, real_accel, theta_err, biasTempDistAccel) (1 + k)*real_accel(t) + a.biasStabDistAccel  + biasTempDistAccel(t)  + g*sin(theta_err(t));

a.gyro_drift = @(t, real_ang_rate, biasTempDistGyro) (1 + k)*real_ang_rate(t) + a.biasStabDistGyro + biasTempDistGyro(t);

figure
subplot(2,1,1)
plot(t, accel_curve); hold on
plot(t, a.accel_drift_vert(t, a.real_accel, a.theta_err, a.biasTempDistAccel)); hold on
plot(t, a.accel_drift_horz(t, a.real_accel, a.theta_err, a.biasTempDistAccel));
xlabel('Time (s)')
ylabel('Bias Offset (dps)')
title('Accelerometer Bias Drift')
legend('Real Acceleration Signal', 'Drifting Vertical Accelerometer Signal', 'Drifting Horizontal Accelerometer Signal')

subplot(2,1,2)
plot(t, gyro_curve); hold on
plot(t, a.gyro_drift(t, a.real_ang_rate, a.biasTempDistGyro));
xlabel('Time (s)')
ylabel('Bias Offset (m/s^{2})')
title('Gyroscope Bias Drift')
legend('Real Angular Rate Signal', 'Drifting Angular Rate Signal')


%% Integrating to show bias integration

% Columns represent each variable, so [accel gyro]
% Rows represents timesteps, max limit of 5 timesteps kept track of
integrator_type = 2;
state_record_bias = zeros(5, 3);
deriv_record_bias = zeros(5, 3);
state_int_bias = zeros(length(t), 3);

for i = 1:length(t)

    time = t(i)

    val_dot_bias = [a.accel_drift_vert(t, a.real_accel, a.theta_err, a.biasTempDistAccel) a.accel_drift_horz(t, a.real_accel, a.theta_err, a.biasTempDistAccel) a.gyro_drift(t, a.real_ang_rate, a.biasTempDistGyro)];
    deriv_record_bias = insertVector(deriv_record_bias, val_dot_bias);

    state_vec_bias = fdm_integrator(state_record_bias, deriv_record_bias, dt, integrator_type);

    state_record_bias = insertVector(state_record_bias, state_vec_bias);

    state_int_bias(i, :) = state_vec_bias;
end

vert_vel_bias = state_int_bias(:, 1);
horz_vel_bias = state_int_bias(:, 2);
angle_bias = state_int_bias(:, 3);

figure
subplot(2,1,1)
plot(t, accel_curve); hold on
plot(t, a.accel_drift_vert(t, a.real_accel, a.theta_err, a.biasTempDistAccel)); hold on
plot(t, a.accel_drift_horz(t, a.real_accel, a.theta_err, a.biasTempDistAccel));
xlabel('Time (s)')
ylabel('Linear Acceleration with Bias (m/s^2)')
title('Biased Linear Acceleration Curve')
legend('Real Accel Curve', 'Accel Biased Vertical Curve', 'Accel Biased Horizontal Curve')

subplot(2,1,2)
plot(t, vert_vel_bias);
hold on
plot(t, horz_vel_bias);
hold on
plot(t, vel);
xlabel('Time (sec)')
ylabel('Linear Velocity (m/s)')
title('Integrated Velocity with Bias (from Linear Acceleration)')
legend('Vertical Velocity (integrated) Biased Curve', 'Horizontal Velocity (integrated) Biased Curve', 'Real Velocity Curve')

figure
subplot(2,1,1)
xlabel('Time (s)')
plot(t, a.real_ang_rate(t)); hold on
plot(t, a.gyro_drift(t, a.real_ang_rate, a.biasTempDistGyro))
ylabel('Angular Velocity with Bias (dps)')
title('Biased Angular Velocity Curve')
legend('Real Angular Rate Curve', 'Angular Rate Biased Curve')

subplot(2,1,2)
plot(t, angle_bias);
hold on
plot(t, angle);
xlabel('Time (sec)')
ylabel('Angle (deg)')
title('Integrated Angle with Bias (from Angular Velocity)')
legend('Angle (integrated) Biased Curve', 'Real Angle Curve')

% Drift Difference

vert_vel_bias_diff = abs(vel - vert_vel_bias);
horz_vel_bias_diff = abs(vel - horz_vel_bias);
angle_bias_diff = abs(angle - angle_bias);

figure
subplot(2,1,1)
scatter(t, vert_vel_bias_diff, sz, 'filled')
hold on
scatter(t, horz_vel_bias_diff, sz, 'filled')
xlabel('Time (sec)')
ylabel('Velocity Error (m/s)')
title('Velocity Difference due to Bias')
legend('Vertical', 'Horizontal')

subplot(2,1,2)
scatter(t, angle_bias_diff, sz, 'filled')
xlabel('Time (sec)')
ylabel('Angle Error (deg)')
title('Angle Difference due to Bias')


%% Combining bias drift, temp bias, turn-on bias, noise, and quantization

accel_turn_on_bias_offset = normrnd(0, accel_noise_std);
gyro_turn_on_bias_offset = normrnd(0, gyro_noise_std);

% o_d_n_a_c_v: offset drifting noisy accel_curve vert
% o_d_n_a_c_h: offset drifting noisy accel curve horz
% o_d_n_g_c: offset drifting noisy gyro curve

a.o_d_n_a_c_v = @(t, biasTempDistAccel, accel_drift_vert, noiseDistAccel, real_accel, theta_err) accel_drift_vert(t, real_accel, theta_err, biasTempDistAccel) + noiseDistAccel(t) + accel_turn_on_bias_offset;
a.o_d_n_a_c_h = @(t, biasTempDistAccel, accel_drift_horz, noiseDistAccel, real_accel, theta_err) accel_drift_horz(t, real_accel, theta_err, biasTempDistAccel) + noiseDistAccel(t) + accel_turn_on_bias_offset;
a.o_d_n_g_c = @(t, biasTempDistGyro, gyro_drift, noiseDistGyro, real_ang_rate) gyro_drift(t, real_ang_rate, biasTempDistGyro) + noiseDistGyro(t) + gyro_turn_on_bias_offset;

a.measured_accel_vert = @(t, o_d_n_a_c_v, biasTempDistAccel, accel_drift_vert, noiseDistAccel, real_accel, theta_err) accel_resolution*floor(o_d_n_a_c_v(t, biasTempDistAccel, accel_drift_vert, noiseDistAccel, real_accel, theta_err)/accel_resolution);
a.measured_accel_horz = @(t, o_d_n_a_c_h, biasTempDistAccel, accel_drift_horz, noiseDistAccel, real_accel, theta_err) accel_resolution*floor(o_d_n_a_c_h(t, biasTempDistAccel, accel_drift_horz, noiseDistAccel, real_accel, theta_err)/accel_resolution);
a.measured_gyro = @(t, o_d_n_g_c, biasTempDistGyro, gyro_drift, noiseDistGyro, real_ang_rate) gyro_resolution*floor(o_d_n_g_c(t, biasTempDistGyro, gyro_drift, noiseDistGyro, real_ang_rate)/gyro_resolution);

measured_accel_vert = a.measured_accel_vert(t, a.o_d_n_a_c_v, a.biasTempDistAccel, a.accel_drift_vert, a.noiseDistAccel, a.real_accel, a.theta_err);
measured_accel_horz = a.measured_accel_horz(t, a.o_d_n_a_c_h, a.biasTempDistAccel, a.accel_drift_horz, a.noiseDistAccel, a.real_accel, a.theta_err);
measured_gyro = a.measured_gyro(t, a.o_d_n_g_c, a.biasTempDistGyro, a.gyro_drift, a.noiseDistGyro, a.real_ang_rate);

figure
scatter(t, accel_curve, sz, 'filled')
hold on
plot(t, measured_accel_vert)
hold on
plot(t, measured_accel_horz)
xlabel('Time (s)')
ylabel('Acceleration (m/s^2)')
title('True Acceleration vs Realistic Accelerometer Signal')
legend('Real Acceleration', 'Measured Vertical Acceleration', 'Measured Horizontal Acceleration')

figure
scatter(t, gyro_curve, sz, 'filled')
hold on
plot(t, measured_gyro)
xlabel('Time (s)')
ylabel('Angular Velocity (dps)')
title('True Angular Velocity vs Realistic Gyroscope Signal')
legend('Real Angular Velocity', 'Measured Angular Velocity')

%% Integrating to show integration with ALL Error

% Columns represent each variable, so [accel gyro]
% Rows represents timesteps, max limit of 5 timesteps kept track of
integrator_type = 2;
state_record_error = zeros(5, 3);
deriv_record_error = zeros(5, 3);
state_int_error = zeros(length(t), 3);

for i = 1:length(t)

    time = t(i)

    val_dot_error = [measured_accel_vert(i) measured_accel_horz(i) measured_gyro(i)];
    deriv_record_error = insertVector(deriv_record_error, val_dot_error);

    state_vec_error = fdm_integrator(state_record_error, deriv_record_error, dt, integrator_type);

    state_record_error = insertVector(state_record_error, state_vec_error);

    state_int_error(i, :) = state_vec_error;
end

vert_vel_error = state_int_error(:, 1);
horz_vel_error = state_int_error(:, 2);
angle_error = state_int_error(:, 3);

figure
subplot(2,1,1)
plot(t, measured_accel_vert)
hold on 
plot(t, measured_accel_horz)
hold on 
plot(t, a.real_accel(t));
xlabel('Time (s)')
ylabel('Linear Acceleration (m/s^2)')
title('Linear Acceleration with All Error')
legend('Vertical Accel with Error', 'Horizontal Accel with Error', 'Real Accel Curve')

subplot(2,1,2)
plot(t, vert_vel_error);
hold on
plot(t, horz_vel_error);
hold on
plot(t, vel);
xlabel('Time (sec)')
ylabel('Linear Velocity (m/s)')
title('Integrated Velocity with All Error (from Linear Acceleration)')
legend('Vertical Velocity (integrated) with Error', 'Horizontal Velocity (integrated) with Error', 'Real Velocity Curve')

figure
subplot(2,1,1)
xlabel('Time (s)')
plot(t, measured_gyro)
hold on
plot(t, a.real_ang_rate(t));
xlabel('Time (sec)')
ylabel('Angular Velocity (dps)')
title('Angular Velocity with All Error')
legend('Angular Rate with Error', 'Real Angular Rate Curve')

subplot(2,1,2)
plot(t, angle_error);
hold on
plot(t, angle);
xlabel('Time (sec)')
ylabel('Angle (deg)')
title('Integrated Angle with All Error (from Angular Velocity)')
legend('Angle (integrated) with Error', 'Real Angle Curve')

% Drift Difference

vert_vel_error_walk_diff = abs(vel - vert_vel_error);
horz_vel_error_walk_diff = abs(vel - horz_vel_error);
angle_error_walk_diff = abs(angle - angle_error);

figure
subplot(2,1,1)
scatter(t, vert_vel_error_walk_diff, sz, 'filled')
hold on
scatter(t, horz_vel_error_walk_diff, sz, 'filled')
xlabel('Time (sec)')
ylabel('Velocity Error (m/s)')
title('Velocity Difference due to All Error')
legend('Vertical', 'Horizontal')

subplot(2,1,2)
scatter(t, angle_error_walk_diff, sz, 'filled')
xlabel('Time (sec)')
ylabel('Angle Error (deg)')
title('Angle Difference due to All Error')

%% Error Compensation

StatesOverTime_measured = zeros(length(t), 6);
StatesOverTime_measured(:, 1) = measured_accel_horz;
StatesOverTime_measured(:, 2) = measured_accel_horz;
StatesOverTime_measured(:, 3) = measured_accel_vert;
StatesOverTime_measured(:, 4) = measured_gyro;
StatesOverTime_measured(:, 5) = measured_gyro;
StatesOverTime_measured(:, 6) = measured_gyro;
StatesOverTime_measured(1,:) = [];

StatesOverTime_corrected = zeros(length(t)-1, 6);

specs.g = a.g;

errorValues = zeros(length(t)-1, 6);

for iter = 1:(length(t)-1)
    measuredState = StatesOverTime_measured(iter, :);
    time = t(iter);
    error_compensation = compensateError(measuredState, specs, time);
    StatesOverTime_corrected(iter, :) = error_compensation;
end

t = t(2:end);

figure

subplot(1,3,1)
plot(t, StatesOverTime_measured(:, 1))
hold on
plot(t, StatesOverTime_corrected(:, 1))
hold on
plot(t, a.real_accel(t), color='black', LineWidth=1)

xlabel('Time (s)')
ylabel('Acceleration (m/s^2)')
title('Drifting Horizontal Accelerometer Signal')
legend('Raw Measurement', 'Measurement Correction', 'Expected Calculation')


subplot(1,3,2)
plot(t, StatesOverTime_measured(:, 3))
hold on
plot(t, StatesOverTime_corrected(:, 3))
hold on
plot(t, a.real_accel(t), color='black', LineWidth=1)

xlabel('Time (s)')
ylabel('Acceleration (m/s^2)')
title('Drifting Vertical Accelerometer Signal')
legend('Raw Measurement', 'Measurement Correction', 'Expected Calculation')


subplot(1,3,3)
plot(t, StatesOverTime_measured(:, 4))
hold on
plot(t, StatesOverTime_corrected(:, 4))
hold on
plot(t, a.real_ang_rate(t), color='black', LineWidth=1)

xlabel('Time (s)')
ylabel('Angular Velocity (deg/s)')
title('Drifting Gyroscope Signal')
legend('Raw Measurement', 'Measurement Correction', 'Expected Calculation')

%% Applying Kalman Filter
StatesOverTime = [StatesOverTime_corrected(:, 1), StatesOverTime_corrected(:, 3), StatesOverTime_corrected(:, 4)];

k_a_scale = 1.2;
k_b_scale = 0.6;
noise_std = [k_a_scale*(accel_noise_std+accel_temp_bias+b_a); k_a_scale*(accel_noise_std+accel_temp_bias+b_a); k_b_scale*(gyro_noise_std+gyro_temp_bias+b_g)];
noise_scaled = k.*noise_std;
q = [dt; dt; dt];
filteredStates = KalmanFilter(t, StatesOverTime', noise_scaled, q);

filteredAccelHorz = filteredStates(1, :);
filteredAccelVert = filteredStates(2, :);
filteredGyro = filteredStates(3, :);

figure
plot(t, StatesOverTime_measured(:, 1))
hold on
plot(t, filteredAccelHorz)
hold on
plot(t, a.real_accel(t))
xlabel('Time (s)')
ylabel('Acceleration (m/s^2)')
title('Horizontal Accelerometer Measured Signal vs Filtered Signal')
legend('Measured Signal', 'Filtered Signal', 'True Signal')

figure
plot(t, StatesOverTime_measured(:, 3))
hold on
plot(t, filteredAccelVert)
hold on
plot(t, a.real_accel(t))
xlabel('Time (s)')
ylabel('Acceleration (m/s^2)')
title('Vertical Accelerometer Measured Signal vs Filtered Signal')
legend('Measured Signal', 'Filtered Signal', 'True Signal')

figure
plot(t, StatesOverTime_measured(:, 4))
hold on
plot(t, filteredGyro)
hold on
plot(t, a.real_ang_rate(t))
xlabel('Time (s)')
ylabel('Angular Velocity (dps)')
title('Gyroscope Measured Signal vs Filtered Signal')
legend('Measured Signal', 'Filtered Signal', 'True Signal')

%% Functions

function min_idx = findNearest(t_ref, t_compare)
    min_diff = 1;
    min_idx = 0;
    for idx=1:length(t_compare)
        curr_t = t_compare(idx);
        if abs(curr_t-t_ref) < min_diff
            min_diff = abs(curr_t-t_ref);
            min_idx = idx;
        end
    end
end

function states = KalmanFilter(t, signal, noise_std, q)
    n = length(t);
    
    dim = size(signal, 1);

    err_measure = noise_std;
    err_estimate = err_measure;

    states = zeros(size(signal));

    last_estimate = signal(:, 1);

    for i = 1:n
        mea = signal(:, i);
        
        E = err_measure + err_estimate;

        kalman_gain = err_estimate ./ E;
        K = diag(kalman_gain);

        curr_estimate = last_estimate + K * (mea - last_estimate);

        diff = diag(abs(last_estimate - curr_estimate));
        err_estimate = (diag(ones(dim,1)) - K)*err_estimate + diff*q;

        last_estimate = curr_estimate;
        states(:, i) = curr_estimate;
    end
end

function state = compensateError(measuredState, specs, time)
    k = specs.k;
    dk = specs.dk;
    b_a = specs.b_a;
    b_g = specs.b_g;
    ARW = specs.ARW;
    VRW = specs.VRW;
    accel_noiseDensity = specs.accel_noiseDensity;
    gyro_noiseDensity = specs.gyro_noiseDensity;
    accel_resolution = specs.accel_resolution;
    gyro_resolution = specs.gyro_resolution;
    g = specs.g;
    
    measured_accel = measuredState(1:3)';
    measured_gyro = measuredState(4:6)';

    % Scale Factor Error
    S_x = k;
    S_y = S_x;
    S_z = S_x;
    
    % Scale Factor Instability
    dS_x = dk;
    dS_y = dS_x;
    dS_z = dS_x;
    
    % Setting misalignment to zero for now.
    % Can be manually inputed for real situation.
    M_xy = 0;
    M_xz = 0;
    M_yx = 0;
    M_yz = 0;
    M_zx = 0;
    M_zy = 0;

    % G-dependent bias instability
    B_gx = 0;
    B_gy = B_gx;
    B_gz = B_gx;

    % b_a = 0;
    % b_g = 0;
    % 
    % ARW = 0;
    % VRW = 0;
    
    % Bias

    ACCEL_BIAS = [b_a b_a b_a]';
    GYRO_BIAS = [b_g b_g b_g]';

    theta_err = b_g/2*time + ARW*sqrt(time);
    % theta_err = b_g*time + ARW*sqrt(time);

    a_adjusted_x_y = measured_accel(1) - ACCEL_BIAS(1) - g*sin(theta_err);
    a_adjusted_z = measured_accel(3) - ACCEL_BIAS(3) - g*(1-cos(theta_err));
    a_adjusted = [a_adjusted_x_y; a_adjusted_x_y; a_adjusted_z];

    % Adding calibration error of around 2 bits of accuracy
    a_adjusted = a_adjusted + 2*accel_resolution*randn(1);
    
    A_FIX = inv([1+S_x+dS_x  M_xy       M_xz
                M_yx         1+S_y+dS_y M_yz
                M_zx         M_zy       1+S_z+dS_z]);

    corrected_a = A_FIX * a_adjusted;


    G_DEP_BIAS = [B_gx  0     0
                  0     B_gy  0
                  0     0     B_gz];

    g_adjusted = measured_gyro - GYRO_BIAS - G_DEP_BIAS*corrected_a;
    
    % Adding calibration error of around 2 bits of accuracy
    g_adjusted = g_adjusted + 2*gyro_resolution*randn(1);

    corrected_g = A_FIX * g_adjusted;

    state(1:3) = corrected_a';
    state(4:6) = corrected_g';
    
end

function output_vec = insertVector(originalVector, addVector)

    % Assumes the addVector is a row
    % Removes last row of originalVector
    % Addes as first row of originalVector

    output_vec = [addVector; originalVector(1:end-1, :)];

end

function state_dot = NoError_FixedInt(t, a, prev_state)

    % rigidArmControl is the EOM for the 1 DOF model of the inertially
    % stabilized platform. It uses inertial acceleration control when the
    % platform is close to the center of the operation region, and uses PID
    % control on the relative position as the platform goes closer to the
    % operational bounderies
    %
    % Inputs:   t    = current time
    %           a    = structure containing environmental constants and gain
    %                  values
    %           s    = vector of states
    %                = [pi; p_dot; p_err_accum; pm; pm_dot; pm_ddot] 
    %                  - pi: platform inertial position 
    %                  - p_dot: platform inertial velocity
    %                  - p_err_accum: platform relative position correction accumulation
    %                  - pm: measured platform relative position
    %                  - pm_dot: measured platform inertial velocity 
    %                  - pm_ddot: measured platform intertial acceleration
    % Outputs:  sdot = time derivative of input state vector
    %                = [p_dot; p_ddot; pr_err] where pdot is the inertial 
    %                  velocity of the platform, p_ddot is the inertial 
    %                  acceleration of the platform,and pr_err is the error in 
    %                  the relative position of the platform
    
    % Current states
    pi = prev_state(1);
    p_dot = prev_state(2);
    p_err_accum = prev_state(3);
    pm = prev_state(4);
    pm_dot = prev_state(5);
    pm_ddot = prev_state(6);
    
    % Error in relative position (distance to center of operation region)
    p = pm-a.d(t);
    p_err = p-a.pr_d;
    
    % For testing gains without mixing proportions
    % ka = a.ka; % Acceleration [kg]
    % kv = a.kv; % Velocity     [kg/s]
    % ks = a.ks; % Position     [kg*s^-2]
    % kp = a.kp; % Proportional [kg*s^-2]
    % kd = a.kd; % Derivative   [kg/s]
    % ki = a.ki; % Integral     [kg*s^-3]
    
    % Control gain proportions
    
    I = a.I(p); % Proportion of inertial stability control to apply

    ka = a.ka*I; % Acceleration [kg]
    kv = a.kv*I; % Velocity     [kg/s]
    ks = a.ks*I; % Position     [kg*s^-2]
    
    % k = a.K(p);     % Proportion of relative position control to apply
    k_h = a.K_h(p);     % Proportion of relative position control to apply

    kp = a.kp*k_h;     % Proportional [kg*s^-2]
    kd = a.kd*k_h;     % Derivative   [kg/s]
    ki = a.ki*k_h;     % Integral     [kg*s^-3]
    
    % Derivative of states
    state_dot = zeros(6,1);
    
    % Derivative of position
    state_dot(1) = p_dot;  % inertial velocity
    state_dot(4) = pm_dot; % measured inertial velocity
    
    % Control Law
    
    % Inertial stability control force
    c_i = a.initial_scale(t); % Initial scale of gains
    f_i = -(ka*pm_ddot + kv*pm_dot + ks*pm)*c_i;
    % Relative position control force
    f_pr = -(kp*p_err + ki*p_err_accum + kd*(pm_dot-a.d_dot(t)));
    
    % Platform EOM
    p_ddot = (f_i+f_pr) / a.m;
    
    % Derivative of velocity
    state_dot(2) = p_ddot; % Inertial acceleration
    state_dot(5) = pm_ddot; % Measured inertial acceleration
    
    % Derivative of measured inertial acceleration
    state_dot(6) = a.omega*(p_ddot - pm_ddot);
    
    % Error in relative position
    state_dot(3) = p_err;
    
    err_v=abs(p_dot-pm_dot);
    err_a=abs(p_ddot-pm_ddot);

end

function state_dot = rigidArmControl_FixedInt(t, a, prev_state)
    % rigidArmControl is the EOM for the 1 DOF model of the inertially
    % stabilized platform. It uses inertial acceleration control when the
    % platform is close to the center of the operation region, and uses PID
    % control on the relative position as the platform goes closer to the
    % operational bounderies
    %
    % Inputs:   t    = current time
    %           s    = vector of states
    %                = [p; p_dot; pr_err_accum] where p is the inertial position 
    %                  of the platform, pdot is the inertial velocity of the 
    %                  platform, and pr_err_accum is the integral of the error 
    %                  in the relative position of the platform
    %           a    = structure containing environmental constants and gain
    %                  values
    % Outputs:  sdot = time derivative of input state vector
    %                = [p_dot; p_ddot; pr_err] where pdot is the inertial 
    %                  velocity of the platform, p_ddot is the inertial 
    %                  acceleration of the platform,and pr_err is the error in 
    %                  the relative position of the platform
    
    % Current states
    p = prev_state(1);
    p_dot = prev_state(2);
    pr_err_accum = prev_state(3);
    pm = prev_state(4);
    pm_dot = prev_state(5);
    pm_ddot = prev_state(6);
    p_theta = prev_state(7);
    
    specs = a.specs;
    
    p_thetadot = 1/(pm_dot^2 + 1);
    
    
    
    % Inserting measured accel manually ---------------------------------------
    sz = 3;
    if t > specs.accel_resolution
    
        measured_a_x = a.measured_accel_horz(t, a.o_d_n_a_c_h, a.biasTempDistAccel, a.accel_drift_horz, a.noiseDistAccel, a.real_accel, a.theta_err);
        measured_a_y = measured_a_x;
        measured_a_z = a.measured_accel_vert(t, a.o_d_n_a_c_v, a.biasTempDistAccel, a.accel_drift_vert, a.noiseDistAccel, a.real_accel, a.theta_err);
    
        measured_g_x = a.measured_gyro(t, a.o_d_n_g_c, a.biasTempDistGyro, a.gyro_drift, a.noiseDistGyro, a.real_ang_rate);
        measured_g_y = measured_g_x;
        measured_g_z = measured_g_x;
    
        measuredState = [measured_a_x measured_a_y measured_a_z measured_g_x measured_g_y measured_g_z];
        corrected_a = compensateError(measuredState, specs, t);
    
        pm_ddot_error = pm_ddot - corrected_a(3);
        p_thetadot_error = p_thetadot - corrected_a(4);
    
        pm_ddot = pm_ddot + (pm_ddot_error/1);
        p_thetadot = p_thetadot + (p_thetadot_error/1);
    
        % pm_ddot = corrected_a(3);
        % p_thetadot = corrected_a(4);
    
        % scatter(a.fi1, t, a_error, "red")
        % hold on
        % 
        % scatter(a.fi1, t, pm_ddot, "blue")
        % hold on
        % 
        % scatter(a.fi1, t, corrected_a(3), "green")
        % hold on
    
    end
    
    
    % scatter(a.fi1, t, corrected_a(1), sz, "filled", color='blue')
    % hold on
    % scatter(a.fi1, t, corrected_a(3), sz, color='blue')
    % hold on
    
    % 
    % error = abs(corrected_a(3) - p_ddot);
    % scatter(a.fi1, t, error, sz, "red")
    % hold on
    % scatter(a.fi3, t, corrected_a(4), sz, "filled", color='blue')
    % hold on
    
    % ----------------------------------------------------------------------------------------------
    
    % Error in relative position (distance to center of operation region)
    pr = pm-a.d(t);
    pr_err = pr-a.pr_d;
    
    % For testing gains without mixing proportions
    % ka = a.ka; % Acceleration [kg]
    % kv = a.kv; % Velocity     [kg/s]
    % ks = a.ks; % Position     [kg*s^-2]
    % kp = a.kp; % Proportional [kg*s^-2]
    % kd = a.kd; % Derivative   [kg/s]
    % ki = a.ki; % Integral     [kg*s^-3]
    
    % Control gain proportions
    
    I = a.I(pr); % Proportion of inertial stability control to apply
    % I = 0.5;
    
    ka = a.ka*I; % Acceleration [kg]
    kv = a.kv*I; % Velocity     [kg/s]
    ks = a.ks*I; % Position     [kg*s^-2]
    
    k = a.K(pr);     % Proportion of relative position control to apply
    k_h = a.K_h(pr);
    % k_h =  0.5;
    
    kp = a.kp*k_h;     % Proportional [kg*s^-2]
    kd = a.kd*k_h;     % Derivative   [kg/s]
    ki = a.ki*k_h;     % Integral     [kg*s^-3]
    
    % Derivative of states
    state_dot = zeros(7,1);
    
    % Derivative of position
    state_dot(1) = p_dot;  % inertial velocity
    state_dot(4) = pm_dot; % measured inertial velocity
    
    % Control Law
    
    % Inertial stability control force
    c_i = a.initial_scale(t); % Initial scale of gains
    f_i = -(ka*pm_ddot + kv*pm_dot + ks*pm)*c_i;
    % Relative position control force
    f_pr = -(kp*pr_err + ki*pr_err_accum + kd*(pm_dot-a.d_dot(t)));
    
    % Platform EOM
    p_ddot = (f_i+f_pr) / a.m;
    
    % scatter(a.fi1, t, p_ddot, sz, "blue")
    % hold on
    
    % Derivative of velocity
    state_dot(2) = p_ddot; % Inertial acceleration
    state_dot(5) = pm_ddot; % Measured inertial acceleration
    
    % Derivative of measured inertial acceleration
    state_dot(6) = a.omega*(p_ddot - pm_ddot);
    
    % Error in relative position
    state_dot(3) = pr_err;
    
    % Derivative of angle
    state_dot(7) = p_thetadot; % Angular Rate
    
    err_v=abs(p_dot-pm_dot);
    err_a=abs(p_ddot-pm_ddot);

    % t

end

function s_dot = noError_ODEFunc(t, s, a)
% rigidArmControl is the EOM for the 1 DOF model of the inertially
% stabilized platform. It uses inertial acceleration control when the
% platform is close to the center of the operation region, and uses PID
% control on the relative position as the platform goes closer to the
% operational bounderies
%
% Inputs:   t    = current time
%           s    = vector of states
%                = [p; p_dot; pr_err_accum] where p is the inertial position 
%                  of the platform, pdot is the inertial velocity of the 
%                  platform, and pr_err_accum is the integral of the error 
%                  in the relative position of the platform
%           a    = structure containing environmental constants and gain
%                  values
% Outputs:  sdot = time derivative of input state vector
%                = [p_dot; p_ddot; pr_err] where pdot is the inertial 
%                  velocity of the platform, p_ddot is the inertial 
%                  acceleration of the platform,and pr_err is the error in 
%                  the relative position of the platform

% Current states
p = s(1);
p_dot = s(2);
pr_err_accum = s(3);
pm = s(4);
pm_dot = s(5);
pm_ddot = s(6);

% Error in relative position (distance to center of operation region)
pr = pm-a.d(t);
pr_err = pr-a.pr_d;

% For testing gains without mixing proportions
% ka = a.ka; % Acceleration [kg]
% kv = a.kv; % Velocity     [kg/s]
% ks = a.ks; % Position     [kg*s^-2]
% kp = a.kp; % Proportional [kg*s^-2]
% kd = a.kd; % Derivative   [kg/s]
% ki = a.ki; % Integral     [kg*s^-3]

% Control gain proportions

I = a.I(pr); % Proportion of inertial stability control to apply
ka = a.ka*I; % Acceleration [kg]
kv = a.kv*I; % Velocity     [kg/s]
ks = a.ks*I; % Position     [kg*s^-2]

k = a.K(pr);     % Proportion of relative position control to apply
k_h = a.K_h(pr);
kp = a.kp*k_h;     % Proportional [kg*s^-2]
kd = a.kd*k_h;     % Derivative   [kg/s]
ki = a.ki*k_h;     % Integral     [kg*s^-3]

% Derivative of states
s_dot = zeros(6,1);

% Derivative of position
s_dot(1) = p_dot;  % inertial velocity
s_dot(4) = pm_dot; % measured inertial velocity

% Control Law

% Inertial stability control force
c_i = a.initial_scale(t); % Initial scale of gains
f_i = -(ka*pm_ddot + kv*pm_dot + ks*pm)*c_i;
% Relative position control force
f_pr = -(kp*pr_err + ki*pr_err_accum + kd*(pm_dot-a.d_dot(t)));

% Platform EOM
p_ddot = (f_i+f_pr) / a.m;

% Derivative of velocity
s_dot(2) = p_ddot; % Inertial acceleration
s_dot(5) = pm_ddot; % Measured inertial acceleration

% Derivative of measured inertial acceleration
s_dot(6) = a.omega*(p_ddot - pm_ddot) ;

% Error in relative position
s_dot(3) = pr_err;

% s_dot

err_v=abs(p_dot-pm_dot);
err_a=abs(p_ddot-pm_ddot);

end

function s_dot = rigidArmControl_ODEFunc(t, s, a)
% rigidArmControl is the EOM for the 1 DOF model of the inertially
% stabilized platform. It uses inertial acceleration control when the
% platform is close to the center of the operation region, and uses PID
% control on the relative position as the platform goes closer to the
% operational bounderies
%
% Inputs:   t    = current time
%           s    = vector of states
%                = [p; p_dot; pr_err_accum] where p is the inertial position 
%                  of the platform, pdot is the inertial velocity of the 
%                  platform, and pr_err_accum is the integral of the error 
%                  in the relative position of the platform
%           a    = structure containing environmental constants and gain
%                  values
% Outputs:  sdot = time derivative of input state vector
%                = [p_dot; p_ddot; pr_err] where pdot is the inertial 
%                  velocity of the platform, p_ddot is the inertial 
%                  acceleration of the platform,and pr_err is the error in 
%                  the relative position of the platform

% Current states
p = s(1);
p_dot = s(2);
pr_err_accum = s(3);
pm = s(4);
pm_dot = s(5);
pm_ddot = s(6);
p_theta = s(7);

specs = a.specs;

p_thetadot = 1/(pm_dot^2 + 1);



% Inserting measured accel manually ---------------------------------------
sz = 3;
if t > specs.accel_resolution

    measured_a_x = a.measured_accel_horz(t, a.o_d_n_a_c_h, a.biasTempDistAccel, a.accel_drift_horz, a.noiseDistAccel, a.real_accel, a.theta_err);
    measured_a_y = measured_a_x;
    measured_a_z = a.measured_accel_vert(t, a.o_d_n_a_c_v, a.biasTempDistAccel, a.accel_drift_vert, a.noiseDistAccel, a.real_accel, a.theta_err);

    measured_g_x = a.measured_gyro(t, a.o_d_n_g_c, a.biasTempDistGyro, a.gyro_drift, a.noiseDistGyro, a.real_ang_rate);
    measured_g_y = measured_g_x;
    measured_g_z = measured_g_x;

    measuredState = [measured_a_x measured_a_y measured_a_z measured_g_x measured_g_y measured_g_z];
    corrected_a = compensateError(measuredState, specs, t);

    pm_ddot_error = pm_ddot - corrected_a(3);
    p_thetadot_error = p_thetadot - corrected_a(4);

    pm_ddot = pm_ddot + (pm_ddot_error/1);
    p_thetadot = p_thetadot + (p_thetadot_error/1);

    % pm_ddot = corrected_a(3);
    % p_thetadot = corrected_a(4);

    % scatter(a.fi1, t, a_error, "red")
    % hold on
    % 
    % scatter(a.fi1, t, pm_ddot, "blue")
    % hold on
    % 
    % scatter(a.fi1, t, corrected_a(3), "green")
    % hold on

end


% scatter(a.fi1, t, corrected_a(1), sz, "filled", color='blue')
% hold on
% scatter(a.fi1, t, corrected_a(3), sz, color='blue')
% hold on

% 
% error = abs(corrected_a(3) - p_ddot);
% scatter(a.fi1, t, error, sz, "red")
% hold on
% scatter(a.fi3, t, corrected_a(4), sz, "filled", color='blue')
% hold on

% ----------------------------------------------------------------------------------------------

% Error in relative position (distance to center of operation region)
pr = pm-a.d(t);
pr_err = pr-a.pr_d;

% For testing gains without mixing proportions
% ka = a.ka; % Acceleration [kg]
% kv = a.kv; % Velocity     [kg/s]
% ks = a.ks; % Position     [kg*s^-2]
% kp = a.kp; % Proportional [kg*s^-2]
% kd = a.kd; % Derivative   [kg/s]
% ki = a.ki; % Integral     [kg*s^-3]

% Control gain proportions

I = a.I(pr); % Proportion of inertial stability control to apply
% I = 0.5;

ka = a.ka*I; % Acceleration [kg]
kv = a.kv*I; % Velocity     [kg/s]
ks = a.ks*I; % Position     [kg*s^-2]

k = a.K(pr);     % Proportion of relative position control to apply
k_h = a.K_h(pr);
% k_h =  0.5;

kp = a.kp*k_h;     % Proportional [kg*s^-2]
kd = a.kd*k_h;     % Derivative   [kg/s]
ki = a.ki*k_h;     % Integral     [kg*s^-3]

% Derivative of states
s_dot = zeros(7,1);

% Derivative of position
s_dot(1) = p_dot;  % inertial velocity
s_dot(4) = pm_dot; % measured inertial velocity

% Control Law

% Inertial stability control force
c_i = a.initial_scale(t); % Initial scale of gains
f_i = -(ka*pm_ddot + kv*pm_dot + ks*pm)*c_i;
% Relative position control force
f_pr = -(kp*pr_err + ki*pr_err_accum + kd*(pm_dot-a.d_dot(t)));

% Platform EOM
p_ddot = (f_i+f_pr) / a.m;

% scatter(a.fi1, t, p_ddot, sz, "blue")
% hold on

% Derivative of velocity
s_dot(2) = p_ddot; % Inertial acceleration
s_dot(5) = pm_ddot; % Measured inertial acceleration

% Derivative of measured inertial acceleration
s_dot(6) = a.omega*(p_ddot - pm_ddot);

% Error in relative position
s_dot(3) = pr_err;

% Derivative of angle
s_dot(7) = p_thetadot; % Angular Rate

err_v=abs(p_dot-pm_dot);
err_a=abs(p_ddot-pm_ddot);

end