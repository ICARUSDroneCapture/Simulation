clc; clear; close all;

simulationParameters;

close all;

% Redefining acceleration/gyro curves for clarity

real_pos = @(t) alpha*sin(beta*t) + hdeck;
real_vel = @(t) beta*alpha*cos(beta*t);
real_accel = @(t) -beta^2*alpha*sin(beta*t); % [m*s^-2]
real_ang = @(t) atan(beta*alpha*cos(beta*t)); % [deg]
real_ang_rate = @(t) 180/pi*(-(alpha*beta^2*sin(beta*t))./(alpha^2*beta^2*(cos(beta*t).^2)+1)); % [deg/s]

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
% tspan = [0 60*10]; % [s]
tspan = [0 10]; % [s]
dt = 0.01; % [s]
t = (tspan(1):dt:tspan(2))';

imx_5_specs
% gx5_specs

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
plot(t, a.noisyAccelCurve(t, real_accel)/a.g)
xlabel('Time (s)')
ylabel('Noise Offset (g)')
title('Accelerometer Noise Offset Distribution')
subplot(1,2,2)
xlabel('Time (s)')
plot(t, a.noisyGyroCurve(t, real_ang_rate))
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

    val_dot = [real_accel(time) real_ang_rate(time)];
    deriv_record = insertVector(deriv_record, val_dot);

    state_vec = fdm_integrator(state_record, deriv_record, dt, integrator_type);
    
    state_record = insertVector(state_record, state_vec);

    state_int(i, :) = state_vec;
end

vel = state_int(:, 1);
angle = state_int(:, 2);

figure
subplot(2,1,1)
plot(t, real_accel(t));
xlabel('Time (sec)')
ylabel('Linear Velocity (m/s)')
title('Linear Acceleration')

subplot(2,1,2)
plot(t, vel);
xlabel('Time (sec)')
ylabel('Linear Acceleration (m/s^2)')
title('Integrated Angle (from angular velocity)')

figure
subplot(2,1,1)
plot(t, real_ang_rate(t));
xlabel('Time (sec)')
ylabel('Angular Velocity (rad/s)')
title('Angular Velocity')

subplot(2,1,2)
plot(t, angle);
xlabel('Time (sec)')
ylabel('Angle (deg)')
title('Integrated Angle (from angular velocity)')

%% Integrate Noisy Signal to show Random Walk

% Columns represent each variable, so [accel gyro]
% Rows represents timesteps, max limit of 5 timesteps kept track of
integrator_type = 2;
state_record_noise = zeros(5, 2);
deriv_record_noise = zeros(5, 2);
state_int_noise = zeros(length(t), 2);

for i = 1:length(t)

    time = t(i);

    val_dot_noise = [a.noisyAccelCurve(time, real_accel) a.noisyGyroCurve(time, real_ang_rate)];
    deriv_record_noise = insertVector(deriv_record_noise, val_dot_noise);

    state_vec_noise = fdm_integrator(state_record_noise, deriv_record_noise, dt, integrator_type);
    
    state_record_noise = insertVector(state_record_noise, state_vec_noise);

    state_int_noise(i, :) = state_vec_noise;
end

vel_noise = state_int_noise(:, 1);
angle_noise = state_int_noise(:, 2);

figure
subplot(2,1,1)
plot(t, a.noisyAccelCurve(t, real_accel))
hold on 
plot(t, real_accel(t));
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
plot(t, a.noisyGyroCurve(t, real_ang_rate))
hold on
plot(t, real_ang_rate(t));
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

subplot(2,1,2)
scatter(t, angle_noise_walk_diff, sz, 'filled')
hold on
plot(t, n_g(t))
xlabel('Time (sec)')
ylabel('Angle Error (deg)')
title('Angle Random Walk')

%% Showing drift for flat signals

accel_curve = real_accel(t);
gyro_curve = real_ang_rate(t);

gyro_tau = 150; % [s]
accel_tau = 20; % [s]

% driftPeriod = 5 * 60;  % drift changes every 5 minutes
driftPeriod = 30;  % drift changes every 5 seconds
driftAlterations = floor(tspan(2)/driftPeriod)+1; % drift changes every 5 minutes

% b_a_drift_vals = b_a.*randn(driftAlterations,1);
% b_g_drift_vals = b_g.*randn(driftAlterations,1);

b_a_drift_vals = b_a.*rand(driftAlterations,1);
b_g_drift_vals = b_g.*rand(driftAlterations,1);

bias_indeces = floor(t./(length(t)/driftAlterations)*100)+1;

a.biasStabDistAccel = b_a_drift_vals(bias_indeces);
a.biasStabDistGyro = b_g_drift_vals(bias_indeces);

% a.biasStabDistGyro = @(t) b_g.*randn(length(t),1);

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
title('Gyroscope Bias Offset Model')

a.biasedAccelCurve = @(t, real_accel) real_accel(t) + b_a_drift_vals(bias_indeces);
a.biasedGyroCurve = @(t, real_ang_rate) real_ang_rate(t) + b_g_drift_vals(bias_indeces);

figure
subplot(2,1,1)
plot(t, a.biasedAccelCurve(t, real_accel));
xlabel('Time (sec)')
ylabel('Acceleration (m/s^2)')
title('Accelerometer Bias Offset Model')
subplot(2,1,2)
plot(t, a.biasedGyroCurve(t, real_ang_rate));
xlabel('Time (sec)')
ylabel('Angular Rate (m/s^2)')
title('Gyroscope Bias Offset Model')

% a.biasStabDistAccel = @(t) b_a*(1-exp(-t/accel_tau));
% a.biasStabDistGyro = @(t) b_g*(1-exp(-t/gyro_tau));

% a.biasStabDistAccel = @(t) b_a*(1-exp(t));
% a.biasStabDistGyro = @(t) b_g*(1-exp(t));

% Drawing offset lines

biasTime = [a.biasStabDistAccel, a.biasStabDistGyro];

b_a_scale = 10;
b_g_scale = 10;
bias_std = [b_a_scale*b_a; b_g_scale*b_g];
bias_std_scaled = k.*bias_std;
q = [dt; dt];
filteredBias = KalmanFilter(t, biasTime', bias_std_scaled, q);

filteredAccelBias = filteredBias(1, :);
filteredGyroBias = filteredBias(2, :);

figure
subplot(2,1,1)
plot(t, a.biasStabDistAccel);
hold on
yline(b_a)
hold on
yline(-b_a)
hold on
plot(t, filteredAccelBias, 'linewidth', 3);
xlabel('Time (sec)')
ylabel('Acceleration (m/s^2)')
title('Accelerometer Bias Offset Model')
subplot(2,1,2)
plot(t, a.biasStabDistGyro);
hold on
yline(b_g)
hold on
yline(-b_g)
hold on
plot(t, filteredGyroBias,'linewidth', 3);
xlabel('Time (sec)')
ylabel('Angular Rate (m/s^2)')
title('Gyroscope Bias Offset Model')

a.theta_err = @(t) a.biasStabDistAccel.*t + ARW.*sqrt(t);

figure
plot(t, a.theta_err(t))

a.accel_drift_vert = @(t, real_accel, theta_err) (1 + k)*real_accel(t) + a.biasStabDistAccel + g*(1-cos(theta_err(t)));
a.accel_drift_horz = @(t, real_accel, theta_err) (1 + k)*real_accel(t) + a.biasStabDistAccel + g*sin(theta_err(t));

a.gyro_drift = @(t, real_ang_rate) (1 + k)*real_ang_rate(t) + a.biasStabDistGyro;

% flat = zeros(length(t), 1);
% gyro_drift_flat = zeros(size(flat));
% accel_drift_flat = zeros(size(flat));
% 
% for i=2:length(t)
%     w_k = flat(i);
% 
%     x_k_1 = gyro_drift_flat(i-1);
%     x_k = exp(-dt/gyro_tau)*x_k_1+w_k;
% 
%     gyro_drift_flat(i) = x_k;
% 
% 
%     x_k_1 = accel_drift_flat(i-1);
%     x_k = exp(-dt/accel_tau)*x_k_1+w_k;
% 
%     accel_drift_flat(i) = x_k;
% end

figure
subplot(2,1,1)
yline(b_a); hold on
plot(t, a.biasStabDistAccel)
xlabel('Time (s)')
ylabel('Bias Offset (dps)')
title('Accelerometer Bias Instability Offset')

subplot(2,1,2)
yline(b_g); hold on
plot(t, a.biasStabDistGyro)
xlabel('Time (s)')
ylabel('Bias Offset (m/s^{2})')
title('Gyroscope Bias Instability Offset')


figure
subplot(2,1,1)
plot(t, accel_curve); hold on
plot(t, a.accel_drift_vert(t, real_accel, a.theta_err)); hold on
plot(t, a.accel_drift_horz(t, real_accel, a.theta_err));
xlabel('Time (s)')
ylabel('Bias Offset (dps)')
title('Accelerometer Bias Drift')
legend('Real Acceleration Signal', 'Drifting Vertical Signal', 'Drifting Horizontal Signal')

subplot(2,1,2)
plot(t, gyro_curve); hold on
plot(t, a.gyro_drift(t, real_ang_rate));
xlabel('Time (s)')
ylabel('Bias Offset (m/s^{2})')
title('Gyroscope Bias Drift')
legend('Real Angular Rate Signal', 'Drifting Angular Rate Signal')

% Bias Error over Temp

a.biasTempDistAccel = @(t) accel_temp_bias*randn(length(t),1);
a.biasTempDistGyro = @(t) gyro_temp_bias*randn(length(t),1);

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


%% Integrating to show drift integration

a.accel_drift_no_g = @(t, real_accel) (1 + k)*real_accel(t) + a.biasStabDistAccel;
a.gyro_drift_no_g = @(t, real_ang_rate) (1 + k)*real_ang_rate(t) + a.biasStabDistGyro;

% Columns represent each variable, so [accel gyro]
% Rows represents timesteps, max limit of 5 timesteps kept track of
integrator_type = 2;
state_record_drift = zeros(5, 2);
deriv_record_drift = zeros(5, 2);
state_int_drift = zeros(length(t), 2);

for i = 1:length(t)

    time = t(i)

    val_dot_drift = [a.accel_drift_no_g(time, real_accel) a.gyro_drift_no_g(time, real_ang_rate)];
    deriv_record_drift = insertVector(deriv_record_drift, val_dot_drift);

    state_vec_drift = fdm_integrator(state_record_drift, deriv_record_drift, dt, integrator_type);
    
    state_record_drift = insertVector(state_record_drift, state_vec_drift);

    state_int_drift(i, :) = state_vec_drift;
end

vel_drift = state_int_drift(:, 1);
angle_drift = state_int_drift(:, 2);

figure
subplot(2,1,1)
plot(t, a.accel_drift_no_g(t, real_accel))
hold on 
plot(t, real_accel(t));
xlabel('Time (s)')
ylabel('Linear Acceleration with Drift (m/s^2)')
title('Drifting Linear Acceleration Curve')
legend('Accel Drifting Curve', 'Real Accel Curve')

subplot(2,1,2)
plot(t, vel_drift);
hold on
plot(t, vel);
xlabel('Time (sec)')
ylabel('Linear Velocity (m/s)')
title('Integrated Velocity with Drift (from Linear Acceleration)')
legend('Velocity (integrated) Drifting Curve', 'Real Velocity Curve')

figure
subplot(2,1,1)
xlabel('Time (s)')
plot(t, a.gyro_drift_no_g(t, real_ang_rate))
hold on
plot(t, real_ang_rate(t));
ylabel('Angular Velocity with Noise (dps)')
title('Drifting Angular Velocity Curve')
legend('Angular Rate Drifting Curve', 'Real Angular Rate Curve')

subplot(2,1,2)
plot(t, angle_drift);
hold on
plot(t, angle);
xlabel('Time (sec)')
ylabel('Angle (deg)')
title('Integrated Angle with Drift (from Angular Velocity)')
legend('Angle (integrated) Drifting Curve', 'Real Angle Curve')

% Drift Difference

vel_drift_walk_diff = abs(vel - vel_drift);
angle_drift_walk_diff = abs(angle - angle_drift);

figure
subplot(2,1,1)
scatter(t, vel_drift_walk_diff, sz, 'filled')
hold on
plot(t, n_a(t))
xlabel('Time (sec)')
ylabel('Velocity Error (m/s)')
title('Velocity Error due to Drift')

subplot(2,1,2)
scatter(t, angle_drift_walk_diff, sz, 'filled')
hold on
plot(t, n_g(t))
xlabel('Time (sec)')
ylabel('Angle Error (deg)')
title('Angle Error due to Drift')


%% Combining bias drift, temp bias, turn-on bias, noise, and quantization

accel_turn_on_bias_offset = normrnd(0, accel_noise_std);
gyro_turn_on_bias_offset = normrnd(0, gyro_noise_std);


a.accel_drift_vert = @(t, real_accel, theta_err) (1 + k)*real_accel(t) + a.biasStabDistAccel + accel_noise_std*randn(length(t),1);
a.accel_drift_horz = @(t, real_accel, theta_err) (1 + k)*real_accel(t) + a.biasStabDistAccel + accel_noise_std*randn(length(t),1);

a.gyro_drift = @(t, real_ang_rate) (1 + k)*real_ang_rate(t) + a.biasStabDistGyro + gyro_noise_std*randn(length(t),1);

% o_d_n_a_c_v: offset_drifting_noisy_accel_curve_vert
% o_d_n_a_c_h: offset_drifting_noisy_accel_curve_horz
% o_d_n_g_c: offset_drifting_noisy_gyro_curve

a.o_d_n_a_c_v = @(t, biasTempDistAccel, accel_drift_vert, noiseDistAccel, real_accel, theta_err) accel_drift_vert(t, real_accel, theta_err) + biasTempDistAccel(t) + noiseDistAccel(t) + accel_turn_on_bias_offset;
a.o_d_n_a_c_h = @(t, biasTempDistAccel, accel_drift_horz, noiseDistAccel, real_accel, theta_err) accel_drift_horz(t, real_accel, theta_err) + biasTempDistAccel(t) + noiseDistAccel(t) + accel_turn_on_bias_offset;
a.o_d_n_g_c = @(t, biasTempDistGyro, gyro_drift, noiseDistGyro, real_ang_rate) gyro_drift(t, real_ang_rate) + biasTempDistGyro(t) + noiseDistGyro(t) + gyro_turn_on_bias_offset;

a.measured_accel_vert = @(t, o_d_n_a_c_v, biasTempDistAccel, accel_drift_vert, noiseDistAccel, real_accel, theta_err) accel_resolution*floor(o_d_n_a_c_v(t, biasTempDistAccel, accel_drift_vert, noiseDistAccel, real_accel, theta_err)/accel_resolution);
a.measured_accel_horz = @(t, o_d_n_a_c_h, biasTempDistAccel, accel_drift_horz, noiseDistAccel, real_accel, theta_err) accel_resolution*floor(o_d_n_a_c_h(t, biasTempDistAccel, accel_drift_horz, noiseDistAccel, real_accel, theta_err)/accel_resolution);
a.measured_gyro = @(t, o_d_n_g_c, biasTempDistGyro, gyro_drift, noiseDistGyro, real_ang_rate) gyro_resolution*floor(o_d_n_g_c(t, biasTempDistGyro, gyro_drift, noiseDistGyro, real_ang_rate)/gyro_resolution);

measured_accel_vert = a.measured_accel_vert(t, a.o_d_n_a_c_v, a.biasTempDistAccel, a.accel_drift_vert, a.noiseDistAccel, real_accel, a.theta_err);
measured_accel_horz = a.measured_accel_horz(t, a.o_d_n_a_c_h, a.biasTempDistAccel, a.accel_drift_horz, a.noiseDistAccel, real_accel, a.theta_err);
measured_gyro = a.measured_gyro(t, a.o_d_n_g_c, a.biasTempDistGyro, a.gyro_drift, a.noiseDistGyro, real_ang_rate);

figure
scatter(t, accel_curve, sz, 'filled')
hold on
plot(t, measured_accel_vert)
hold on
plot(t, measured_accel_horz)
xlabel('Time (s)')
ylabel('Acceleration (m/s^2)')
title('Accelerometer Real vs Realistic Signal')
legend('Real Acceleration', 'Measured Vertical Acceleration', 'Measured Horizontal Acceleration')

figure
scatter(t, gyro_curve, sz, 'filled')
hold on
plot(t, measured_gyro)
xlabel('Time (s)')
ylabel('Angular Velocity (dps)')
title('Gyroscope Real vs Realistic Signal')
legend('Real Angular Velocity', 'Measured Angular Velocity')


%% Integrating to show drift integration

% Columns represent each variable, so [accel gyro]
% Rows represents timesteps, max limit of 5 timesteps kept track of
integrator_type = 2;
state_record_error = zeros(5, 2);
deriv_record_error = zeros(5, 2);
state_int_error = zeros(length(t), 2);

for i = 1:length(t)

    time = t(i)

    val_dot_error = [measured_accel_vert(i) measured_gyro(i)];
    deriv_record_error = insertVector(deriv_record_error, val_dot_error);

    state_vec_error = fdm_integrator(state_record_error, deriv_record_error, dt, integrator_type);
    
    state_record_error = insertVector(state_record_error, state_vec_error);

    state_int_error(i, :) = state_vec_error;
end

vel_error = state_int_error(:, 1);
angle_error = state_int_error(:, 2);

figure
subplot(2,1,1)
plot(t, measured_accel_vert)
hold on 
plot(t, real_accel(t));
xlabel('Time (s)')
ylabel('Linear Acceleration with All Error (m/s^2)')
title('Error Linear Acceleration Curve')
legend('Accel with Error Curve', 'Real Accel Curve')

subplot(2,1,2)
plot(t, vel_error);
hold on
plot(t, vel);
xlabel('Time (sec)')
ylabel('Linear Velocity (m/s)')
title('Integrated Velocity with All Error (from Linear Acceleration)')
legend('Velocity (integrated) with Error Curve', 'Real Velocity Curve')

figure
subplot(2,1,1)
xlabel('Time (s)')
plot(t, measured_gyro)
hold on
plot(t, real_ang_rate(t));
ylabel('Angular Velocity with All Error (dps)')
title('Error Angular Velocity Curve')
legend('Angular Rate with Error Curve', 'Real Angular Rate Curve')

subplot(2,1,2)
plot(t, angle_error);
hold on
plot(t, angle);
xlabel('Time (sec)')
ylabel('Angle (deg)')
title('Integrated Angle with All Error (from Angular Velocity)')
legend('Angle (integrated) with Error Curve', 'Real Angle Curve')

% Drift Difference

vel_drift_walk_diff = abs(vel - vel_error);
angle_drift_walk_diff = abs(angle - angle_error);

figure
subplot(2,1,1)
scatter(t, vel_drift_walk_diff, sz, 'filled')
hold on
plot(t, n_a(t))
xlabel('Time (sec)')
ylabel('Velocity Error (m/s)')
title('Velocity Error due to All Error')

subplot(2,1,2)
scatter(t, angle_drift_walk_diff, sz, 'filled')
hold on
plot(t, n_g(t))
xlabel('Time (sec)')
ylabel('Angle Error (deg)')
title('Angle Error due to All Error')


%% Applying Kalman Filter
StatesOverTime = [measured_accel_vert, measured_accel_horz, measured_gyro];

k_a_scale = 1;
k_b_scale = 1;
noise_std = [k_a_scale*(accel_noise_std+accel_temp_bias); k_a_scale*(accel_noise_std+accel_temp_bias); k_b_scale*(gyro_noise_std+gyro_temp_bias)];
noise_scaled = k.*noise_std;
q = [dt; dt; dt];
filteredStates = KalmanFilter(t, StatesOverTime', noise_scaled, q);

filteredAccelVert = filteredStates(1, :);
filteredAccelHorz = filteredStates(2, :);
filteredGyro = filteredStates(3, :);

figure
plot(t, measured_accel_vert)
hold on
plot(t, filteredAccelVert)
xlabel('Time (s)')
ylabel('Acceleration (m/s^2)')
title('Vertical Accelerometer Real vs Filtered Signal')
legend('Measured Vertical Acceleration Signal', 'Filtered Vertical Acceleration Signal')

figure
plot(t, measured_accel_horz)
hold on
plot(t, filteredAccelHorz)
xlabel('Time (s)')
ylabel('Acceleration (m/s^2)')
title('Horizontal Accelerometer Real vs Filtered Signal')
legend('Measured Horizontal Acceleration Signal', 'Filtered Horizontal Acceleration Signal')

figure
plot(t, measured_gyro)
hold on
plot(t, filteredGyro)
xlabel('Time (s)')
ylabel('Angular Velocity (dps)')
title('Gyroscope Real vs Filtered Signal')
legend('Measured Angular Velocity Signal', 'Filtered Angular Velocity Signal')

%% Error Compensation

g_sub_accel_horz = measured_accel_horz + g.*(1-cos(angle_error));
g_sub_accel_vert = measured_accel_vert + g.*sin(angle_error);

figure
plot(t, real_accel(t))
hold on
plot(t, g_sub_accel_horz)
hold on
plot(t, g_sub_accel_vert)
xlabel('Time (s)')
ylabel('Acceleration (m/s^2)')
title('G-removed Accelerometer Signal vs Real Signal')
legend('Real Acceleration Signal', 'Horizontal Acceleration Signal', 'Vertical Acceleration Signal')



StatesOverTime_measured = zeros(length(t), 6);
StatesOverTime_measured(:, 1) = g_sub_accel_horz;
StatesOverTime_measured(:, 2) = g_sub_accel_horz;
StatesOverTime_measured(:, 3) = g_sub_accel_vert;
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
    error_compensation = compensateError(measuredState, specs, time, angle_error(iter));
    StatesOverTime_corrected(iter, :) = error_compensation;
end

t = t(2:end);

figure

subplot(1,3,1)
plot(t, StatesOverTime_measured(:, 1))
hold on
plot(t, StatesOverTime_corrected(:, 1))
hold on
plot(t, real_accel(t), color='black', LineWidth=1)

xlabel('Time (s)')
ylabel('Acceleration (m/s^2)')
title('Drifting Horizontal Accelerometer Signal')
legend('Raw Measurement', 'Measurement Correction', 'Expected Calculation')


subplot(1,3,2)
plot(t, StatesOverTime_measured(:, 3))
hold on
plot(t, StatesOverTime_corrected(:, 3))
hold on
plot(t, real_accel(t), color='black', LineWidth=1)

xlabel('Time (s)')
ylabel('Acceleration (m/s^2)')
title('Drifting Vertical Accelerometer Signal')
legend('Raw Measurement', 'Measurement Correction', 'Expected Calculation')


subplot(1,3,3)
plot(t, StatesOverTime_measured(:, 4))
hold on
plot(t, StatesOverTime_corrected(:, 4))
hold on
plot(t, real_ang_rate(t), color='black', LineWidth=1)

xlabel('Time (s)')
ylabel('Angular Velocity (deg/s)')
title('Drifting Gyroscope Signal')
legend('Raw Measurement', 'Measurement Correction', 'Expected Calculation')

%% Functions

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

function state = compensateError(measuredState, specs, time, theta_err)
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

    % theta_err = b_g/2*time + ARW*sqrt(time);

    a_adjusted_x_y = measured_accel(1) - ACCEL_BIAS(1) - g*sin(theta_err);
    a_adjusted_z = measured_accel(3) - ACCEL_BIAS(3) - g*(1-cos(theta_err));
    a_adjusted = [a_adjusted_x_y; a_adjusted_x_y; a_adjusted_z];

    % Adding calibration error of up to 2 bits of accuracy
    a_adjusted = a_adjusted - 2*accel_resolution;
    
    A_FIX = inv([1+S_x+dS_x  M_xy       M_xz
                M_yx         1+S_y+dS_y M_yz
                M_zx         M_zy       1+S_z+dS_z]);

    corrected_a = A_FIX * a_adjusted;


    G_DEP_BIAS = [B_gx  0     0
                  0     B_gy  0
                  0     0     B_gz];

    g_adjusted = measured_gyro - GYRO_BIAS - G_DEP_BIAS*corrected_a;
    
    % Adding calibration error of up to 2 bits of accuracy
    g_adjusted = g_adjusted - 2*gyro_resolution;

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