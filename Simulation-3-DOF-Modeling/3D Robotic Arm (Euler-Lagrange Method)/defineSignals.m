% Script to define sensor specs and create function handles for inserting
% error into deck disturbance signals

imx_5_specs
% imx_5_specs_deg

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

imu_accel_sampling_rate = 3*specs.accel_bandwidth;
imu_gyro_sampling_rate = 3*specs.gyro_bandwidth;

imu_rate = max(imu_accel_sampling_rate, imu_gyro_sampling_rate);

%% Standard Noise Visualization

accel_noise_std = specs.accel_noiseDensity * sqrt(specs.accel_bandwidth); % Noise standard deviation (m/s^2)
gyro_noise_std = specs.gyro_noiseDensity * sqrt(specs.gyro_bandwidth); % Noise standard deviation (rad/s)

% Getting noise values based on standard distribution
a.noiseDistAccel = accel_noise_std*randn(length(t), 1);
a.noiseDistGyro = gyro_noise_std*randn(length(t), 1);

% Noisy deck signal
a.noisyAccelCurve = @(real_accel) real_accel(t) + accel_noise_std*randn(length(t),1);
a.noisyGyroCurve = @(real_ang_rate) real_ang_rate(t) + gyro_noise_std*randn(length(t),1);

%% Adding bias instability

% But first, turn-on bias
a.accel_turn_on_bias_offset = 0.1; % m/s^2
a.gyro_turn_on_bias_offset = 0.1; % rad

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

a.biasStabDistAccel = accel_bias_norm(floor(t./dt)+1);
a.biasStabDistGyro = gyro_bias_norm(floor(t./dt)+1);

a.biasTempDistAccel = accel_bias_temp_norm(floor(t./dt)+1);
a.biasTempDistGyro = gyro_bias_temp_norm(floor(t./dt)+1);

a.theta_err = a.biasStabDistGyro.*t + ARW.*sqrt(t);

a.accel_drift_vert = @(accel, real_ang, biasStabDistAccel, biasTempDistAccel) ((1 + k)*accel + biasStabDistAccel(t)  + biasTempDistAccel(t) + accel_turn_on_bias_offset) * cos(real_ang(t));
a.accel_drift_horz = @(accel, real_ang, biasStabDistAccel, biasTempDistAccel) -((1 + k)*accel + biasStabDistAccel(t)  + biasTempDistAccel(t) + accel_turn_on_bias_offset + a.g) * cos(real_ang(t));

a.accel_drift = @(accel, biasStabDistAccel, biasTempDistAccel) (1 + k)*accel + biasStabDistAccel  + biasTempDistAccel(t) + accel_turn_on_bias_offset;
a.gyro_drift = @(ang_rate, biasStabDistGyro, biasTempDistGyro) (1 + k)*ang_rate + a.biasStabDistGyro + biasTempDistGyro(t) + gyro_turn_on_bias_offset;

%% Combining bias drift, temp bias, turn-on bias, noise, and quantization

% o_d_n_a_c_v: offset drifting noisy accel_curve vert
% o_d_n_a_c_h: offset drifting noisy accel curve horz
% o_d_n_a_c: offset drifting noisy accel curve
% o_d_n_g_c: offset drifting noisy gyro curve

a.o_d_n_a_c_v = @(biasStabDistAccel, biasTempDistAccel, accel_drift_vert, noiseDistAccel, real_accel, real_ang) accel_drift_vert(t, real_accel, real_ang, biasStabDistAccel, biasTempDistAccel) + noiseDistAccel(t);
a.o_d_n_a_c_h = @(t, biasStabDistAccel, biasTempDistAccel, accel_drift_horz, noiseDistAccel, real_accel, real_ang) accel_drift_horz(t, real_accel, real_ang, biasStabDistAccel, biasTempDistAccel) + noiseDistAccel(t);

a.o_d_n_a_c = @(biasStabDistAccel, biasTempDistAccel, accel_drift, noiseDistAccel, real_accel) accel_drift(t, real_accel, biasStabDistAccel, biasTempDistAccel) + noiseDistAccel(t);
a.o_d_n_g_c = @(biasStabDistGyro, biasTempDistGyro, gyro_drift, noiseDistGyro, real_ang_rate) gyro_drift(t, real_ang_rate, biasStabDistGyro, biasTempDistGyro) + noiseDistGyro(t);

a.measured_accel_vert = @(o_d_n_a_c_v, biasStabDistAccel, biasTempDistAccel, accel_drift_vert, noiseDistAccel, real_accel, real_ang) accel_resolution*floor(o_d_n_a_c_v(t, biasStabDistAccel, biasTempDistAccel, accel_drift_vert, noiseDistAccel, real_accel, real_ang)/accel_resolution);
a.measured_accel_horz = @(o_d_n_a_c_h, biasStabDistAccel, biasTempDistAccel, accel_drift_horz, noiseDistAccel, real_accel, real_ang) accel_resolution*floor(o_d_n_a_c_h(t, biasStabDistAccel, biasTempDistAccel, accel_drift_horz, noiseDistAccel, real_accel, real_ang)/accel_resolution);

a.measured_accel = @(o_d_n_a_c, biasStabDistAccel, biasTempDistAccel, accel_drift, noiseDistAccel, real_accel) accel_resolution*floor(o_d_n_a_c(t, biasStabDistAccel, biasTempDistAccel, accel_drift, noiseDistAccel, real_accel)/accel_resolution);
a.measured_gyro = @(o_d_n_g_c, biasStabDistGyro, biasTempDistGyro, gyro_drift, noiseDistGyro, real_ang_rate) gyro_resolution*floor(o_d_n_g_c(t, biasStabDistGyro, biasTempDistGyro, gyro_drift, noiseDistGyro, real_ang_rate)/gyro_resolution);

% measured_accel_vert = a.measured_accel_vert(t, a.o_d_n_a_c_v, a.biasStabDistAccel, a.biasTempDistAccel, a.accel_drift_vert, a.noiseDistAccel, a.real_accel, a.real_ang);
% measured_accel_horz = a.measured_accel_horz(t, a.o_d_n_a_c_h, a.biasStabDistAccel, a.biasTempDistAccel, a.accel_drift_horz, a.noiseDistAccel, a.real_accel, a.real_ang);
% measured_gyro = a.measured_gyro(t, a.o_d_n_g_c, a.biasStabDistGyro, a.biasTempDistGyro, a.gyro_drift, a.noiseDistGyro, a.real_ang_rate);

% 3D Versions of all the function handles, take vectors
a.accel_drift_3D = @(a, curr_accel) (1 + k).*curr_accel + a.biasStabDistAccel  + a.biasTempDistAccel + a.accel_turn_on_bias_offset;
a.gyro_drift_3D = @(a, curr_ang_rate) (1 + k).*curr_ang_rate + a.biasStabDistGyro + a.biasTempDistGyro + a.gyro_turn_on_bias_offset;

a.o_d_n_a_c_3D = @(a, curr_accel) a.accel_drift_3D(a, curr_accel) + a.noiseDistAccel;
a.o_d_n_g_c_3D = @(a, curr_ang_rate) a.gyro_drift_3D(a, curr_ang_rate) + a.noiseDistGyro;

a.measured_accel_3D = @(a, curr_accel) accel_resolution.*floor(a.o_d_n_a_c_3D(a, curr_accel)./accel_resolution);
a.measured_gyro_3D = @(a, curr_ang_rate) gyro_resolution.*floor(a.o_d_n_g_c_3D(a, curr_ang_rate)./gyro_resolution);
