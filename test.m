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

real_vel = @(t, y) 1/3*t.^3;
real_ang = @(t, y) 1/6*t.^4; % [deg]

real_accel = @(t, y) t.^2; % [m*s^-2]
real_ang_rate = @(t, y) 2/3*t.^3; % [deg/s]

dynamics = @(t, y) [ t.^2; 2/3*t.^3 ]; % [ real_accel real_ang_rate]

%% Sensor Model Aspects

% imx_5_specs
imx_5_specs_deg

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
finishTime = 20;
tspan = [startTime finishTime]; % [s]

%solving the system
% freq = 160; % Hz
% dt = 1 / freq; % Timestep (s)

imu_accel_sampling_rate = 3*specs.accel_bandwidth;
imu_gyro_sampling_rate = 3*specs.gyro_bandwidth;

imu_rate = max(imu_accel_sampling_rate, imu_gyro_sampling_rate);

% dt = 1/imu_rate;  % [s]
dt = 0.0001;
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

a.theta_err = @(t) a.biasStabDistGyro(t).*t + ARW.*sqrt(t);

a.accel_drift_vert = @(t, accel, theta_err, biasStabDistAccel, biasTempDistAccel) (1 + k)*accel + biasStabDistAccel(t)  + biasTempDistAccel(t) - a.g + accel_turn_on_bias_offset;
a.accel_drift_horz = @(t, accel, theta_err, biasStabDistAccel, biasTempDistAccel) (1 + k)*accel + biasStabDistAccel(t)  + biasTempDistAccel(t) + g*sin(theta_err(t)) + accel_turn_on_bias_offset;

a.gyro_drift = @(t, ang_rate, biasStabDistGyro, biasTempDistGyro) (1 + k)*ang_rate + a.biasStabDistGyro(t) + biasTempDistGyro(t) + gyro_turn_on_bias_offset;

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


%% Run Control Dynamics Integration

fprintf('\nStarting Integration with NO Sensor Error.')
fprintf("\nTime: ")

% % Initial States
% p0 =  a.d(tspan(1))+a.pr_d;   % Platform position [m]
% p_dot0 = 0;   % Platform velocity [m/s]
% pr_err_accum0 = 0;            % Integral of relative position error [m*s]
% pm0 = p0;                     % Platform inetegrated position [m]
% pm_dot = p_dot0;              % Platform integrated velocity [m/s]
% pm_ddot = a.d_ddot(tspan(1)); % Platform measured acceleration [m*s^-2]
% p_theta0 = a.real_ang(tspan(1)); % Platform inertial angle [deg]

% Initial States
p0 =  a.d(tspan(1))+a.pr_d;   % Platform position [m]
p_dot0 = 0;   % Platform velocity [m/s]
pr_err_accum0 = 0;            % Integral of relative position error [m*s]
pm0 = p0;                     % Platform inetegrated position [m]
pm_dot = p_dot0;              % Platform integrated velocity [m/s]
pm_ddot = 0; % Platform measured acceleration [m*s^-2]
p_theta0 = 0; % Platform inertial angle [deg]
p_theta_err_accum0 = 0;

s0 = [p0 p_dot0 pr_err_accum0 pm0 pm_dot pm_ddot];

control_dynamics = @(t, state) NoError_FixedInt(t, a, state);

[t_control, sol_control]= rk4_solver(control_dynamics, tspan, s0, dt);


fprintf('\nFinished Integration with NO Sensor Error.\n')

%% Running Control Law Simulation WITH Sensor Error

fprintf('\nStarting Integration WITH Sensor Error.')
fprintf("\nTime: ")

% Running simulation with sensor error

s0 = [p0 p_dot0 pr_err_accum0 pm0 pm_dot pm_ddot p_theta0 p_theta_err_accum0];
control_dynamics_err = @(t, state) rigidArmControl_FixedInt(t, a, state);

% % Initial States
% p0 =  a.pr_d;   % Platform position [m]
% p_dot0 = 0;   % Platform velocity [m/s]
% pr_err_accum0 = 0;            % Integral of relative position error [m*s]
% pm0 = p0;                     % Platform inetegrated position [m]
% pm_dot = p_dot0;              % Platform integrated velocity [m/s]
% pm_ddot = 0; % Platform measured acceleration [m*s^-2]
% p_theta0 = 0; % Platform inertial angle [deg]
% p_theta_err_accum0 = 0;
% 
% % Relative Position Control
% a.kp = 8;  % Proportional [N/m]
% a.kd = 1;  % Derivative [Ns/m]    
% a.ki = 0;  % Integral [N/ms]
% 
% % % Inertial Stabilization Control
% a.ka =  10;  % Acceleration Control [kg]
% a.kv = 0;  % Velocity Control [kg/s]
% a.ks = 0;  % Position Control [kg*s^-2]
% 
% % Sensor Drift Control
% a.kw = 0.9; % 
% a.kt = 0.9; % 
% 
% s0 = [p_dot0 p_theta0 p0 p_theta_err_accum0 pr_err_accum0];
% control_dynamics_err = @(t, state) ControlSignal(t, a, state);

% Outputs in sensor (inertial) frame, where zero is defined AT platform
[t_error, sol_error]= rk4_solver(control_dynamics_err, tspan, s0, dt);

% Get platform position in inertial frame, with deck as reference zero
% plat_pos = sol_error(:, 3) + a.d(t);
plat_pos = sol_error(:,1);

fprintf('\nFinished Integration WITH Sensor Error.\n')

% %% Plotting other states
% 
% figure
% plot(t_error, sol_error(:, 1))
% xlabel('Time (sec)')
% ylabel('Velocity (m/s)')
% title('Platform Inertial Velocity')
% 
% figure
% plot(t_error, sol_error(:, 2))
% xlabel('Time (sec)')
% ylabel('Angle (rad)')
% title('Platform Angle')

%% Plotting Platform Position

figure;
plot(t_control, sol_control(:,1))
hold on
plot(t_error, plat_pos)
hold on
plot(t,a.d(t))
hold on
plot(t, a.d(t)+1)
hold on
plot(t, a.d(t)+0.09, '--')
hold on
plot(t, a.d(t)+0.5+0.41, '--')
title('Platform Inertial Position vs Time')
xlabel('Time (s)')
ylabel('Position (m)')
title('Platform Inertial Position over Time')
legend('Fixed-Step (without sensor error) Integration', 'Fixed-Step (with sensor error) Integration')


%% Getting and Plotting error

pos_err = plat_pos - sol_control(:, 1);

sz = 2;
plot_scale = 0.001;

figure
scatter(t_error, pos_err*100, sz, 'filled', displayName="Positional Error")
% hold on
% plot(t,a.d(t)/200, displayName="Deck Disturbance")
% hold on
% plot(t,plot_scale*a.d(t))
% hold on
% plot(t, plot_scale*(a.d(t)+1))
% hold on
% plot(t, plot_scale*(a.d(t)+0.09), '--')
% hold on
% plot(t, plot_scale*(a.d(t)+0.5+0.41), '--')
title('Worst Case Relative Position Error vs Time')
xlabel('Time (s)')
ylabel('Error (cm)')
legend

plot_scale = 0.00001;
growth = diff(pos_err);

% figure
% scatter(t_error(2:end), growth*100, sz, 'filled', displayName="Positional Error")
% hold on
% plot(t,plot_scale*a.d(t))
% hold on
% plot(t, plot_scale*(a.d(t)+1))
% hold on
% plot(t, plot_scale*(a.d(t)+0.09), '--')
% hold on
% plot(t, plot_scale*(a.d(t)+0.5+0.41), '--')
% % hold on
% % plot(t,a.d(t)/200, displayName="Deck Disturbance")
% title('Error Growth over Time')
% xlabel('Time (s)')
% ylabel('Error (cm)')
% legend
% 
% 
