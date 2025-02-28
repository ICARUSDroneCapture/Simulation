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

%% To Demo Stable Control

% alpha = 0.45; % wave amplitdue [m]
% 
% % Relative Position Control
% a.kp = 10;  % Proportional [N/m]
% a.kd = 2;  % Derivative [Ns/m]    
% a.ki = 0;  % Integral [N/ms]
% 
% % Inertial Stabilization Control
% a.ka =  0.98;  % Acceleration Control [kg]
% a.kv = 2.71;  % Velocity Control [kg/s]
% a.ks = 0;  % Position Control [kg*s^-2]
% 
% finishTime = 20; % Time can be changed to anything

% Alternative Gains, to show integrator strength

% % Relative Position Control
% a.kp = 80;  % Proportional [N/m]
% a.kd = 10;  % Derivative [Ns/m]    
% a.ki = 0;  % Integral [N/ms]
% 
% % Inertial Stabilization Control
% a.ka =  1.2;  % Acceleration Control [kg]
% a.kv = 2.8;  % Velocity Control [kg/s]
% a.ks = 0;  % Position Control [kg*s^-2]

% Relative Position Control
a.kp = 8;  % Proportional [N/m]
a.kd = 1;  % Derivative [Ns/m]    
a.ki = 0;  % Integral [N/ms]

% Inertial Stabilization Control
a.ka =  1.5;  % Acceleration Control [kg]
a.kv = 2.5;  % Velocity Control [kg/s]
a.ks = 0;  % Position Control [kg*s^-2]

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
finishTime = 20; %
tspan = [startTime finishTime]; % [s]

% 1 & 3 work at 0.0001 seconds
% 2 & 4 & 5 work at 0.00001 seconds
dt_base = 0.01;  % [s]
dt_big = 0.0001;  % [s]
dt_small = 0.00005;  % [s]
control_integrator_type = 3;

t = (tspan(1):dt_base:tspan(2))';

t_count = length(t);

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

max_vel = max(a.real_vel(t));

%% Standard Noise Visualization

accel_noise_std = specs.accel_noiseDensity * sqrt(specs.accel_bandwidth); % Noise standard deviation (m/s^2)
gyro_noise_std = specs.gyro_noiseDensity * sqrt(specs.gyro_bandwidth); % Noise standard deviation (dps)

% Plotting noise distributions

a.noiseDistAccel = @(t) accel_noise_std*randn(length(t),1);
a.noiseDistGyro = @(t) gyro_noise_std*randn(length(t),1);

% Noisy signal
a.noisyAccelCurve = @(t, real_accel) real_accel(t) + accel_noise_std*randn(length(t),1);
a.noisyGyroCurve = @(t, real_ang_rate) real_ang_rate(t) + gyro_noise_std*randn(length(t),1);

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

% bias_indeces = floor(t./(length(t)/driftAlterations)*100)+1;

a.biasStabDistAccel = @(t) b_a_drift_vals(floor(t./(t_count/driftAlterations)*100)+1);
a.biasStabDistGyro = @(t) b_g_drift_vals(floor(t./(t_count/driftAlterations)*100)+1);

a.biasedAccelCurve = @(t, real_accel) real_accel(t) + b_a_drift_vals(floor(t./(length(t)/driftAlterations)*100)+1);
a.biasedGyroCurve = @(t, real_ang_rate) real_ang_rate(t) + b_g_drift_vals(floor(t./(length(t)/driftAlterations)*100)+1);

% Bias Error over Temp

% a.biasTempDistAccel = @(t) accel_temp_bias*randn(length(t),1);
% a.biasTempDistGyro = @(t) gyro_temp_bias*randn(length(t),1);

b_a_temp_vals = b_a.*rand(driftAlterations,1);
b_g_temp_vals = b_g.*rand(driftAlterations,1);

a.biasTempDistAccel = @(t) b_a_temp_vals(floor(t./(t_count/driftAlterations)*100)+1);
a.biasTempDistGyro = @(t) b_g_temp_vals(floor(t./(t_count/driftAlterations)*100)+1);

a.theta_err = @(t) a.biasStabDistGyro(t).*t + ARW.*sqrt(t);

a.accel_drift_vert = @(t, real_accel, theta_err, biasStabDistAccel, biasTempDistAccel) (1 + k)*real_accel(t) + a.biasStabDistAccel(t)  + biasTempDistAccel(t) + g*(1-cos(theta_err(t)));
a.accel_drift_horz = @(t, real_accel, theta_err, biasStabDistAccel, biasTempDistAccel) (1 + k)*real_accel(t) + a.biasStabDistAccel(t)  + biasTempDistAccel(t)  + g*sin(theta_err(t));

a.gyro_drift = @(t, real_ang_rate, biasStabDistGyro, biasTempDistGyro) (1 + k)*real_ang_rate(t) + a.biasStabDistGyro(t) + biasTempDistGyro(t);

%% Combining bias drift, temp bias, turn-on bias, noise, and quantization

accel_turn_on_bias_offset = normrnd(0, accel_noise_std);
gyro_turn_on_bias_offset = normrnd(0, gyro_noise_std);

% o_d_n_a_c_v: offset drifting noisy accel_curve vert
% o_d_n_a_c_h: offset drifting noisy accel curve horz
% o_d_n_g_c: offset drifting noisy gyro curve

a.o_d_n_a_c_v = @(t, biasStabDistAccel, biasTempDistAccel, accel_drift_vert, noiseDistAccel, real_accel, theta_err) accel_drift_vert(t, real_accel, theta_err, biasStabDistAccel, biasTempDistAccel) + noiseDistAccel(t) + accel_turn_on_bias_offset;
a.o_d_n_a_c_h = @(t, biasStabDistAccel, biasTempDistAccel, accel_drift_horz, noiseDistAccel, real_accel, theta_err) accel_drift_horz(t, real_accel, theta_err, biasStabDistAccel, biasTempDistAccel) + noiseDistAccel(t) + accel_turn_on_bias_offset;
a.o_d_n_g_c = @(t, biasStabDistGyro, biasTempDistGyro, gyro_drift, noiseDistGyro, real_ang_rate) gyro_drift(t, real_ang_rate, biasStabDistGyro, biasTempDistGyro) + noiseDistGyro(t) + gyro_turn_on_bias_offset;

a.measured_accel_vert = @(t, o_d_n_a_c_v, biasStabDistAccel, biasTempDistAccel, accel_drift_vert, noiseDistAccel, real_accel, theta_err) accel_resolution*floor(o_d_n_a_c_v(t, biasStabDistAccel, biasTempDistAccel, accel_drift_vert, noiseDistAccel, real_accel, theta_err)/accel_resolution);
a.measured_accel_horz = @(t, o_d_n_a_c_h, biasStabDistAccel, biasTempDistAccel, accel_drift_horz, noiseDistAccel, real_accel, theta_err) accel_resolution*floor(o_d_n_a_c_h(t, biasStabDistAccel, biasTempDistAccel, accel_drift_horz, noiseDistAccel, real_accel, theta_err)/accel_resolution);
a.measured_gyro = @(t, o_d_n_g_c, biasStabDistGyro, biasTempDistGyro, gyro_drift, noiseDistGyro, real_ang_rate) gyro_resolution*floor(o_d_n_g_c(t, biasStabDistGyro, biasTempDistGyro, gyro_drift, noiseDistGyro, real_ang_rate)/gyro_resolution);

measured_accel_vert = a.measured_accel_vert(t, a.o_d_n_a_c_v, a.biasStabDistAccel, a.biasTempDistAccel, a.accel_drift_vert, a.noiseDistAccel, a.real_accel, a.theta_err);
measured_accel_horz = a.measured_accel_horz(t, a.o_d_n_a_c_h, a.biasStabDistAccel, a.biasTempDistAccel, a.accel_drift_horz, a.noiseDistAccel, a.real_accel, a.theta_err);
measured_gyro = a.measured_gyro(t, a.o_d_n_g_c, a.biasStabDistGyro, a.biasTempDistGyro, a.gyro_drift, a.noiseDistGyro, a.real_ang_rate);

%% Running Fixed-Integration Simulations

% Initial States
p0 =  a.d(tspan(1))+a.pr_d;   % Platform position [m]
p_dot0 = max_vel;   % Platform velocity [m/s]
pr_err_accum0 = 0;            % Integral of relative position error [m*s]
pm0 = p0;                     % Platform inetegrated position [m]
pm_dot = p_dot0;              % Platform integrated velocity [m/s]
pm_ddot = a.d_ddot(tspan(1)); % Platform measured acceleration [m*s^-2]
p_theta0 = a.real_ang(tspan(1)); % Platform inertial angle [deg]

% --------------------- Running no error simulation -----------------------

fprintf('\nBeginning Fixed-Step Integration (without sensor error).\n\n');

t_fix_control = (tspan(1):dt_big:tspan(2))';

prev_state = [p0; p_dot0; pr_err_accum0; pm0; pm_dot; pm_ddot];

% Columns represent each variable, so [position, velocity, relative position err, ...]
% Rows represents timesteps, max limit of 5 timesteps kept track of
state_record_control = zeros(5, 6);
deriv_record_control = zeros(5, 6);
s_fix_control = zeros(length(t_fix_control), 6);

state_record_control = insertVector(state_record_control, prev_state');

fprintf("Time: ")

for i = 1:length(t_fix_control)

    if i == 1
        prev_state(2) = 0;
    end

    time = t_fix_control(i);

    val_dot_control = NoError_FixedInt(time, a, prev_state)';

    deriv_record_control = insertVector(deriv_record_control, val_dot_control);

    state_vec_control = fdm_integrator(state_record_control, deriv_record_control, dt_big, control_integrator_type);

    state_record_control = insertVector(state_record_control, state_vec_control);

    s_fix_control(i, :) = state_vec_control;
    prev_state = state_vec_control;

    if mod(time, 1) == 0
        fprintf("%i  ", time)
    end
end

fprintf('\n\nFixed Integration (without sensor error) finished.\nBeginning fixed integration (with sensor error).\n\n');

figure;
plot(t_fix_control, s_fix_control(:,1))
title('Platform Inertial Position vs Time')
xlabel('Time (s)')
ylabel('Position (m)')
title('Platform Inertial Position over Time')
legend('Fixed-Step (without sensor error) Integration')

% --------------------- Running with error simulation -----------------------

t_fix_error = (tspan(1):dt_small:tspan(2))';

prev_state = [p0; p_dot0; pr_err_accum0; pm0; pm_dot; pm_ddot; p_theta0];

% Columns represent each variable, so [position, velocity, relative position err, ...]
% Rows represents timesteps, max limit of 5 timesteps kept track of
state_record_error = zeros(5, 7);
deriv_record_error = zeros(5, 7);
s_fix_error = zeros(length(t_fix_error), 7);

state_record_error = insertVector(state_record_error, prev_state');

fprintf("Time: ")

for i = 1:length(t_fix_error)

    if i == 1
        prev_state(2) = 0;
    end

    time = t_fix_error(i);

    val_dot_error = rigidArmControl_FixedInt(time, a, prev_state)';

    deriv_record_error = insertVector(deriv_record_error, val_dot_error);

    state_vec_error = fdm_integrator(state_record_error, deriv_record_error, dt_small, control_integrator_type);

    state_record_error = insertVector(state_record_error, state_vec_error);

    s_fix_error(i, :) = state_vec_error;
    prev_state = state_vec_error;

    if mod(time, 1) == 0
        fprintf("%i  ", time)
    end
end

fprintf('\n\nFixed Integration (with sensor error) finished.\n\n');

figure;
plot(t_fix_control, s_fix_control(:,1))
hold on
plot(t_fix_error, s_fix_error(:,1))
title('Platform Inertial Position vs Time')
xlabel('Time (s)')
ylabel('Position (m)')
title('Platform Inertial Position over Time')
legend('Fixed-Step (without sensor error) Integration', 'Fixed-Step (with sensor error) Integration')

%% Plotting

% With error ------------------------------------------------------------

% Plotting Position vs Time and Acceleration vs Time
figure;

% Position
plot(t_fix_control, s_fix_control(:,1))
hold on
plot(t_fix_error, s_fix_error(:,1))
hold on
plot(t,a.d(t))
hold on
plot(t, a.d(t)+1)
hold on
plot(t, a.d(t)+0.09, '--')
hold on
plot(t, a.d(t)+0.5+0.41, '--')
title('Inertial Position vs Time')
xlabel('Time (s)')
ylabel('Position (m)')
legend('With Sensor Error', 'Without Sensor Error', 'Deck','Location','southeast')

% Plotting relative position
figure;
plot(t_fix_error,s_fix_error(:,1)-a.d(t_fix_error))
hold on
% Plotting inertial control region
x = [tspan, flip(tspan)];
yf = [a.pr_d-a.r_g, a.pr_d-a.r_g, a.pr_d+a.r_g, a.pr_d+a.r_g];
fill(x,yf,'y','FaceAlpha',0.2,'EdgeColor','none')
% Plotting Relative position control region
x = [tspan, flip(tspan)];
yta = [a.pr_d+a.r_k, a.pr_d+a.r_k, 1, 1];
ytb = [0, 0, a.pr_d-a.r_k, a.pr_d-a.r_k];
fill(x,yta,'b','FaceAlpha',0.2,'EdgeColor','none')
fill(x,ytb,'b','FaceAlpha',0.2,'EdgeColor','none')
yline(a.pr_d,'--','Label','$p_{rd}$','Interpreter','latex','FontSize',15)
title('Relative Position vs Time')
xlabel('Time (s)')
ylabel('Position (m)')
legend('','Full Inertial','Relative Position','')



%% Getting error

% Plotting time
sz = 2;
figure
scatter(1:length(t_fix_error), t_fix_error, sz, 'filled', displayName="With Error")
hold on
scatter(1:length(t_fix_control), t_fix_control, sz, 'filled', displayName="Without Error")
ylabel('Time Values')
title('Timesteps used in Integration')
legend

% Precision (number of decimals) of interpolation
err_round = 3;

% Getting the more precise time matrix
if size(t_fix_control) > size(t_fix_error)
    t_precise = t_fix_control;
    p_precise = s_fix_control(:,1);
    t_compare = t_fix_error;
    p_compare = s_fix_error(:,1);
else
    t_precise = t_fix_error;
    p_precise = s_fix_error(:,1);
    t_compare = t_fix_control;
    p_compare = s_fix_control(:,1);
end

finishIndex = length(t_precise);
pos_err = zeros(1, finishIndex);

for i=1:finishIndex

    time_precise = t_precise(i);
    pos_precise = p_precise(i);

    t_ref = round(time_precise, err_round);
    time_interp_idx = findNearest(t_ref, t_compare);

    time_interp = t_compare(time_interp_idx);
    pos_interp = p_compare(time_interp_idx);

    pos_err(i) = abs(pos_precise - pos_interp);
end

%% Plotting Error

sz = 2;

figure
scatter(t_precise, pos_err*100, sz, 'filled', displayName="Positional Error")
hold on
plot(t,a.d(t)/200, displayName="Deck Disturbance")
title('Worst Case Relative Position Error vs Time')
xlabel('Time (s)')
ylabel('Error (cm)')
legend

growth = diff(pos_err);

figure
scatter(t_precise(2:end), growth*100, sz, 'filled', displayName="Positional Error")
hold on
plot(t,a.d(t)/200, displayName="Deck Disturbance")
title('Error Growth over Time')
xlabel('Time (s)')
ylabel('Error (cm)')
legend