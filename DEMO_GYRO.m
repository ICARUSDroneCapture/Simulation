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

fprintf('\nIntegrating Angular Rate (no Sensor Error)')
fprintf("\nTime: ")

% Initial Angle
p_theta0 = a.real_ang(tspan(1)); % Platform inertial angle [deg]

[t, theta_base]= rk4_solver(a.real_ang_rate, tspan, p_theta0, dt);

% Integration Verification Plot
figure
plot(t, 180/pi*theta_base)
hold on
plot(t, 180/pi*a.real_ang(t))
xlabel('Time (sec)')
ylabel('Angle (deg)')
legend('Integrated Angle', 'Angle Equation')
title('Basic Integrated Angle Plot')

fprintf('\nFinished Integration with NO Sensor Error.\n')
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

control_dynamics = @(t_i, state) a.measured_gyro(t_i, a.o_d_n_g_c, a.biasStabDistGyro, a.biasTempDistGyro, a.gyro_drift, a.noiseDistGyro, a.real_ang_rate(t_i));

[t_error, theta_error]= rk4_solver(control_dynamics, tspan, p_theta0, dt);

fprintf('\nFinished Integration WITH Sensor Error.\n')

%% Plotting Platform Angle

figure
plot(t, 180/pi*theta_base)
hold on
plot(t_error, 180/pi*theta_error)
xlabel('Time (sec)')
ylabel('Angle (deg)')
legend('NO Sensor Error', 'WITH Sensor Error')
title('Basic Integrated Angle Plot')
ylim([-30 30])

%% Getting and Plotting Angle Error

pos_err = theta_error(:, 1) - theta_base(:, 1);

sz = 2;
plot_scale = 10;

figure
scatter(t_error, pos_err*180/pi, sz, 'filled', displayName="Angle Error")
title('Integrated Angle Error vs Time')
xlabel('Time (s)')
ylabel('Error (deg)')

%% Control Law Drift Compensation

fprintf('\nStarting Integration WITH Sensor Compensation Control Law.')
fprintf("\nTime: ")

control_dynamics = @(t_i, state) GyroDriftCorrection(t_i, a, state);

theta_0 = 0;
theta_dot_0 = 0;

s0 = [theta_0 theta_dot_0];

[t_corr, theta_corr]= rk4_solver(control_dynamics, tspan, s0, dt);

fprintf('\nFinished Integration WITH Sensor Compensation Control Law.\n')

theta_der = diff(theta_corr) / dt;
figure
plot(t(2:end), theta_der(:, 1))
hold on
plot(t, a.real_ang_rate(t))
xlabel('Time (sec)')
ylabel('Angular Rate (rad/s)')
legend('Controlled Angular Rate', 'Real Deck Angular Rate')

figure
plot(t, 180/pi*theta_base)
hold on
plot(t_error, 180/pi*theta_error)
hold on
plot(t_corr, 180/pi*theta_corr(:, 1))
xlabel('Time (sec)')
ylabel('Angle (deg)')
legend('NO Sensor Error', 'WITH Sensor Error', 'Corrected Sensor Error')
title('Basic Integrated Angle Plot')
ylim([-30 30])


function s_dot = GyroDriftCorrection(time_i, a, prev_state)

    % Current states
    theta = prev_state(1);
    theta_err_accum = prev_state(2);

    specs = a.specs;

    kw = a.kw;
    kt = a.kt;

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
