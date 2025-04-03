close all; clear; clc;

rng(1,"twister");

% Use same random seed

set(groot,'DefaultLineLineWidth',1)

simulationParameters;

close all;

% Redefining acceleration/gyro curves for clarity

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
% defineSignalsNoNoise

scale_w = 1;
scale_t = 1;

% a.kw = scale_w*a.beta_min; % 
% a.kt = scale_t*a.beta_max; % 

T_x = 2*pi/(beta/2); % wave frequency [rad/s] (15 sec)
T_y = 2*pi/(beta/4); % wave frequency [rad/s] (30 sec)
T_z = 2*pi/(beta); % wave frequency [rad/s] (7.5 sec)

t_min_x = 0.1;
t_max_x = 18;

t_min_y = 0.1;
t_max_y = 30.1;

t_min_z = 7.1;
t_max_z = 30;

a.beta_min_x = 1/t_min_x;
a.beta_min_y = 1/t_min_y;
a.beta_min_z = 1/t_min_z;

a.beta_max_x = 1/t_max_x;
a.beta_max_y = 1/t_max_y;
a.beta_max_z = 1/t_max_z;

a.kw = [scale_t*a.beta_max_x; scale_t*a.beta_max_y; scale_t*a.beta_max_z; scale_t*a.beta_max_x; scale_t*a.beta_max_y; scale_t*a.beta_max_z];
a.kt = [scale_w*a.beta_min_x; scale_w*a.beta_min_y; scale_w*a.beta_min_z; scale_w*a.beta_min_x; scale_w*a.beta_min_y; scale_w*a.beta_min_z];

%% 3D motion equations

a.real_pos_xI = @(t) 0.4/(beta^2)*sin(beta*t/2);
a.real_pos_yI = @(t) 1.6/(beta^2)*sin(beta*t/4);
a.real_pos_zI = @(t) alpha*sin(beta*t) + hdeck;

a.real_vel_xI = @(t) 0.2/beta*cos(beta/2*t);
a.real_vel_yI = @(t) 0.4/beta*cos(beta/4*t);
a.real_vel_zI = @(t) beta*alpha*cos(beta*t);

a.real_accel_xI = @(t) -0.1*sin(beta/2*t); % [m*s^-2]
a.real_accel_yI = @(t) -0.1*sin(beta/4*t); % [m*s^-2]
a.real_accel_zI = @(t) -beta^2*alpha*sin(beta*t) - 9.81; % [m*s^-2]

a.theta = @(t) -atan(beta*alpha*cos(beta*t)); % [rad]
a.phi = @(t) atan(0.2/beta*cos(beta/2*t)); % [rad]
a.psi = @(t) atan(0.4/beta*cos(beta/4*t)); % [rad]

% a.theta = @(t) atan(0.4/beta*cos(beta/4*t)); % [rad]
% a.phi = @(t) -atan(beta*alpha*cos(beta*t)); % [rad]
% a.psi = @(t) atan(0.2/beta*cos(beta/2*t)); % [rad]

% a.theta_dot = @(t) ((alpha*beta^2*sin(beta*t))./(alpha^2*beta^2*(cos(beta*t).^2)+1)); % [rad/s]
% a.phi_dot = @(t) ((-0.1*sin(beta/2*t))/(((0.04*(cos(beta*t/2).^2))/(beta^2))+1)); % [rad/s]
% a.psi_dot = @(t) ((-0.1*sin(beta/4*t))/(((0.16*(cos(beta*t/4).^2))/(beta^2))+1)); % [rad/s]

% Angular rate needs to be taken manually since the derivative equation is +/-
theta_vals = a.theta(t);
phi_vals = a.phi(t);
psi_vals = a.psi(t);

theta_dot = zeros(1, length(t));
phi_dot = zeros(1, length(t));
psi_dot = zeros(1, length(t));

theta_dot(2:end) = diff(theta_vals)/dt;
phi_dot(2:end) = diff(phi_vals)/dt;
psi_dot(2:end) = diff(psi_vals)/dt;

a.theta_dot_eq = @(t) theta_dot(floor(t./dt)+1);
a.phi_dot_eq = @(t) phi_dot(floor(t./dt)+1);
a.psi_dot_eq = @(t) psi_dot(floor(t./dt)+1);

s = @(x) sin(x);
c = @(x) cos(x);

%% Simplifying to just z-direction equations for now
a.real_pos = @(t) alpha*sin(beta*t) + hdeck;
a.real_vel = @(t) beta*alpha*cos(beta*t);
a.real_accel = @(t, y) -beta^2*alpha*sin(beta*t) - 9.81; % [m*s^-2]
a.real_ang = @(t) atan(0.4/beta*cos(beta/4*t)); % [rad]
a.real_ang_rate = @(t, y) theta_dot(floor(t./dt)+1);

%% Sensor Frame Accelerations

a_S_true = zeros(3, length(t));

for i = 1:length(t)

    t_i = t(i);

    theta = a.theta(t_i);
    phi = a.phi(t_i);
    psi = a.psi(t_i);

    a_I = [a.real_accel_xI(t_i); a.real_accel_yI(t_i); a.real_accel_zI(t_i)];

    a_S = Rotate_I_S(a_I, theta, phi, psi);

    a_S_true(:, i) = a_S;

end

%% Get Sensor Frame Measurements

% Sensor frame values determined using REAL angle, with sensor error added

a_S_measured = zeros(3, length(t));

for i = 1:length(t)

    t_i = t(i);

    theta = a.theta(t_i);
    phi = a.phi(t_i);
    psi = a.psi(t_i);

    a_I = [a.real_accel_xI(t_i); a.real_accel_yI(t_i); a.real_accel_zI(t_i)];

    a_S = Rotate_I_S(a_I, theta, phi, psi);

    a_Sx_m = a.measured_accel(t_i, a.o_d_n_a_c, a.biasStabDistAccel, a.biasTempDistAccel, a.accel_drift, a.noiseDistAccel, a_S(1));
    a_Sy_m = a.measured_accel(t_i, a.o_d_n_a_c, a.biasStabDistAccel, a.biasTempDistAccel, a.accel_drift, a.noiseDistAccel, a_S(2));
    a_Sz_m = a.measured_accel(t_i, a.o_d_n_a_c, a.biasStabDistAccel, a.biasTempDistAccel, a.accel_drift, a.noiseDistAccel, a_S(3));
   
    a_S_measured(:, i) = [a_Sx_m; a_Sy_m; a_Sz_m];

end

figure
subplot(3,1,1)
plot(t, a.real_accel_xI(t))
hold on
plot(t, a_S_measured(1, :))
xlabel('Time (sec)')
ylabel('Acceleration (m/s^2)')
title('Measured Sensor Acceleration X (m/s^2)')
legend('Inertial Frame Acceleration', 'Sensor Frame Acceleration')

subplot(3,1,2)
plot(t, a.real_accel_yI(t))
hold on
plot(t, a_S_measured(2, :))
xlabel('Time (sec)')
ylabel('Acceleration (m/s^2)')
title('Measured Sensor Acceleration Y (m/s^2)')
legend('Inertial Frame Acceleration', 'Sensor Frame Acceleration')

subplot(3,1,3)
plot(t, a.real_accel_zI(t))
hold on
plot(t, a_S_measured(3, :))
xlabel('Time (sec)')
ylabel('Acceleration (m/s^2)')
title('Measured Sensor Acceleration Z (m/s^2)')
legend('Inertial Frame Acceleration', 'Sensor Frame Acceleration')

%% Convert Error Signal Back to Inertial Without Offset Correction

a_I_ref = zeros(3, length(t));

for i = 1:length(t)

    t_i = t(i);

    a_S = a_S_measured(:, i);

    theta = a.theta(t_i);
    phi = a.phi(t_i);
    psi = a.psi(t_i);

    a_I = Rotate_S_I(a_S, theta, phi, psi);

    a_I_ref(:, i) = a_I;

end

figure
subplot(3,1,1)
plot(t, a.real_accel_xI(t))
hold on
plot(t, a_I_ref(1, :))
xlabel('Time (sec)')
ylabel('Acceleration (m/s^2)')
title('Inertial Acceleration X Converted From Sensor (m/s^2)')
legend('Inertial Frame Acceleration', 'Sensor Frame Acceleration')

subplot(3,1,2)
plot(t, a.real_accel_yI(t))
hold on
plot(t, a_I_ref(2, :))
xlabel('Time (sec)')
ylabel('Acceleration (m/s^2)')
title('Inertial Acceleration Y Converted From Sensor (m/s^2)')
legend('Inertial Frame Acceleration', 'Sensor Frame Acceleration')

subplot(3,1,3)
plot(t, a.real_accel_zI(t))
hold on
plot(t, a_I_ref(3, :))
xlabel('Time (sec)')
ylabel('Acceleration (m/s^2)')
title('Inertial Acceleration Z Converted From Sensor (m/s^2)')
legend('Inertial Frame Acceleration', 'Sensor Frame Acceleration')

% %% Get Error Signal
% 
% fprintf('\nStarting Integration WITH Sensor Error.')
% fprintf("\nTime: ")
% 
% % Running simulation with sensor error
% 
% error_signals = zeros(t_count, 6);
% 
% for i = 1:length(t)
%     time_i = t(i);
% 
%     ang_rate_i = a.real_ang_rate(time_i);
%     accel_i = a.real_accel(time_i);
% 
%     measured_gyro = a.measured_gyro(time_i, a.o_d_n_g_c, a.biasStabDistGyro, a.biasTempDistGyro, a.gyro_drift, a.noiseDistGyro, ang_rate_i);
%     measured_accel_v = a.measured_accel_vert(time_i, a.o_d_n_a_c_v, a.biasStabDistAccel, a.biasTempDistAccel, a.accel_drift_vert, a.noiseDistAccel, accel_i, a.theta_err);
%     measured_accel_h = a.measured_accel_horz(time_i, a.o_d_n_a_c_h, a.biasStabDistAccel, a.biasTempDistAccel, a.accel_drift_horz, a.noiseDistAccel, accel_i, a.theta_err);
% 
%     measuredState = [measured_accel_h measured_accel_h measured_accel_v measured_gyro measured_gyro measured_gyro];
% 
%     error_signals(i, :) = measuredState;
% end
% 
% theta0_x = 0;
% theta0_y = 0;
% theta0_z = 0;
% vel0_x = 0;
% vel0_y = 0;
% vel0_z = 0;
% 
% s0 = [theta0_x theta0_y theta0_z vel0_x vel0_y vel0_z];
% 
% control_dynamics = @(t_i, state) error_signals(floor(t_i./dt)+1)';
% 
% [t_error, int_signals_error]= rk4_solver(control_dynamics, tspan, s0, dt);
% 
% fprintf('\nFinished Integration WITH Sensor Error.\n')
% 
% figure
% plot(t_error, error_signals(:, 1))
% hold on
% plot(t_error, error_signals(:, 3))
% hold on
% plot(t, a.real_accel(t))
% xlabel('Time (sec)')
% ylabel('Acceleration (m/s^2)')
% legend('WITH Sensor Error (horizontal)', 'WITH Sensor Error (vertical)', 'NO Sensor Error')
% title('Measured Accelerometer Signal')
% % ylim([-5 5])
% 
% figure
% plot(t_error, error_signals(:, 4))
% hold on
% plot(t, a.real_ang_rate(t))
% xlabel('Time (sec)')
% ylabel('Angular Rate (deg/s)')
% legend('WITH Sensor Error', 'NO Sensor Error')
% title('Measured Gyroscope Signal')
% % ylim([-5 5])
% 
% %% Error Compensation
% 
% fprintf('\nGetting Corrected Signal')
% 
% corrected_signals = zeros(length(t_count), 6);
% 
% for i = 1:length(t)
%     time_i = t(i);
% 
%     accel_i = a.real_accel(time_i);
%     ang_rate_i = a.real_ang_rate(time_i);
%     curr_angle = a.real_ang(time_i);
% 
%     measured_gyro = a.measured_gyro(time_i, a.o_d_n_g_c, a.biasStabDistGyro, a.biasTempDistGyro, a.gyro_drift, a.noiseDistGyro, ang_rate_i);
% 
%     measured_accel_v = a.measured_accel_vert(time_i, a.o_d_n_a_c_v, a.biasStabDistAccel, a.biasTempDistAccel, a.accel_drift_vert, a.noiseDistAccel, accel_i, a.theta_err);
%     measured_accel_h = a.measured_accel_horz(time_i, a.o_d_n_a_c_h, a.biasStabDistAccel, a.biasTempDistAccel, a.accel_drift_horz, a.noiseDistAccel, accel_i, a.theta_err);
%     measured_accel_h = measured_accel_h / cos(curr_angle);
%     measured_accel_v = -measured_accel_v / cos(curr_angle) - a.g;
% 
%     measuredState = [measured_accel_h measured_accel_h measured_accel_v measured_gyro measured_gyro measured_gyro];
%     corrected_state = compensateError(measuredState, specs, time_i);
% 
%     corrected_signals(i, :) = corrected_state;
% end
% 
% fprintf('\nFinished Corrected Signal\n')
% 
% % Initial States
% p0 =  a.d(tspan(1))+a.pr_d;   % Platform position [m]
% p_dot0 = a.d_dot(tspan(1));   % Platform velocity [m/s]
% pr_err_accum0 = 0;            % Integral of relative position error [m*s]
% pm0 = p0;                     % Platform inetegrated position [m]
% pm_dot = p_dot0;              % Platform integrated velocity [m/s]
% pm_ddot = a.d_ddot(tspan(1)); % Platform measured acceleration [m*s^-2]
% p_theta0 = a.real_ang(tspan(1)); % Platform inertial angle [deg]
% 
% s0 = [p0; p_dot0; pr_err_accum0; pm0; pm_dot; pm_ddot];
% 
% fprintf('\nStarting Integration WITH Sensor Error.')
% fprintf("\nTime: ")
% 
% % Integrate corrected gyro signal
% 
% control_dynamics = @(t_i, state) corrected_signals((floor(t_i./dt)+1), :)';
% 
% [t_corr, int_sig_corr]= rk4_solver(control_dynamics, tspan, s0, dt);
% 
% fprintf('\nFinished Integration WITH Sensor Error.\n')
% 
% figure
% plot(t_error, corrected_signals(:, 1))
% hold on
% plot(t_error, corrected_signals(:, 3))
% hold on
% plot(t, a.real_accel(t))
% xlabel('Time (sec)')
% ylabel('Acceleration (m/s^2)')
% legend('Corrected Sensor Error (horizontal)', 'WITH Sensor Error (vertical)', 'NO Sensor Error')
% title('Measured Accelerometer Signal')
% % ylim([-5 5])
% 
% figure
% plot(t_error, corrected_signals(:, 4))
% hold on
% plot(t, a.real_ang_rate(t))
% xlabel('Time (sec)')
% ylabel('Angular Rate (deg/s)')
% legend('Corrected Sensor Error', 'NO Sensor Error')
% title('Measured Gyroscope Signal')
% % ylim([-5 5])
% 
% figure
% plot(t_corr, int_sig_corr(:, 1))
% hold on
% plot(t_corr, int_sig_corr(:, 3))
% hold on
% plot(t, a.real_vel(t))
% xlabel('Time (sec)')
% ylabel('Velocity (m/s)')
% legend('Corrected Sensor Error (horizontal)', 'Corrected Sensor Error (vertical)', 'NO Sensor Error')
% title('Integrated Velocity Signal')
% % ylim([-5 5])
% 
% figure
% plot(t_corr, int_sig_corr(:, 4))
% hold on
% plot(t, a.real_ang(t))
% xlabel('Time (sec)')
% ylabel('Angle (deg)')
% legend('Corrected Sensor Error', 'NO Sensor Error')
% title('Integrated Angle Signal')
% % ylim([-5 5])

%% Control Law Drift Compensation

fprintf('\nStarting Integration WITH Sensor Compensation Control Law.')
fprintf("\nTime: ")

control_dynamics = @(t_i, state) IMUDriftCorrection(t_i, a, state);

vel0_x = 0;
vel0_y = 0;
vel0_z = 0;
theta0 = 0;
phi0 = 0;
psi0 = 0;
p0_x = 0;
p0_y = 0;
p0_z = 0;
a_err0_x = 0;
a_err0_y = 0;
a_err0_z = 0;

s0 = [vel0_x vel0_y vel0_z theta0 phi0 psi0 p0_x p0_y p0_z a_err0_x a_err0_y a_err0_z];

[t_corr, int_state_corr]= rk4_solver(control_dynamics, tspan, s0, dt);

fprintf('\nFinished Integration WITH Sensor Compensation Control Law.\n')


%% Plotting

figure
subplot(3,1,1)
plot(t, a.theta(t))
hold on
plot(t_corr, int_state_corr(:, 4))
xlabel('Time (sec)')
ylabel('Angle (deg)')
legend('NO Sensor Error', 'Controlled Error')
title('Controlled Integrated Angle')
% ylim([-5 5])

subplot(3,1,2)
plot(t, a.phi(t))
hold on
plot(t_corr, int_state_corr(:, 5))
xlabel('Time (sec)')
ylabel('Angle (deg)')
legend('NO Sensor Error', 'Controlled Error')
title('Controlled Integrated Angle')
% ylim([-5 5])

subplot(3,1,3)
plot(t, a.psi(t))
hold on
plot(t_corr, int_state_corr(:, 6))
xlabel('Time (sec)')
ylabel('Angle (deg)')
legend('NO Sensor Error', 'Controlled Error')
title('Controlled Integrated Angle')
% ylim([-5 5])


figure
subplot(3,1,1)
plot(t, a.real_vel_xI(t))
hold on
plot(t_corr, int_state_corr(:, 1))
xlabel('Time (sec)')
ylabel('Velocity (m/s)')
legend('NO Sensor Error', 'Controlled Error Integrated')
title('Controlled Integrated (vertical) Velocity')
% ylim([-5 5])

subplot(3,1,2)
plot(t, a.real_vel_yI(t))
hold on
plot(t_corr, int_state_corr(:, 2))
xlabel('Time (sec)')
ylabel('Velocity (m/s)')
legend('NO Sensor Error', 'Controlled Error Integrated')
title('Controlled Integrated (vertical) Velocity')
% ylim([-5 5])

subplot(3,1,3)
plot(t, a.real_vel_zI(t))
hold on
plot(t_corr, int_state_corr(:, 3))
xlabel('Time (sec)')
ylabel('Velocity (m/s)')
legend('NO Sensor Error', 'Controlled Error Integrated')
title('Controlled Integrated (vertical) Velocity')
% ylim([-5 5])

%% Plotting Acceleration


accel_m_controlled_x = diff(int_state_corr(:, 1))/dt;
accel_m_controlled_y = diff(int_state_corr(:, 2))/dt;
accel_m_controlled_z = diff(int_state_corr(:, 3))/dt;

figure
subplot(3,1,1)
plot(t, a.real_accel_xI(t))
hold on
plot(t(2:end), accel_m_controlled_x)
xlabel('Time (sec)')
ylabel('Acceleration (m/s^2)')
legend('NO Sensor Error', 'Controlled Error Integrated')
title('Controlled Measured Acceleration X')
% ylim([-5 5])

subplot(3,1,2)
plot(t, a.real_accel_yI(t))
hold on
plot(t(2:end), accel_m_controlled_y)
xlabel('Time (sec)')
ylabel('Acceleration (m/s^2)')
legend('NO Sensor Error', 'Controlled Error Integrated')
title('Controlled Measured Acceleration Y')
% ylim([-5 5])

subplot(3,1,3)
plot(t, a.real_accel_zI(t) + 9.81)
hold on
plot(t(2:end), accel_m_controlled_z)
xlabel('Time (sec)')
ylabel('Acceleration (m/s^2)')
legend('NO Sensor Error', 'Controlled Error Integrated')
title('Controlled Measured Acceleration Z')
% ylim([-5 5])


function s_dot = IMUDriftCorrection(time_i, a, prev_state)

    % Current states
    vel_x = prev_state(1);
    vel_y = prev_state(2);
    vel_z = prev_state(3);

    theta = prev_state(4);
    phi = prev_state(5);
    psi = prev_state(6);

    p_x = prev_state(7);
    p_y = prev_state(8);
    p_z = prev_state(9);

    theta_err_accum = prev_state(10);
    phi_err_accum = prev_state(11);
    psi_err_accum = prev_state(12);

    specs = a.specs;

    state = [vel_x; vel_y; vel_z; theta; phi; psi];
    state_err_accum = [p_x; p_y; p_z; theta_err_accum; phi_err_accum; psi_err_accum];

    kw = a.kw;
    kt = a.kt;

    state_0 = 0;
    state_dot_0 = 0;
    
    % Get measured acceleration based on real angle
    theta_real = a.theta(time_i);
    phi_real = a.phi(time_i);
    psi_real = a.psi(time_i);

    a_I = [a.real_accel_xI(time_i); a.real_accel_yI(time_i); a.real_accel_zI(time_i)];

    a_S = Rotate_I_S(a_I, theta_real, phi_real, psi_real);

    a_Sx_m = a.measured_accel(time_i, a.o_d_n_a_c, a.biasStabDistAccel, a.biasTempDistAccel, a.accel_drift, a.noiseDistAccel, a_S(1));
    a_Sy_m = a.measured_accel(time_i, a.o_d_n_a_c, a.biasStabDistAccel, a.biasTempDistAccel, a.accel_drift, a.noiseDistAccel, a_S(2));
    a_Sz_m = a.measured_accel(time_i, a.o_d_n_a_c, a.biasStabDistAccel, a.biasTempDistAccel, a.accel_drift, a.noiseDistAccel, a_S(3));
   
    accel_S = [a_Sx_m; a_Sy_m; a_Sz_m];

    theta_d = a.theta_dot_eq(time_i);
    phi_d = a.phi_dot_eq(time_i);
    psi_d = a.psi_dot_eq(time_i);

    theta_dot_m = a.measured_gyro(time_i, a.o_d_n_g_c, a.biasStabDistGyro, a.biasTempDistGyro, a.gyro_drift, a.noiseDistGyro, theta_d);
    phi_dot_m = a.measured_gyro(time_i, a.o_d_n_g_c, a.biasStabDistGyro, a.biasTempDistGyro, a.gyro_drift, a.noiseDistGyro, phi_d);
    psi_dot_m = a.measured_gyro(time_i, a.o_d_n_g_c, a.biasStabDistGyro, a.biasTempDistGyro, a.gyro_drift, a.noiseDistGyro, psi_d);

    accel_I = Rotate_S_I(accel_S, theta, phi, psi);
    accel_I(3) = accel_I(3) + a.g;
    
    % Remove gravity using trig
    % corr_a = AccelRemoveGrav(accel_I, curr_angle, a);
    % 
    % measured_a_h = corr_a(1);
    % measured_a_v = corr_a(3);

    measuredState = [accel_I', theta_dot_m, phi_dot_m, psi_dot_m];
    corrected_state = compensateError(measuredState, specs, time_i);
    
    state_dot_m = corrected_state';

    state_control = state - state_0;

    state_dot_comp = kt .* state_err_accum + kw .* state_control;

    state_dot = state_dot_m + state_dot_0 - state_dot_comp;

    s_dot = zeros(12,1);

    s_dot(1:6) = state_dot;
    s_dot(7:12) = state_control;

end