rng(1,"twister");

% Use same random seed

set(groot,'DefaultLineLineWidth',1)

% ---------------------------- Define sim time ----------------------------

% Simulation time
startTime = 0;
finishTime = 60*3 + 60;
tspan = [startTime finishTime]; % [s]

a.finishCalibrationTime = 60; % seconds

freq = 160; % Hz
dt = 1/freq;
a.dt = dt;
t = (tspan(1):dt:tspan(2))';
t_count = length(t);
indeces = @(t) floor(t/dt)+1;

% -------------------------------------------------------------------------

simulationParameters

%% Sensor Frame Accelerations

a_I_over_time = [a.real_accel_xI(t), a.real_accel_yI(t), a.real_accel_zI(t)];
omega_over_time = [a.theta_dot(t), a.phi_dot(t), a.psi_dot(t)];

theta_over_time = a.theta(t);
phi_over_time = a.phi(t);
psi_over_time = a.psi(t);

angles_over_time = [theta_over_time, phi_over_time, psi_over_time];

a_S_true = zeros(length(t), 3);
a_S_measured = zeros(length(t), 3);
a_I_ref = zeros(length(t), 3);

for i = 1:length(t)

    t_i = t(i);

    theta_i = theta_over_time(i);
    phi_i = phi_over_time(i);
    psi_i = psi_over_time(i);

    a_I_i = a_I_over_time(i, :)';
    
    % Convert inertial frame accelerations to sensor frame
    a_S_i = Rotate_I_S(a_I_i, theta_i, phi_i, psi_i);
    
    % Get measured sensor frame accelerations
    a_S_m = a.measured_accel_3D(a, t_i, a_S_i);
    
    % Convert Measured Accelerations back to Inertial without any Correction
    a_I_i = Rotate_S_I(a_S_m, theta_i, phi_i, psi_i);

    a_S_true(i, :) = a_S_i';
    a_S_measured(i, :) = a_S_m';
    a_I_ref(i, :) = a_I_i';

end

%% Control Law Drift Compensation

fprintf('\nStarting Integration WITH Sensor Compensation Control Law.')
fprintf("\nTime: ")

control_dynamics = @(t_i, state) IMUDriftCorrection(t_i, a, state, a.finishCalibrationTime, a_I_over_time, omega_over_time, angles_over_time);

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

% Get acceleration
accel_m_controlled_x = diff(int_state_corr(:, 1))/dt;
accel_m_controlled_y = diff(int_state_corr(:, 2))/dt;
accel_m_controlled_z = diff(int_state_corr(:, 3))/dt;

% Get angular velocity
gyro_m_controlled_theta = diff(int_state_corr(:, 4))/dt;
gyro_m_controlled_phi = diff(int_state_corr(:, 5))/dt;
gyro_m_controlled_psi = diff(int_state_corr(:, 6))/dt;


%% Plot Verification

figure
subplot(3,1,1)
plot(t, a.real_accel_xI(t))
hold on
plot(t, a_S_measured(:, 1))
xlabel('Time (sec)')
ylabel('Acceleration (m/s^2)')
title('Measured Sensor Acceleration X (m/s^2)')
legend('Inertial Frame Acceleration', 'Sensor Frame Acceleration')

subplot(3,1,2)
plot(t, a.real_accel_yI(t))
hold on
plot(t, a_S_measured(:, 2))
xlabel('Time (sec)')
ylabel('Acceleration (m/s^2)')
title('Measured Sensor Acceleration Y (m/s^2)')
legend('Inertial Frame Acceleration', 'Sensor Frame Acceleration')

subplot(3,1,3)
plot(t, a.real_accel_zI(t))
hold on
plot(t, a_S_measured(:, 3))
xlabel('Time (sec)')
ylabel('Acceleration (m/s^2)')
title('Measured Sensor Acceleration Z (m/s^2)')
legend('Inertial Frame Acceleration', 'Sensor Frame Acceleration')


figure
subplot(3,1,1)
plot(t, a.real_accel_xI(t))
hold on
plot(t, a_I_ref(:, 1))
xlabel('Time (sec)')
ylabel('Acceleration (m/s^2)')
title('Inertial Acceleration X Converted From Sensor (m/s^2)')
legend('Inertial Frame Acceleration', 'Sensor Frame Acceleration')

subplot(3,1,2)
plot(t, a.real_accel_yI(t))
hold on
plot(t, a_I_ref(:, 3))
xlabel('Time (sec)')
ylabel('Acceleration (m/s^2)')
title('Inertial Acceleration Y Converted From Sensor (m/s^2)')
legend('Inertial Frame Acceleration', 'Sensor Frame Acceleration')

subplot(3,1,3)
plot(t, a.real_accel_zI(t))
hold on
plot(t, a_I_ref(:, 3))
xlabel('Time (sec)')
ylabel('Acceleration (m/s^2)')
title('Inertial Acceleration Z Converted From Sensor (m/s^2)')
legend('Inertial Frame Acceleration', 'Sensor Frame Acceleration')



figure
subplot(3,1,1)
plot(t, a.theta(t))
hold on
plot(t_corr, int_state_corr(:, 4))
xlabel('Time (sec)')
ylabel('Angle (deg)')
legend('NO Sensor Error', 'Controlled Error')
title('Controlled Integrated Angle')
ylim([-0.3 0.3])

subplot(3,1,2)
plot(t, a.phi(t))
hold on
plot(t_corr, int_state_corr(:, 5))
xlabel('Time (sec)')
ylabel('Angle (deg)')
legend('NO Sensor Error', 'Controlled Error')
title('Controlled Integrated Angle')
ylim([-0.3 0.3])

subplot(3,1,3)
plot(t, a.psi(t))
hold on
plot(t_corr, int_state_corr(:, 6))
xlabel('Time (sec)')
ylabel('Angle (deg)')
legend('NO Sensor Error', 'Controlled Error')
title('Controlled Integrated Angle')
ylim([-0.3 0.3])


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

figure
subplot(3,1,1)
plot(t, a.real_accel_xI(t))
hold on
plot(t(2:end), accel_m_controlled_x)
hold on
xline(a.finishCalibrationTime, 'b--')
xlabel('Time (sec)')
ylabel('Acceleration (m/s^2)')
legend('NO Sensor Error', 'Controlled Error Integrated')
title('Controlled Measured Acceleration X')
% ylim([-0.4 0.4])

subplot(3,1,2)
plot(t, a.real_accel_yI(t))
hold on
plot(t(2:end), accel_m_controlled_y)
hold on
xline(a.finishCalibrationTime, 'b--')
xlabel('Time (sec)')
ylabel('Acceleration (m/s^2)')
legend('NO Sensor Error', 'Controlled Error Integrated')
title('Controlled Measured Acceleration Y')
% ylim([-0.4 0.4])

subplot(3,1,3)
plot(t, a.real_accel_zI(t) + a.g)
hold on
plot(t(2:end), accel_m_controlled_z)
hold on
xline(a.finishCalibrationTime, 'b--')
xlabel('Time (sec)')
ylabel('Acceleration (m/s^2)')
legend('NO Sensor Error', 'Controlled Error Integrated')
title('Controlled Measured Acceleration Z')
% ylim([-0.4 0.4]) 

%% Plotting Angular Velocity


gyro_m_controlled_theta = diff(int_state_corr(:, 4))/dt;
gyro_m_controlled_phi = diff(int_state_corr(:, 5))/dt;
gyro_m_controlled_psi = diff(int_state_corr(:, 6))/dt;

figure
subplot(3,1,1)
plot(t, a.theta_dot(t))
hold on
plot(t(2:end), gyro_m_controlled_theta)
hold on
xline(a.finishCalibrationTime, 'b--')
xlabel('Time (sec)')
ylabel('Angular Velocity (rad/s)')
legend('NO Sensor Error', 'Controlled Error Integrated')
title('Controlled Measured Angular Velocity Theta')
% ylim([-5 5])

subplot(3,1,2)
plot(t, a.phi_dot(t))
hold on
plot(t(2:end), gyro_m_controlled_phi)
hold on
xline(a.finishCalibrationTime, 'b--')
xlabel('Time (sec)')
ylabel('Angular Velocity (rad/s)')
legend('NO Sensor Error', 'Controlled Error Integrated')
title('Controlled Measured Angular Velocity Phi')
% ylim([-5 5])

subplot(3,1,3)
plot(t, a.psi_dot(t))
hold on
plot(t(2:end), gyro_m_controlled_psi)
hold on
xline(a.finishCalibrationTime, 'b--')
xlabel('Time (sec)')
ylabel('Angular Velocity (rad/s)')
legend('NO Sensor Error', 'Controlled Error Integrated')
title('Controlled Measured Angular Velocity Psi')
% ylim([-5 5])

%% Functions

function s_dot = IMUDriftCorrection(time_i, a, prev_state, finishCalibrationTime, a_I_over_time, omega_over_time, angles_over_time)

    % Current states
    vel_x = prev_state(1);
    vel_y = prev_state(2);
    vel_z = prev_state(3);

    angle_theta = prev_state(4);
    angle_phi = prev_state(5);
    angle_psi = prev_state(6);

    p_x = prev_state(7);
    p_y = prev_state(8);
    p_z = prev_state(9);

    theta_err_accum = prev_state(10);
    phi_err_accum = prev_state(11);
    psi_err_accum = prev_state(12);

    specs = a.specs;

    state = [vel_x; vel_y; vel_z; angle_theta; angle_phi; angle_psi];
    state_err_accum = [p_x; p_y; p_z; theta_err_accum; phi_err_accum; psi_err_accum];

    kw = a.kw;
    kt = a.kt;

    state_0 = 0;
    state_dot_0 = 0;
    
    % Get real inertial accelerations
    % a_I = [a.real_accel_xI(time_i); a.real_accel_yI(time_i); a.real_accel_zI(time_i)];
    % ang_rate_i = [a.theta_dot(time_i); a.phi_dot(time_i); a.psi_dot(time_i)];

    a_I = a_I_over_time(floor(time_i./a.dt)+1, :)';
    ang_rate_i = omega_over_time(floor(time_i./a.dt)+1, :)';
    
    % Get real angles
    angles_i = angles_over_time(floor(time_i./a.dt)+1, :);

    theta_real = angles_i(1);
    phi_real = angles_i(2);
    psi_real = angles_i(3);

    % Get real sensor frame accelerations
    a_S = Rotate_I_S(a_I, psi_real, theta_real, phi_real);
    
    accel_m = a.measured_accel_3D(a, time_i, a_S);
    gyro_m = a.measured_gyro_3D(a, time_i, ang_rate_i);

    % accel_m = a_S;
    % gyro_m = ang_rate_i;

    theta_use = theta_real;
    phi_use = phi_real;
    psi_use = psi_real;

    if time_i > finishCalibrationTime
        theta_use = angle_theta;
        phi_use = angle_phi;
        psi_use = angle_psi;
    end
    
    accel_I = Rotate_S_I(accel_m, psi_use, theta_use, phi_use);

    accel_state = [accel_I(1) accel_I(2) accel_I(3)+a.g];

    measuredState = [accel_state, gyro_m'];
    corrected_state = compensateError(measuredState, specs, time_i);
    
    state_dot_m = corrected_state';

    state_control = state - state_0;

    state_dot_comp = kt .* state_err_accum + kw .* state_control;

    state_dot = state_dot_m + state_dot_0 - state_dot_comp;

    s_dot = zeros(12,1);

    s_dot(1:6) = state_dot;
    s_dot(7:12) = state_control;

end
