clc; clear; close all;

simulationParameters;

% Redefining acceleration/gyro curves for clarity

real_pos = @(t) alpha*sin(beta*t) + hdeck;
real_vel = @(t) beta*alpha*cos(beta*t);
real_accel = @(t) -beta^2*alpha*sin(beta*t); % [m*s^-2]
real_ang_rate = @(t) 180/pi*(-(alpha*beta^2*sin(beta*t))./(alpha^2*beta^2*(cos(beta*t).^2)+1)); % [deg/s]

%% Sensor Specs and Error

% ARW/VRW is our noise measurement
% Bias is our drift measurement, note:
%   measurement of how the bias will drift during operation over time at a constant temperature
% Resolution doesn't matter, as per this:
%   https://www.vectornav.com/resources/inertial-navigation-primer/specifications--and--error-budgets/specs-imuspecs

% Sensor
%   IMX-5: https://docs.inertialsense.com/datasheets/IMX-5_IMU_AHRS_GNSS-INS_Datasheet.pdf

% Simulation time
tspan = [0 300]; % [s]
steps = 2000;
timestep = tspan(1)/steps;
t = linspace(0, tspan(2), steps);

k = 0; % Scale Factor Error, measured as percentage FSR
dk = 0.02; % Scale Factor Nonlinearity, %FS
V_err_0 = 0; % Initial Velocity Error
P_err_0 = 0; % Initial Position Error

accel_resolution = 0.122 / 1000 * a.g; % m/s
accel_samplingRate = 4000; % Hz
accel_noiseDensity = 60 * 10^-6 * a.g; % m/s^2/sqrt(Hz)

gyro_resolution = 0.0076; % deg/s
gyro_samplingRate = 8000; % Hz
gyro_noiseDensity = 5 * 10^-3; % dps/sqrt(Hz)

% Accel Specs
b_a = 0.019; % Time Varying Bias (mg)
VRW = 0.02 / 60; % Velocity Random Walk (m/s/sqrt(s))

% Gyro Specs
b_g = 1.5 / 60; % Time Varying Bias (deg/s)
ARW = 0.16 / 60; % Angle Random Walk (deg/sqrt(s))

% b_a = 0;
% b_g = 0;
% VRW = 0;
% ARW = 0;

% Velocity Error
V_err = @(t, real_accel) V_err_0 + k*(real_accel(t) * timestep) + b_a*t + VRW*sqrt(t) + a.g*(0.5*b_g*t.^2 + 2/3*ARW*t.^(3/2));

% Position Error
P_err = @(t, real_vel) P_err_0 + k*(real_vel(t) * timestep) + V_err_0*t + 0.5*b_a*t.^2 + 2/3*VRW*t.^(3/2) + a.g*(1/6*b_g*t.^3 + 4/15*ARW*t.^(5/2));

% Angular Error
theta_err = @(t) b_g*t + ARW*sqrt(t);

figure(1)

subplot(1,2,1)
plot(t, real_vel(t))
xlabel("Time (sec)")
ylabel("Velocity (m/s)")
title("Deck Velocity over Time")

subplot(1,2,2)
plot(t, real_pos(t))
xlabel("Time (sec)")
ylabel("Position (m)")
title("Deck Position over Time")


figure(2)

subplot(1,3,1)
plot(t, V_err(t, real_accel))
xlabel("Time (s)")
ylabel("Velocity Error (m/s)")
title("Velocity Error over Time")

subplot(1,3,2)
plot(t, P_err(t, real_vel))
xlabel("Time (s)")
ylabel("Position Error (m)")
title("Position Error over Time")


subplot(1,3,3)
plot(t, theta_err(t))
xlabel("Time (s)")
ylabel("Angular Error (deg)")
title("Angular Error over Time")

%% Plotting Real Accel/Gyro Signals

figure(3)

subplot(1,2,1)
plot(t, real_accel(t))
xlabel('Time (s)')
ylabel('Acceleration (m/s^2)')
title('Deck (Real) Acceleration')

subplot(1,2,2)
plot(t, real_ang_rate(t))
xlabel('Time (s)')
ylabel('Angular Rate (deg/s)')
title('Deck (Real) Angular Rate')


%% Creating Realistic Accelerometer Signal

figure(4)

% Resolution/Quantization

accel_quantized = @(t, real_accel) accel_resolution*floor(real_accel(t)/accel_resolution);

subplot(1,3,1)
plot(t,accel_quantized(t, real_accel))
xlabel('Time (s)')
ylabel('Acceleration (m/s^2)')
title('Quantized Accelerometer Signal')

% Noise

accel_noise_std = accel_noiseDensity * sqrt(accel_samplingRate); % Noise standard deviation (microg)
accel_randomNoise = accel_noise_std*randn(length(t), 1);
i = floor(t*steps/tspan(2)) + 1;
i(end) = i(end) - 1;

accel_noise = @(t, i, real_accel) real_accel(t) + accel_randomNoise(i)';

subplot(1,3,2)
plot(t,accel_noise(t,i, real_accel))
xlabel('Time (s)')
ylabel('Acceleration (m/s^2)')
title('Noisy Accelerometer Signal')

% Error Bias

% This is verticle acceleration. For horizontal, replace
% a.g*(1-cos(theta_err)) with a.g*sin(theta_err)
n_a = @(t) 0.5*VRW*t.^(-0.5);
accel_drift_vert = @(t, n_a, real_accel, theta_err) (1 + k)*real_accel(t) + b_a + n_a(t) + a.g*(1-cos(theta_err(t)));
accel_drift_horz = @(t, n_a, real_accel_horz, theta_err) (1 + k)*real_accel_horz(t) + b_a + n_a(t) + a.g*sin(theta_err(t));

subplot(1,3,3)
plot(t,accel_drift_vert(t, n_a, real_accel, theta_err))
hold on
plot(t,accel_drift_horz(t, n_a, real_accel, theta_err))
hold on
plot(t, real_accel(t))
xlabel('Time (s)')
ylabel('Acceleration (m/s^2)')
title('Drifting Accelerometer Signal')
legend('Drifting Accelerometer Vertical Measurements','Drifting Accelerometer Horizontal Measurements', 'Real Acceleration')


%% Creating Realistic Gyroscope Signal

figure(5)

% Resolution/Quantization

gyro_quantized = @(t, real_ang_rate) gyro_resolution*floor(real_ang_rate(t)/gyro_resolution);

subplot(1,3,1)
plot(t,gyro_quantized(t, real_ang_rate))
xlabel('Time (s)')
ylabel('Angular Velocity (deg/s)')
title('Quantized Gyroscope Signal')

% Noise

gyro_noise_std = gyro_noiseDensity * sqrt(gyro_samplingRate); % Noise standard deviation (microg)
gyro_randomNoise = gyro_noise_std*randn(length(t), 1);
i = floor(t*steps/tspan(2)) + 1;
i(end) = i(end) - 1;

gyro_noise = @(t, i, real_ang_rate) real_ang_rate(t) + gyro_randomNoise(i)';

subplot(1,3,2)
plot(t,gyro_noise(t,i, real_ang_rate))
xlabel('Time (s)')
ylabel('Angular Velocity (deg/s)')
title('Noisy Gyroscope Signal')

% Error Bias
n_g = @(t) 0.5*ARW*t.^(-0.5);
gyro_drift = @(t, n_g, real_ang_rate) (1 + k)*real_ang_rate(t) + b_g + n_g(t);

subplot(1,3,3)
plot(t, gyro_drift(t, n_g, real_ang_rate))
hold on
plot(t, real_ang_rate(t))
xlabel('Time (s)')
ylabel('Angular Velocity (deg/s)')
title('Drifting Gyroscope Signal')
legend('Drifting Gyro Measurements','Real Gyro Measurements')

%% Combined Realistic Accel/Gyro Signals

% This assumes quantization, followed by noise, followed by drift

figure(6)

quant_noise_accel = @(t, i, real_accel) accel_quantized(t, real_accel) + accel_randomNoise(i)';

drift_error_accel_vert = @(t, n_a, real_accel, theta_err, i) (1 + k)*quant_noise_accel(t, i, real_accel) + b_a + n_a(t) + a.g*(1-cos(theta_err(t)));
% drift_error_accel_horz = @(t, n_a, real_accel, theta_err, i) (1 + k)*quant_noise_accel(t, i, real_accel) + b_a + n_a(t) + a.g*sin(theta_err(t));
accel_measured_vert = @(t, n_a, real_accel, theta_err, i) drift_error_accel_vert(t, n_a, real_accel, theta_err, i);
% accel_measured_horz = @(t, n_a, real_accel, theta_err, i) drift_error_accel_horz(t, n_a, real_accel, theta_err, i);

quant_noise_gyro = @(t, i, real_ang_rate) gyro_quantized(t, real_ang_rate) + gyro_randomNoise(i)';
drift_error_gyro = @(t, n_g, real_ang_rate, theta_err, i) (1 + k)*quant_noise_gyro(t, i, real_ang_rate) + b_g + n_g(t);
gyro_measured = @(t, n_g, real_ang_rate, theta_err, i) drift_error_gyro(t, n_g, real_ang_rate, theta_err, i);

subplot(1,2,1)
plot(t, accel_measured_vert(t, n_a, real_accel, theta_err, i))
hold on
% plot(t, accel_measured_horz(t, n_a, real_accel, theta_err, i))
% hold on
plot(t, real_accel(t))
xlabel('Time (s)')
ylabel('Acceleration (m/s^2)')
title('Drifting Accelerometer Signal')

subplot(1,2,2)
plot(t, gyro_measured(t, n_g, real_ang_rate, theta_err, i))
hold on
plot(t, real_ang_rate(t))
xlabel('Time (s)')
ylabel('Angular Velocity (deg/s)')
title('Drifting Gyroscope Signal')

StatesOverTime_measured = zeros(length(t), 6);
StatesOverTime_measured(:, 1) = accel_measured_vert(t, n_a, real_accel, theta_err, i);
StatesOverTime_measured(:, 2) = accel_measured_vert(t, n_a, real_accel, theta_err, i);
StatesOverTime_measured(:, 3) = accel_measured_vert(t, n_a, real_accel, theta_err, i);
StatesOverTime_measured(:, 4) = gyro_measured(t, n_a, real_ang_rate, theta_err, i);
StatesOverTime_measured(:, 5) = gyro_measured(t, n_a, real_ang_rate, theta_err, i);
StatesOverTime_measured(:, 6) = gyro_measured(t, n_a, real_ang_rate, theta_err, i);
StatesOverTime_measured(1,:) = [];


%% Error Compensation

StatesOverTime_corrected = zeros(length(t)-1, 6);

constants.k = k;
constants.dk = dk;
constants.b_a = b_a;
constants.b_g = b_g;
constants.ARW = ARW;
constants.VRW = VRW;
constants.accel_noiseDensity = accel_noiseDensity;
constants.gyro_noiseDensity = gyro_noiseDensity;
constants.g = a.g;

errorValues = zeros(length(t)-1, 6);
accumError = zeros(length(t), 6);

for i = 1:(length(t)-1)
    measuredState = StatesOverTime_measured(i, :);
    time = t(i);
    error_compensation = compensateError(measuredState, constants, time);
    StatesOverTime_corrected(i, :) = error_compensation;
    accumError(i+1, :) = error_compensation + accumError(i, :);
end

figure(7)

subplot(1,2,1)
plot(t(2:end), StatesOverTime_measured(:, 1))
hold on
plot(t(2:end), StatesOverTime_corrected(:, 1))
hold on
plot(t, real_accel(t))

xlabel('Time (s)')
ylabel('Acceleration (m/s^2)')
title('Drifting Accelerometer Signal')
legend('Raw Measurement', 'Measurement Correction', 'Expected Calculation')


subplot(1,2,2)
plot(t(2:end), StatesOverTime_measured(:, 4))
hold on
plot(t(2:end), StatesOverTime_corrected(:, 4))
hold on
plot(t, real_ang_rate(t))

xlabel('Time (s)')
ylabel('Angular Velocity (deg/s)')
title('Drifting Gyroscope Signal')
legend('Raw Measurement', 'Measurement Correction', 'Expected Calculation')

figure(8)

subplot(1,2,1)
plot(t(2:end), errorValues(:, 1))
hold on
plot(t, accumError(:, 1))
legend('Error Values', 'Accumulated Error')

subplot(1,2,2)
plot(t(2:end), errorValues(:, 4))
hold on
plot(t, accumError(:, 4))
legend('Error Values', 'Accumulated Error')

close(figure(8))

function state = compensateError(measuredState, constants, time)
    k = constants.k;
    dk = constants.dk;
    b_a = constants.b_a;
    b_g = constants.b_g;
    ARW = constants.ARW;
    VRW = constants.VRW;
    accel_noiseDensity = constants.accel_noiseDensity;
    gyro_noiseDensity = constants.gyro_noiseDensity;
    g = constants.g;
    
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
    
    % Bias

    ACCEL_BIAS = [b_a b_a b_a]';
    GYRO_BIAS = [b_g b_g b_g]';

    % Noise

    ACCEL_NOISE = [accel_noiseDensity accel_noiseDensity accel_noiseDensity]';
    GYRO_NOISE = [gyro_noiseDensity gyro_noiseDensity gyro_noiseDensity]';

    theta_err = b_g*time + ARW*sqrt(time);

    a_adjusted = measured_accel - ACCEL_BIAS - ACCEL_NOISE - g*(1-cos(theta_err));
    
    A_FIX = inv([1+S_x+dS_x  M_xy       M_xz
                M_yx         1+S_y+dS_y M_yz
                M_zx         M_zy       1+S_z+dS_z]);

    corrected_a = A_FIX * a_adjusted;


    G_DEP_BIAS = [B_gx  0     0
                  0     B_gy  0
                  0     0     B_gz];

    g_adjusted = measured_gyro - GYRO_BIAS - GYRO_NOISE - G_DEP_BIAS*corrected_a;
    
    corrected_g = A_FIX * g_adjusted;

    state(1:3) = corrected_a';
    state(4:6) = corrected_g';
    
end
