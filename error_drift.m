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
tspan = [0 10]; % [s]
steps = 2000;
timestep = tspan(1)/steps;
t = linspace(0, tspan(2), steps);

imx_5_specs
% gx5_specs

k = specs.k;
dk = specs.dk;
V_err_0 = specs.V_err_0;
P_err_0 = specs.P_err_0;
accel_resolution = specs.accel_resolution;
accel_samplingRate = specs.accel_samplingRate;
accel_noiseDensity = specs.accel_noiseDensity;
gyro_resolution = specs.gyro_resolution;
gyro_samplingRate = specs.gyro_samplingRate;
gyro_noiseDensity = specs.gyro_noiseDensity;
b_a = specs.b_a;
VRW = specs.VRW;
specs.b_g = 0;
b_g = specs.b_g; 
ARW = specs.ARW;
accel_bandwidth = specs.accel_bandwidth;
gyro_bandwidth = specs.gyro_bandwidth;

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

accel_noise_std = accel_noiseDensity * sqrt(accel_bandwidth); % Noise standard deviation (microg)
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

figure(5)

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

figure(6)

% Resolution/Quantization

gyro_quantized = @(t, real_ang_rate) gyro_resolution*floor(real_ang_rate(t)/gyro_resolution);

subplot(1,3,1)
plot(t,gyro_quantized(t, real_ang_rate))
xlabel('Time (s)')
ylabel('Angular Velocity (deg/s)')
title('Quantized Gyroscope Signal')

% Noise

gyro_noise_std = gyro_noiseDensity * sqrt(gyro_bandwidth); % Noise standard deviation (microg)
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

figure(7)

plot(t, gyro_drift(t, n_g, real_ang_rate))
hold on
plot(t, real_ang_rate(t))
xlabel('Time (s)')
ylabel('Angular Velocity (deg/s)')
title('Drifting Gyroscope Signal')
legend('Drifting Gyro Measurements','Real Gyro Measurements')

%% Combined Realistic Accel/Gyro Signals

% This assumes quantization, followed by noise, followed by drift

figure(8)

quant_noise_accel = @(t, i, real_accel) accel_quantized(t, real_accel) + accel_randomNoise(i)';

drift_error_accel_vert = @(t, n_a, real_accel, theta_err, i) (1 + k)*quant_noise_accel(t, i, real_accel) + b_a + n_a(t) + a.g*(1-cos(theta_err(t)));
drift_error_accel_horz = @(t, n_a, real_accel, theta_err, i) (1 + k)*quant_noise_accel(t, i, real_accel) + b_a + n_a(t) + a.g*sin(theta_err(t));
accel_measured_vert = @(t, n_a, real_accel, theta_err, i) drift_error_accel_vert(t, n_a, real_accel, theta_err, i);
accel_measured_horz = @(t, n_a, real_accel, theta_err, i) drift_error_accel_horz(t, n_a, real_accel, theta_err, i);

quant_noise_gyro = @(t, i, real_ang_rate) gyro_quantized(t, real_ang_rate) + gyro_randomNoise(i)';
drift_error_gyro = @(t, n_g, real_ang_rate, theta_err, i) (1 + k)*quant_noise_gyro(t, i, real_ang_rate) + b_g + n_g(t);
gyro_measured = @(t, n_g, real_ang_rate, theta_err, i) drift_error_gyro(t, n_g, real_ang_rate, theta_err, i);

subplot(1,2,1)
% plot(t, accel_measured_vert(t, n_a, real_accel, theta_err, i))
% hold on
plot(t, accel_measured_horz(t, n_a, real_accel, theta_err, i))
hold on
plot(t, real_accel(t))
xlabel('Time (s)')
ylabel('Acceleration (m/s^2)')
title('Real Accelerometer Signal')

subplot(1,2,2)
plot(t, gyro_measured(t, n_g, real_ang_rate, theta_err, i))
hold on
plot(t, real_ang_rate(t))
xlabel('Time (s)')
ylabel('Angular Velocity (deg/s)')
title('Real Gyroscope Signal')

StatesOverTime_measured = zeros(length(t), 6);
StatesOverTime_measured(:, 1) = accel_measured_horz(t, n_a, real_accel, theta_err, i);
StatesOverTime_measured(:, 2) = accel_measured_horz(t, n_a, real_accel, theta_err, i);
StatesOverTime_measured(:, 3) = accel_measured_vert(t, n_a, real_accel, theta_err, i);
StatesOverTime_measured(:, 4) = gyro_measured(t, n_a, real_ang_rate, theta_err, i);
StatesOverTime_measured(:, 5) = gyro_measured(t, n_a, real_ang_rate, theta_err, i);
StatesOverTime_measured(:, 6) = gyro_measured(t, n_a, real_ang_rate, theta_err, i);
StatesOverTime_measured(1,:) = [];


%% Error Compensation

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

figure(9)

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

%% Finding Accuracy of Corrected Data

real_accel_data = real_accel(t)';
real_gyro_data = real_ang_rate(t)';

figure(11)

specs.beta = beta;
fitSolutionAccel = LeastSquares(t', StatesOverTime_corrected(:, 1), specs);
fitSolutionGyro = LeastSquares(t', StatesOverTime_corrected(:, 4), specs);

error_accel = abs(fitSolutionAccel - real_accel_data);
error_gyro = abs(fitSolutionGyro - real_gyro_data);

subplot(1,2,1)
plot(t, fitSolutionAccel, color='red')
hold on
plot(t, real_accel(t), color='black', LineWidth=0.5)
hold on
plot(t, error_accel)

xlabel('Time (s)')
ylabel('Acceleration (m/s^2)')
title('Fit Curve vs Real Curve - Accelerometer')
legend('Fit Line Accel', 'Expected Accel Data', 'Error')

subplot(1,2,2)
plot(t, fitSolutionGyro, color='red')
hold on
plot(t, real_ang_rate(t), color='black', LineWidth=0.5)
hold on
plot(t, error_gyro)

xlabel('Time (s)')
ylabel('Angular Velocity (deg/s)')
title('Fit Curve vs Real Curve - Gyroscope')
legend('Fit Line Gyro', 'Expected Gyro Data', 'Error')

figure(12)

sz = 3;

subplot(1,2,1)
scatter(t, error_accel, sz, "filled")

xlabel('Time (s)')
ylabel("Error (m/s^2)")
title('Accelerometer Error Over Time')
ylim([0 0.15])

subplot(1,2,2)
scatter(t, error_gyro, sz, "filled")

xlabel('Time (s)')
ylabel("Error (deg/s)")
title('Gyroscope Error Over Time')
ylim([0 2])

%% Finding Accuracy

accel_accuracy = 100*(1-min(abs(error_accel), abs(real_accel_data)) ./ max(abs(error_accel), abs(real_accel_data)));
gyro_accuracy = 100*(1-min(abs(error_gyro), abs(real_gyro_data)) ./ max(abs(error_gyro), abs(real_gyro_data)));

accel_LSB = error_accel / accel_resolution;
gyro_LSB = error_gyro / gyro_resolution;

fitAccelAccuracy = LeastSquares(t', accel_accuracy, specs);
fitGyroAccuracy = LeastSquares(t', gyro_accuracy, specs);

AccelAccuracy = mean(fitAccelAccuracy);
GyroAccuracy = mean(fitGyroAccuracy);

figure(13)

subplot(1,2,1)
plot(t, accel_accuracy)
hold on
plot(t, fitAccelAccuracy)
xlabel("Time (s)")
ylabel("Accuracy (%)")
title("Accelerometer Accuracy (with correction)")
legend("Correction Accuracy", append("General Accuracy of ", num2str(AccelAccuracy),"%"))

subplot(1,2,2)
plot(t, gyro_accuracy)
hold on
plot(t, fitGyroAccuracy)
xlabel("Time (s)")
ylabel("Accuracy (%)")
title("Gyroscope Accuracy (with correction)")
legend("Correction Accuracy", append("General Accuracy of ", num2str(GyroAccuracy),"%"))

figure(14)

subplot(1,2,1)
yyaxis left
plot(t, accel_accuracy)
hold on
plot(t, fitAccelAccuracy, color="black")
hold on
xlabel("Time (s)")
ylabel("Accuracy (%)")
title("Accelerometer Accuracy (with correction)")

yyaxis right
plot(t, real_accel(t))
ylim([-0.4 0.4])
ylabel('Deck Disturbance Acceleration (m*s^-2)')
legend("Correction Accuracy", append("General Accuracy of ", num2str(AccelAccuracy),"%"), "Deck Disturbance")

subplot(1,2,2)
yyaxis left
plot(t, gyro_accuracy)
hold on
plot(t, fitGyroAccuracy, color="black")
hold on
xlabel("Time (s)")
ylabel("Accuracy (%)")
title("Gyroscope Accuracy (with correction)")

yyaxis right
plot(t, real_ang_rate(t))
ylim([-25 25])
ylabel('Deck Disturbance Angular Rate (deg/s)')
legend("Correction Accuracy", append("General Accuracy of ", num2str(GyroAccuracy),"%"), "Deck Disturbance")


figure(15)

subplot(1,2,1)
yyaxis left
plot(t, accel_LSB)
hold on
xlabel("Time (s)")
ylabel("Precision Error (Bits)")
title("Number of Acceleration Bits of Error (with correction)")

yyaxis right
plot(t, real_accel(t))
ylim([-0.4 0.4])
ylabel('Deck Disturbance Acceleration (m*s^-2)')

subplot(1,2,2)
yyaxis left
plot(t, gyro_LSB)
hold on
xlabel("Time (s)")
ylabel("Precision Error (Bits)")
title("Number of Gyroscope Bits of Error (with correction)")

yyaxis right
plot(t, real_ang_rate(t))
ylim([-25 25])
ylabel('Deck Disturbance Angular Rate (deg/s)')

%% Find new position/velocity error from acceleration error

accel_reduction = (100 - AccelAccuracy)/100;
gyro_reduction = (100 - GyroAccuracy)/100;

b_a_red = accel_reduction*b_a;
b_g_red = gyro_reduction*b_g;
k_red = accel_reduction*k;
VRW_red = accel_reduction*VRW;
ARW_red = gyro_reduction*ARW;

% Velocity Error
V_err_corrected = @(t, real_accel) V_err_0 + k_red*(real_accel(t) * timestep) + b_a_red*t + VRW_red*sqrt(t) + a.g*(0.5*b_g_red*t.^2 + 2/3*ARW_red*t.^(3/2));

% Position Error
P_err_corrected = @(t, real_vel) P_err_0 + k_red*(real_vel(t) * timestep) + V_err_0*t + 0.5*b_a_red*t.^2 + 2/3*VRW_red*t.^(3/2) + a.g*(1/6*b_g_red*t.^3 + 4/15*ARW_red*t.^(5/2));

% Angular Error
theta_err_corrected = @(t) b_g_red*t + ARW_red*sqrt(t);

figure(16)

subplot(1,3,1)
plot(t, V_err_corrected(t, real_accel))
xlabel("Time (s)")
ylabel("Velocity Error (m/s)")
title("Velocity Error over Time (corrected)")

subplot(1,3,2)
plot(t, P_err_corrected(t, real_vel))
xlabel("Time (s)")
ylabel("Position Error (m)")
title("Position Error over Time (corrected)")


subplot(1,3,3)
plot(t, theta_err_corrected(t))
xlabel("Time (s)")
ylabel("Angular Error (deg)")
title("Angular Error over Time (corrected)")

%% Functions

function state = compensateError(measuredState, specs, time)
    k = specs.k;
    dk = specs.dk;
    b_a = specs.b_a;
    b_g = specs.b_g;
    ARW = specs.ARW;
    VRW = specs.VRW;
    accel_noiseDensity = specs.accel_noiseDensity;
    gyro_noiseDensity = specs.gyro_noiseDensity;
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
    
    % Bias

    ACCEL_BIAS = [b_a b_a b_a]';
    GYRO_BIAS = [b_g b_g b_g]';

    % Noise

    ACCEL_NOISE = [accel_noiseDensity accel_noiseDensity accel_noiseDensity]';
    GYRO_NOISE = [gyro_noiseDensity gyro_noiseDensity gyro_noiseDensity]';

    theta_err = b_g*time + ARW*sqrt(time);

    a_adjusted_x_y = measured_accel(1) - ACCEL_BIAS(1) - g*sin(theta_err);
    a_adjusted_z = measured_accel(3) - ACCEL_BIAS(3) - g*(1-cos(theta_err));
    a_adjusted = [a_adjusted_x_y; a_adjusted_x_y; a_adjusted_z];
    
    A_FIX = inv([1+S_x+dS_x  M_xy       M_xz
                M_yx         1+S_y+dS_y M_yz
                M_zx         M_zy       1+S_z+dS_z]);

    corrected_a = A_FIX * a_adjusted;


    G_DEP_BIAS = [B_gx  0     0
                  0     B_gy  0
                  0     0     B_gz];

    g_adjusted = measured_gyro - GYRO_BIAS - G_DEP_BIAS*corrected_a;
    
    corrected_g = A_FIX * g_adjusted;

    state(1:3) = corrected_a';
    state(4:6) = corrected_g';
    
end

function soln = LeastSquares(time, Datapoints, specs)
    beta = specs.beta;
    time_size = size(time);
    data_size = size(Datapoints);
    if time_size(1) ~= data_size(1)
        disp('Check raw data and time sizes.')
    end

    A = zeros(time_size(1), 2);
    A(:,1) = sin(beta * time);
    A(:,2) = ones(time_size(1), 1);
    
    v = pinv(A)*Datapoints;
    a = v(1);
    b = v(2);
    soln = a*sin(beta * time)+b;
end