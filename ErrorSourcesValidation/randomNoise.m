close all; clear; clc;

rng(1,"twister");
s = RandStream('mt19937ar','Seed',15,'NormalTransform','Polar')

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

%% Sensor Specs and Error

% ARW/VRW is our noise measurement
% Bias is our drift measurement, note:
%   measurement of how the bias will drift during operation over time at a constant temperature
% Resolution doesn't matter, as per this:
%   https://www.vectornav.com/resources/inertial-navigation-primer/specifications--and--error-budgets/specs-imuspecs

% Sensor
%   IMX-5: https://docs.inertialsense.com/datasheets/IMX-5_IMU_AHRS_GNSS-INS_Datasheet.pdf

% See following page for more information on spec relationships
% https://www.analog.com/en/resources/analog-dialogue/articles/low-noise-feedback-control.html

% Simulation time
startTime = 0;
finishTime = 120;
tspan = [startTime finishTime]; % [s]

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

imu_accel_sampling_rate = 2*specs.accel_bandwidth;
imu_gyro_sampling_rate = 2*specs.gyro_bandwidth;

imu_rate = max(imu_accel_sampling_rate, imu_gyro_sampling_rate);

dt = 1/imu_rate;  % [s]
t = (tspan(1):dt:tspan(2))';

accel_noise_std = specs.accel_noiseDensity * sqrt(imu_rate); % Noise standard deviation (m/s^2)
gyro_noise_std = specs.gyro_noiseDensity * sqrt(imu_rate); % Noise standard deviation (dps)

% Plotting noise distributions

a.noiseDistAccel = @(t) accel_noise_std*randn(length(t),1);
a.noiseDistGyro = @(t) gyro_noise_std*randn(length(t),1);

% Noisy signal
a.noisyAccelCurve = @(t, real_accel) real_accel(t) + accel_noise_std*randn(length(t),1);
a.noisyGyroCurve = @(t, real_ang_rate) real_ang_rate(t) + gyro_noise_std*randn(length(t),1);

figure
subplot(1,2,1)
plot(t, a.noisyAccelCurve(t, a.real_accel)/a.g)
xlabel('Time (s)')
ylabel('Noise Offset (g)')
title('Accelerometer Noisy Curve')
subplot(1,2,2)
xlabel('Time (s)')
plot(t, a.noisyGyroCurve(t, a.real_ang_rate))
ylabel('Noise Offset (dps)')
title('Gyroscope Noisy Curve')

% Columns represent each variable, so [accel gyro]
% Rows represents timesteps, max limit of 5 timesteps kept track of
integrator_type = 1;
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

%% Show Random Walk

runs = 5;

accel_bound = @(t) VRW * sqrt(t);
gyro_bound = @(t) ARW * sqrt(t);

sigma_accel = VRW/3 * sqrt(imu_rate);
sigma_gyro = ARW/3 * sqrt(imu_rate);

white_noise_accel = sigma_accel * randn(length(t), runs);
drift_walk_accel = cumsum(white_noise_accel, 1) / imu_rate;

white_noise_gyro = sigma_gyro * randn(length(t), runs);
drift_walk_gyro = cumsum(white_noise_gyro, 1) / imu_rate;

figure
plot(t, drift_walk_accel)
hold on
plot(t, gyro_bound(t))
hold on
plot(t, -gyro_bound(t))
xlabel("Time (sec)")
ylabel("Angle (deg))")
title("Bounded Velocity Random Walk")

figure
plot(t, drift_walk_accel)
hold on
plot(t, accel_bound(t))
hold on
plot(t, -accel_bound(t))
xlabel("Time (sec)")
ylabel("Angle (deg)")
title("Bounded Angle Random Walk")


%% Functions

function output_vec = insertVector(originalVector, addVector)

    % Assumes the addVector is a row
    % Removes last row of originalVector
    % Addes as first row of originalVector

    output_vec = [addVector; originalVector(1:end-1, :)];

end



