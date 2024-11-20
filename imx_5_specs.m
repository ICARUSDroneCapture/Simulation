k = 0.1; % Scale Factor Error, measured as percentage FSR
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
b_a = 0.019 / 1000 * 9.8; % Time Varying Bias (m/s^-2)
VRW = 0.02 / 60; % Velocity Random Walk (m/s/sqrt(s))

% Gyro Specs
b_g = 1.5 / 3600; % Time Varying Bias (deg/s)
ARW = 0.16 / 60; % Angle Random Walk (deg/sqrt(s))