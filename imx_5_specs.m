% specs.k = 0.1; % Scale Factor Error, measured as percentage FSR
specs.k = 0.1; % Scale Factor Error, measured as percentage FSR
specs.dk = 0.02; % Scale Factor Nonlinearity, %FS

specs.V_err_0 = 0; % Initial Velocity Error
specs.P_err_0 = 0; % Initial Position Error

specs.accel_resolution = 0.122 / 1000 * a.g; % m/s
specs.accel_samplingRate = 4000; % Hz
specs.accel_noiseDensity = 60 * 10^-6 * a.g; % m/s^2/sqrt(Hz)
specs.accel_bandwidth = 218; % Hz
specs.accel_temp_bias = 3.7 / 1000 * a.g; % m/s^2 RMS

specs.gyro_resolution = 0.0076 / 180*pi; % rad/s
specs.gyro_samplingRate = 8000; % Hz
specs.gyro_noiseDensity = 5 * (10^-3) / 180*pi; % rad/s/sqrt(Hz)
specs.gyro_bandwidth = 250; % Hz
specs.gyro_temp_bias = 0.3 / 180*pi; % rad/s RMS

% Accel Specs
specs.b_a = 0.019 / 1000 * a.g; % Time Varying Bias (m/s^-2)
specs.VRW = 0.02 / 60; % Velocity Random Walk (m/s/sqrt(s))

% Gyro Specs
specs.b_g = 1.5 / 3600 / 180*pi; % Time Varying Bias (rad/s)
specs.ARW = 0.16 / 60 / 180*pi; % Angle Random Walk (rad/sqrt(s))
