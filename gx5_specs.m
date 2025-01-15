specs.k = 0.1; % Scale Factor Error, measured as percentage FSR
specs.dk = 0.04; % Scale Factor Nonlinearity, %FS
specs.V_err_0 = 0; % Initial Velocity Error
specs.P_err_0 = 0; % Initial Position Error

specs.accel_resolution = 0.002 / 1000 * a.g; % m/s
specs.accel_samplingRate = 1000; % Hz
specs.accel_noiseDensity = 20 * 10^-6 * a.g; % m/s^2/sqrt(Hz)
specs.accel_bandwidth = 225; % Hz

specs.gyro_resolution = 0.003; % deg/s
specs.gyro_samplingRate = 4000; % Hz
specs.gyro_noiseDensity = 5 * 10^-3; % dps/sqrt(Hz)
specs.gyro_bandwidth = 250; % Hz

% Accel Specs
specs.b_a = 0.04 / 1000 * 9.8; % Time Varying Bias (m/s^-2)
specs.VRW = 0; % Velocity Random Walk (m/s/sqrt(s))

% Gyro Specs
specs.b_g = 8 / 3600; % Time Varying Bias (deg/s)
specs.ARW = 0; % Angle Random Walk (deg/sqrt(s))
