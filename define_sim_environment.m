%% Simulation Duration and Timesteps

% Simulation time
startTime = 0;
finishTime = 30;
tspan = [startTime finishTime]; % [s]

finishCalibrationTime = 60; % seconds

% dt = 1/imu_rate;  % [s]
dt = 0.0001;
a.dt = dt;
t = (tspan(1):dt:tspan(2))';
t_count = length(t);
indeces = @(t) floor(t/dt)+1;

%% Sensor Model Aspects

defineSignals
% defineSignalsNoNoise

%% Environment Parameters

alpha = 0.45; % wave amplitdue [m]
hdeck = 1;   % inertial reference deck hight [m] (arbitrary)

% Wave frequency
t_min = 2.5;     % Minimum period [s]
t_max = 9;    % Maximum period [s]

period = 7.5;    % Expected period [s]

period_small = 2; % Small period of 2 seconds for x/y translations
period_big = 7.5; % Big period of 7.5 (maybe 10) seconds for z motion

amp_small = 5;
amp_big = 1.5;

k = 1;
T = period / k;  % Period of deck disturbance [s]
beta = 2*pi/T; % wave frequency [rad/s]
a.beta_min = 1/t_max; % Minimum frequency [Hz]
a.beta_max = 1/t_min; % Maximum frequency [Hz]
 
% Cosine Wave
% a.d = @(t) alpha*cos(beta*t) + hdeck;      % [m]
% a.d_dot = @(t) -beta*alpha*sin(beta*t);    % [m/s]
% a.d_ddot = @(t) -beta^2*alpha*cos(beta*t); % [m*s^-2]

% Square Wave
% N = 3; % Number of terms in the Fourier series
% coefficients = 1:2:(2*N - 1); % Odd harmonics: 1, 3, 5, ..., (2*N-1)
% beta2 = 2*pi/15;
% a.sw = @(t) alpha*((4/pi) * ...
%     sum(arrayfun(@(n) sin(n*beta2*t)/n, coefficients))) + hdeck; % [m]
% a.d = @(t) arrayfun(a.sw, t);
% a.sw_dot = @(t) alpha*((4/pi) * ...
%     sum(arrayfun(@(n) beta2*cos(n*beta2*t), coefficients)));    % [m/s]
% a.d_dot = @(t) arrayfun(a.sw_dot, t);
% a.sw_ddot = @(t) alpha*((4/pi) * ...
%     sum(arrayfun(@(n) -n*beta2^2*sin(n*beta2*t), coefficients))); % [m*s^-2]
% a.d_ddot = @(t) arrayfun(a.sw_ddot, t);

% Stacked sine wave
% a.d = @(t) alpha*(1.5*sin(beta*t/6) + 0.75*sin(beta*t)) + hdeck; % [m]
% a.d_dot = @(t) alpha*(0.25*beta*cos(beta*t/6) ...
%                                     + 0.75*beta*cos(beta*t)); % [m/s]
% a.d_ddot = @(t) -alpha*(0.0417*beta^2*sin(beta*t/6) ...
%                                     + 0.75*beta^2*sin(beta*t)); % [m*s^-2]

%% Sensor Offset Correction

% For quick gain scaling
scale_w = 1;
scale_t = 1;

% % Raising t_min increasing how quichkly it gets there, but introduces
% more offset (moves offset right, increases overshoot)
% % Decreasing t_min makes it get there slower, but has less offset

% % Raising t_max decreases how quickly the average is brought to nominal,
% % reduces offset
% % Decreasing t_max increases how quickly the long term average is brough to
% % nominal, introduces more offset (moves offset left, decreases overshoot)

% -------------------------------------------------------------------------
%                            SET GYRO GAINS
% -------------------------------------------------------------------------


% ------------ solid gains for 30 second average wave period --------------
% 
% t_min = 4;
% t_max = 6000;


% ------- solid gains for 7.5 second wave period, big amplitude -------

t_min_big_ampl = 5.5;
t_max_big_ampl = 32;

% ------- solid gains for 7.5 second wave period, small amplitude ---------

t_min_small_ampl = 6;
t_max_small_ampl = 36;

% --------- solid gains for 2 second wave period, small amplitude ---------
% 
% t_min = 1.95;
% t_max = 3.05;

% Generalized frequency bounds, these allow for slower convergence but
% better long term tracking

t_min_gyro_x = t_min_small_ampl;
t_min_gyro_y = t_min_big_ampl;
t_min_gyro_z = t_min_small_ampl;

t_max_gyro_x = t_max_small_ampl;
t_max_gyro_y = t_max_big_ampl;
t_max_gyro_z = t_max_small_ampl;

% Gyro frequencies
a.beta_min_x_gyro = 1/t_max_gyro_x;
a.beta_min_y_gyro = 1/t_max_gyro_y;
a.beta_min_z_gyro = 1/t_max_gyro_z;

a.beta_max_x_gyro = 1/t_min_gyro_x;
a.beta_max_y_gyro = 1/t_min_gyro_y;
a.beta_max_z_gyro = 1/t_min_gyro_z;


% -------------------------------------------------------------------------
%                            SET ACCEL GAINS
% -------------------------------------------------------------------------


% ------------ solid gains for 30 second average wave period --------------
% 
% t_min = 4;
% t_max = 6000;


% ------- solid gains for 7.5 second wave period, big amplitude -------

t_min_big_ampl = 5.5;
t_max_big_ampl = 32;

% -------- solid gains for 7.5 second wave period, small amplitude --------

t_min_small_ampl = 5.5;
t_max_small_ampl = 60;

% --------- solid gains for 2 second wave period, small amplitude ---------
% 
% t_min = 1.95;
% t_max = 3.05;

% Generalized frequency bounds, these allow for slower convergence but
% better long term tracking

t_min_accel_x = t_min_small_ampl; % X axis will have small amplitude
t_min_accel_y = t_min_small_ampl; % Y axis will have small amplitude
t_min_accel_z = t_min_big_ampl; % Z axis will have big amplitude

t_max_accel_x = t_max_small_ampl;
t_max_accel_y = t_max_small_ampl;
t_max_accel_z = t_max_big_ampl;

% Accel frequencies
a.beta_min_x_accel = 1/t_max_accel_x;
a.beta_min_y_accel = 1/t_max_accel_y;
a.beta_min_z_accel = 1/t_max_accel_z;

a.beta_max_x_accel = 1/t_min_accel_x;
a.beta_max_y_accel = 1/t_min_accel_y;
a.beta_max_z_accel = 1/t_min_accel_z;

% ---------------------- Assign gyro and accel gains ----------------------

a.kw = [scale_t*a.beta_max_x_accel; scale_t*a.beta_max_y_accel; scale_t*a.beta_max_z_accel; scale_t*a.beta_max_x_gyro; scale_t*a.beta_max_y_gyro; scale_t*a.beta_max_z_gyro];
a.kt = [scale_w*a.beta_min_x_accel; scale_w*a.beta_min_y_accel; scale_w*a.beta_min_z_accel; scale_w*a.beta_min_x_gyro; scale_w*a.beta_min_y_gyro; scale_w*a.beta_min_z_gyro];

%% 3D motion equations

% Saved gains for these equations:

% T_x = 2*pi/(beta/2); % wave frequency [rad/s] (15 sec)
% T_y = 2*pi/(beta/4); % wave frequency [rad/s] (30 sec)
% T_z = 2*pi/(beta); % wave frequency [rad/s] (7.5 sec)

% t_min_x = 13;
% t_max_x = 18;
% 
% t_min_y = 29;
% t_max_y = 33;
% 
% t_min_z = 8;
% t_max_z = 9;

% ------------------- Example Equations & Derived Trig --------------------

% a.real_pos_xI = @(t) 0.4/(beta^2)*sin(beta*t/2);
% a.real_pos_yI = @(t) 1.6/(beta^2)*sin(beta*t/4);
% a.real_pos_zI = @(t) alpha*sin(beta*t) + hdeck;
% 
% a.real_vel_xI = @(t) 0.2/beta*cos(beta/2*t);
% a.real_vel_yI = @(t) 0.4/beta*cos(beta/4*t);
% a.real_vel_zI = @(t) beta*alpha*cos(beta*t);
% 
% a.real_accel_xI = @(t) -0.1*sin(beta/2*t); % [m*s^-2]
% a.real_accel_yI = @(t) -0.1*sin(beta/4*t); % [m*s^-2]
% a.real_accel_zI = @(t) -beta^2*alpha*sin(beta*t) - 9.81; % [m*s^-2]
% 
% a.theta = @(t) -atan(beta*alpha*cos(beta*t)); % [rad]
% a.phi = @(t) atan(0.2/beta*cos(beta/2*t)); % [rad]
% a.psi = @(t) atan(0.4/beta*cos(beta/4*t)); % [rad]
% 
% % Angular rate needs to be taken manually since the derivative equation is +/-
% theta_vals = a.theta(t);
% phi_vals = a.phi(t);
% psi_vals = a.psi(t);
% 
% theta_dot_vals = zeros(1, length(t));
% phi_dot_vals = zeros(1, length(t));
% psi_dot_vals = zeros(1, length(t));
% 
% theta_dot_vals(2:end) = diff(theta_vals)/dt;
% phi_dot_vals(2:end) = diff(phi_vals)/dt;
% psi_dot_vals(2:end) = diff(psi_vals)/dt;
% 
% a.theta_dot = @(t) theta_dot_vals(floor(t./dt)+1);
% a.phi_dot = @(t) phi_dot_vals(floor(t./dt)+1);
% a.psi_dot = @(t) psi_dot_vals(floor(t./dt)+1);

% -------------------------------------------------------------------------



% ------------------- Dynamics Equations used in 3 DOF --------------------

a.real_pos_xI = @(t) 0.5/amp_small*cos((2*pi/period_big)*t);
a.real_pos_yI = @(t) 0.5/amp_small*cos((2*pi/period_big)*t);
a.real_pos_zI = @(t) 0.5/amp_big*cos((2*pi/period_big)*t) + hdeck;

a.real_vel_xI = @(t) -0.5/amp_small*(2*pi/period_big)*sin((2*pi/period_big)*t);
a.real_vel_yI = @(t) -0.5/amp_small*(2*pi/period_big)*sin((2*pi/period_big)*t);
a.real_vel_zI = @(t) -0.5/amp_big*(2*pi/period_big)*sin((2*pi/period_big)*t);

a.real_accel_xI = @(t, y) -0.5/amp_small*((2*pi/period_big)^2)*cos((2*pi/period_big)*t); % [m*s^-2]
a.real_accel_yI = @(t, y) -0.5/amp_small*((2*pi/period_big)^2)*cos((2*pi/period_big)*t); % [m*s^-2]
a.real_accel_zI = @(t, y) -0.5/amp_big*((2*pi/period_big)^2)*cos((2*pi/period_big)*t) - a.g; % [m*s^-2]

a.theta = @(t) 10/amp_small*(pi/180)*sin((2*pi/(period_big))*t); % [rad]
a.phi = @(t) 10/amp_big*(pi/180)*sin((2*pi/(period_big))*t); % [rad]
a.psi = @(t) 10/amp_small*(pi/180)*sin((2*pi/(period_big))*t); % [rad]

a.theta_dot = @(t, y) 20/amp_small*(pi/180)*(pi/period_big)*cos((2*pi/period_big)*t);
a.phi_dot = @(t, y) 20/amp_big*(pi/180)*(pi/period_big)*cos((2*pi/period_big)*t);
a.psi_dot = @(t, y) 20/amp_small*(pi/180)*(pi/period_big)*cos((2*pi/period_big)*t);

% -------------------------------------------------------------------------

%% Standardize Naming (obsolete)

% Inertial Position, Velocity, and Acceleration of Deck

% % Sine Wave
% a.d = @(t) alpha*sin(beta*t) + hdeck;      % [m]
% a.d_dot = @(t) beta*alpha*cos(beta*t);    % [m/s]
% a.d_ddot = @(t) -beta^2*alpha*sin(beta*t); % [m*s^-2]
% 
% a.theta = @(t) atan(beta*alpha*cos(beta*t)); % [rad]
% a.theta_dot = @(t) -(alpha*beta^2*sin(beta*t))/(alpha^2*beta^2*(cos(beta*t)^2)+1); % [rad]

%% Rotation Matrix Helpers

s = @(x) sin(x);
c = @(x) cos(x);
