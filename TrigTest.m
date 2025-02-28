close all; clear; clc;

rng(1,"twister");

% Use same random seed

set(groot,'DefaultLineLineWidth',1)

simulationParameters;

close all;

% Redefining acceleration/gyro curves for clarity

a.real_pos = @(t) alpha*sin(beta*t) + hdeck;
a.real_vel = @(t) beta*alpha*cos(beta*t);
a.real_accel = @(t) -beta^2*alpha*sin(beta*t) - a.g; % [m*s^-2]
a.real_ang = @(t) atan(beta*alpha*cos(beta*t)); % [rad]
a.real_ang_rate = @(t, y) (-(alpha*beta^2*sin(beta*t))./(alpha^2*beta^2*(cos(beta*t).^2)+1)); % [rad/s]

% Simulation time
startTime = 0;
finishTime = 3*60;
% finishTime = 30;
tspan = [startTime finishTime]; % [s]

% Set timestep
freq = 160; % Hz
dt = 1 / freq; % Timestep (s)

% Get time range
t = (tspan(1):dt:tspan(2))';
t_count = length(t);
indeces = @(t) floor(t/dt)+1;

%% Plot prelim environment disturbances

dynamics = @(t, y) a.noisyAccelCurve(t, a.real_accel);
[t, vel]= rk4_solver(dynamics, tspan, a.real_vel(tspan(1)), dt);

function s_dot = DriftCorrection1D(time_i, a, state)

    % Current states
    vel = state(1);
    theta = state(2);
    p = state(3);
    theta_err_accum = state(4);

    curr_state = [vel; theta];
    state_err_accum = [p; theta_err_accum];

    specs = a.specs;

    kw = a.kw; % 
    kt = a.kt; % 

    state_0 = [0; 0];
    state_dot_0 = [0; 0];

    accel_i = a.real_accel(time_i);
    ang_rate_i = a.real_ang_rate(time_i);
    
    accel_dot_m_v = a.measured_accel_vert(time_i, a.o_d_n_a_c_v, a.biasStabDistAccel, a.biasTempDistAccel, a.accel_drift_vert, a.noiseDistAccel, accel_i, a.theta_err);
    gyro_dot_m = a.measured_gyro(time_i, a.o_d_n_g_c, a.biasStabDistGyro, a.biasTempDistGyro, a.gyro_drift, a.noiseDistGyro, ang_rate_i);

    state_dot_m = [accel_dot_m_v; gyro_dot_m];

    state_control = curr_state - state_0;

    state_dot_comp = kt * state_err_accum + kw * state_control;

    state_dot = state_dot_m + state_dot_0 - state_dot_comp;

    s_dot = zeros(length(state),1);

    s_dot(1:2) = state_dot;
    s_dot(3:4) = state_control;

end