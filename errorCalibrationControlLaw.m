close all; clear; clc;

rng(1,"twister");

% Use same random seed

set(groot,'DefaultLineLineWidth',1)

simulationParameters;

a.real_pos = @(t) alpha*sin(beta*t) + hdeck;
a.real_vel = @(t) beta*alpha*cos(beta*t);
a.real_accel = @(t) -beta^2*alpha*sin(beta*t); % [m*s^-2]
a.real_ang = @(t) atan(beta*alpha*cos(beta*t)); % [deg]
a.real_ang_rate = @(t) 180/pi*(-(alpha*beta^2*sin(beta*t))./(alpha^2*beta^2*(cos(beta*t).^2)+1)); % [deg/s]

% Environment and Sensor Specs

specs.alpha = 0.4; % wave amplitdue [m]
specs.hdeck = 1;   % inertial reference deck hight [m] (arbitrary)

% Wave frequency
specs.Tmax = 7.5;    % Maximum period [s]
specs.T = Tmax / 1;  % Period of deck disturbance [s]
specs.beta = 2*pi/T; % wave frequency [rad/s]

imx_5_specs

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

%% Defining Measured Signal

% Simulation time
startTime = 0;
finishTime = 12;
tspan = [startTime finishTime]; % [s]
dt = 0.01;  % [s]
t = (tspan(1):dt:tspan(2))';

% Measured Signal Constants

n_a = 0.5*VRW*t.^(-0.5);
n_g = 0.5*ARW*t.^(-0.5);

% n_a_std = VRW*sqrt(t(2));
% n_g_std = ARW*sqrt(t(2));

accel_noise_std = specs.accel_noiseDensity * sqrt(specs.accel_bandwidth); % Noise standard deviation (m/s^2)
gyro_noise_std = specs.gyro_noiseDensity * sqrt(specs.gyro_bandwidth); % Noise standard deviation (dps)

accel_turn_on_bias_offset = normrnd(0, accel_noise_std);
gyro_turn_on_bias_offset = normrnd(0, gyro_noise_std);

gyro_tau = 15; % [s]
accel_tau = 2; % [s]


% driftPeriod = 5 * 60;  % drift changes every 5 minutes
driftPeriod = 1;  % drift changes every 5 seconds
driftAlterations = floor(tspan(2)/driftPeriod)+1; % drift changes every 5 minutes

b_a_drift_vals = b_a.*rand(driftAlterations,1);
b_g_drift_vals = b_g.*rand(driftAlterations,1);

% b_a_drift_vals = b_a.*ones(driftAlterations,1);
% b_g_drift_vals = b_g.*ones(driftAlterations,1);

bias_indeces = floor(t./(length(t)/driftAlterations)*100)+1;

a.biasStabDistAccel = b_a_drift_vals(bias_indeces);
a.biasStabDistGyro = b_g_drift_vals(bias_indeces);

% Measured Signal Function Handles

a.noiseDistAccel = @(t) accel_noise_std*randn(length(t),1);
a.noiseDistGyro = @(t) gyro_noise_std*randn(length(t),1);

n_a = @(t) VRW*t.^(0.5);
n_g = @(t) ARW*t.^(0.5);

a.biasTempDistAccel = @(t) accel_temp_bias*randn(length(t),1);
a.biasTempDistGyro = @(t) gyro_temp_bias*randn(length(t),1);

a.theta_err = @(t) a.biasStabDistAccel.*t + ARW.*sqrt(t);

a.accel_drift_vert = @(t, real_accel, theta_err) (1 + k)*real_accel(t) + a.biasStabDistAccel + g*(1-cos(theta_err(t)));
a.accel_drift_horz = @(t, real_accel, theta_err) (1 + k)*real_accel(t) + a.biasStabDistAccel + g*sin(theta_err(t));

a.gyro_drift = @(t, real_ang_rate) (1 + k)*real_ang_rate(t) + a.biasStabDistGyro;

a.o_d_n_a_c_v = @(t, biasTempDistAccel, accel_drift_vert, noiseDistAccel, real_accel, theta_err) accel_drift_vert(t, real_accel, theta_err) + biasTempDistAccel(t) + noiseDistAccel(t) + accel_turn_on_bias_offset;
a.o_d_n_a_c_h = @(t, biasTempDistAccel, accel_drift_horz, noiseDistAccel, real_accel, theta_err) accel_drift_horz(t, real_accel, theta_err) + biasTempDistAccel(t) + noiseDistAccel(t) + accel_turn_on_bias_offset;
a.o_d_n_g_c = @(t, biasTempDistGyro, gyro_drift, noiseDistGyro, real_ang_rate) gyro_drift(t, real_ang_rate) + biasTempDistGyro(t) + noiseDistGyro(t) + gyro_turn_on_bias_offset;

a.measured_accel_vert = @(t, o_d_n_a_c_v, biasTempDistAccel, accel_drift_vert, noiseDistAccel, real_accel, theta_err) accel_resolution*floor(o_d_n_a_c_v(t, biasTempDistAccel, accel_drift_vert, noiseDistAccel, real_accel, theta_err)/accel_resolution);
a.measured_accel_horz = @(t, o_d_n_a_c_h, biasTempDistAccel, accel_drift_horz, noiseDistAccel, real_accel, theta_err) accel_resolution*floor(o_d_n_a_c_h(t, biasTempDistAccel, accel_drift_horz, noiseDistAccel, real_accel, theta_err)/accel_resolution);
a.measured_gyro = @(t, o_d_n_g_c, biasTempDistGyro, gyro_drift, noiseDistGyro, real_ang_rate) gyro_resolution*floor(o_d_n_g_c(t, biasTempDistGyro, gyro_drift, noiseDistGyro, real_ang_rate)/gyro_resolution);

measured_accel_vert = a.measured_accel_vert(t, a.o_d_n_a_c_v, a.biasTempDistAccel, a.accel_drift_vert, a.noiseDistAccel, a.real_accel, a.theta_err);
measured_accel_horz = a.measured_accel_horz(t, a.o_d_n_a_c_h, a.biasTempDistAccel, a.accel_drift_horz, a.noiseDistAccel, a.real_accel, a.theta_err);
measured_gyro = a.measured_gyro(t, a.o_d_n_g_c, a.biasTempDistGyro, a.gyro_drift, a.noiseDistGyro, a.real_ang_rate);

% Initial States
p0 =  a.d(tspan(1))+a.pr_d;   % Platform position [m]
p_dot0 = a.d_dot(tspan(1));   % Platform velocity [m/s]
pr_err_accum0 = 0;            % Integral of relative position error [m*s]
pm0 = p0;                     % Platform inetegrated position [m]
pm_dot = p_dot0;              % Platform integrated velocity [m/s]
pm_ddot = a.d_ddot(tspan(1)); % Platform measured acceleration [m*s^-2]
p_theta0 = a.real_ang(tspan(1)); % Platform inertial angle [deg]

s0 = [p0; p_dot0; pr_err_accum0; pm0; pm_dot; pm_ddot];

% Plotting readl disturbance vs measured disturbance

figure
a.fi1 = axes;
xlabel('Time (s)')
ylabel('Drift (m/s^-2)')
title('Vertical Accelerometer Signal vs Time')
%% Running simulations

% --------------------- Running no error simulation -----------------------

% tstep = 1/60;
% startTime = 0;
% finishTime = 10;
% t = startTime:tstep:finishTime;
% s = zeros(length(t), length(s0));
% 
% s(1,:) = noError(t(1),s0,a);
% sStep = s(1,:).*tstep+s0;
% s(2,:) = noError(t(2),sStep,a);
% 
% for i = 2:length(t)
%     time = t(i);
%     sStep = s(i,:).*tstep+s(i-1,:);
%     s(i+1,:) = noError(t(i),sStep,a);
% end

tolerance = 6e-6;
% Running Simulation
op = odeset('RelTol',tolerance,'AbsTol',tolerance); % Tolerance options
[t_reg, s_reg] = ode23s(@(t,s)noError(t,s,a),tspan,s0,op);

% % Fixed-step integrator
% dt = 0.01; % s, 10ms
% type = 2;
% s_i_fixed = [p0              0;
%              p_dot0          0;
%              pr_err_accum0   0;
%              pm0             0;
%              pm_dot          0;
%              pm_ddot         0];
% 
% s_fixed = fixedIntegration(s_i_fixed, dt, tspan, type, a);

% --------------------- Running with error simulation ---------------------

s0 = [p0; p_dot0; pr_err_accum0; pm0; pm_dot; pm_ddot; p_theta0];

% Running Simulation
op = odeset('RelTol',tolerance,'AbsTol',tolerance); % Tolerance options
[t, s] = ode23s(@(t,s)rigidArmControl(t,s,a),tspan,s0,op);

% plot(a.fi1, t,a.real_accel(t), color='black')

% Feeding states back through EOM to calculating inertial acceleration of
% the platfor
p_ddot = zeros(size(t));
p_theta = zeros(size(t));
for i = 1:length(t)
    s_dot = rigidArmControl(t(i),s(i,:),a);
    p_ddot(i) = s_dot(6);
    p_theta(i) = s_dot(7);
end

%% Plotting

% With error ------------------------------------------------------------

% Plotting Position vs Time and Acceleration vs Time
figure;
sgtitle('Inertial Stability Performance')

% Position
subplot(1,3,1);
plot(t,s(:,1))
hold on
plot(t_reg, s_reg(:,1))
hold on
plot(t,a.d(t))
title('Inertial Position vs Time')
xlabel('Time (s)')
ylabel('Position (m)')
legend('Platform with Error', 'Platform without Error', 'Deck','Location','southeast')

% Position
subplot(1,3,2);
plot(t,s(:,2))
hold on
plot(t,a.d_dot(t))
title('Inertial Velocity vs Time')
xlabel('Time (s)')
ylabel('Velocity (m/s)')
legend('Platform', 'Deck','Location','southeast')

% Acceleration
subplot(1,3,3);
plot(t,p_ddot)
hold on
plot(t,a.d_ddot(t))
hold on
yline(p_ddot_max,'--')
yline(-p_ddot_max,'--')
title('Inertial Acceleration vs Time')
xlabel('Time (s)')
ylabel('Acceleration (m/s^2)')
legend('Platform', 'Deck','Location','southeast')

figure;
subplot(1,2,1)
plot(t, a.real_ang_rate(t))
xlabel('Time (s)')
ylabel('Angular Rate (deg/s)')
title('Platform Angle Before Control')
subplot(1,2,2)
plot(t, p_theta)
xlabel('Time (s)')
ylabel('Angular Rate (deg/s)')
title('Platform Angle After Control')

figure;
subplot(1,2,1)
plot(t, a.real_accel(t))
xlabel('Time (s)')
ylabel('Acceleration (m/s^2)')
title('Platform Acceleration Before Control')
subplot(1,2,2)
plot(t, p_ddot)
xlabel('Time (s)')
ylabel('Acceleration (m/s^2)')
title('Platform Acceleration After Control')

% Plotting relative position
figure;
plot(t,s(:,1)-a.d(t))
hold on
% Plotting inertial control region
x = [tspan, flip(tspan)];
yf = [a.pr_d-a.r_g, a.pr_d-a.r_g, a.pr_d+a.r_g, a.pr_d+a.r_g];
fill(x,yf,'y','FaceAlpha',0.2,'EdgeColor','none')
% Plotting Relative position control region
x = [tspan, flip(tspan)];
yta = [a.pr_d+a.r_k, a.pr_d+a.r_k, 1, 1];
ytb = [0, 0, a.pr_d-a.r_k, a.pr_d-a.r_k];
fill(x,yta,'b','FaceAlpha',0.2,'EdgeColor','none')
fill(x,ytb,'b','FaceAlpha',0.2,'EdgeColor','none')
yline(a.pr_d,'--','Label','$p_{rd}$','Interpreter','latex','FontSize',15)
title('Relative Position vs Time')
xlabel('Time (s)')
ylabel('Position (m)')
legend('','Full Inertial','Relative Position','')

%% Getting error

% Checking that both simulations finished to desired time
if t_reg(end) ~= finishTime
    fprintf('Simulation (without error introduced) failed to finish.')
    return
elseif t(end) ~= finishTime
    fprintf('Simulation (with error introduced) failed to finish.')
    return
end

% Plotting time
sz = 2;
figure
scatter(1:length(t), t, sz, 'filled', displayName="With Error")
hold on
scatter(1:length(t_reg), t_reg, sz, 'filled', displayName="Without Error")
ylabel('Time Values')
title('Timesteps used in Integration')
legend

% Plotting both positions, with and without error
figure;
plot(t, s(:,1))
hold on
plot(t_reg, s_reg(:,1))
title('Inertial Position vs Time')
xlabel('Time (s)')
ylabel('Position (m)')
legend('Platform Position with Error', 'Platform Position without Error')

% Precision (number of decimals) of interpolation
err_round = 3;

% Getting the more precise time matrix
if size(t_reg) > size(t)
    t_precise = t_reg;
    p_precise = s_reg(:,1);
    t_compare = t;
    p_compare = s(:,1);
else
    t_precise = t;
    p_precise = s(:,1);
    t_compare = t_reg;
    p_compare = s_reg(:,1);
end

finishIndex = length(t_precise);
pos_err = zeros(1, finishIndex);

for i=1:finishIndex
    
    time_precise = t_precise(i);
    pos_precise = p_precise(i);

    t_ref = round(time_precise, err_round);
    time_interp_idx = findNearest(t_ref, t_compare);

    time_interp = t_compare(time_interp_idx);
    pos_interp = p_compare(time_interp_idx);
    
    pos_err(i) = abs(pos_precise - pos_interp);
end

%% Plotting Error

sz = 2;
% figure
% subplot(3,1,1)
% scatter(t_ref, pos_err_worse*100, sz, 'filled', displayName="Positional Error")
% hold on
% plot(t,a.d(t)/50, displayName="Deck Disturbance")
% title('Worst Relative Position Error vs Time')
% xlabel('Time (s)')
% ylabel('Error (cm)')
% legend
% 
% subplot(3,1,2)
% scatter(t_ref, pos_err_worse*100, sz, 'filled', displayName="Positional Error")
% hold on
% plot(t,a.d(t)/50, displayName="Deck Disturbance")
% title('Average Relative Position Error vs Time')
% xlabel('Time (s)')
% ylabel('Error (cm)')
% legend
% 
% subplot(3,1,3)
% scatter(t_ref, pos_err_worse*100, sz, 'filled', displayName="Positional Error")
% hold on
% plot(t,a.d(t)/50, displayName="Deck Disturbance")
% title('Best Relative Position Error vs Time')
% xlabel('Time (s)')
% ylabel('Error (cm)')
% legend

figure
scatter(t_precise, pos_err*100, sz, 'filled', displayName="Positional Error")
hold on
plot(t,a.d(t)/200, displayName="Deck Disturbance")
title('Worst Case Relative Position Error vs Time')
xlabel('Time (s)')
ylabel('Error (cm)')
legend

growth = diff(pos_err);

figure
scatter(t_precise(2:end), growth*100, sz, 'filled', displayName="Positional Error")
hold on
plot(t,a.d(t)/200, displayName="Deck Disturbance")
title('Error Growth over Time')
xlabel('Time (s)')
ylabel('Error (cm)')
legend

%% Functions

function min_idx = findNearest(t_ref, t_compare)
    min_diff = 1;
    min_idx = 0;
    for idx=1:length(t_compare)
        curr_t = t_compare(idx);
        if abs(curr_t-t_ref) < min_diff
            min_diff = abs(curr_t-t_ref);
            min_idx = idx;
        end
    end
end

function s_dot = rigidArmControl(t, s, a)
% rigidArmControl is the EOM for the 1 DOF model of the inertially
% stabilized platform. It uses inertial acceleration control when the
% platform is close to the center of the operation region, and uses PID
% control on the relative position as the platform goes closer to the
% operational bounderies
%
% Inputs:   t    = current time
%           s    = vector of states
%                = [p; p_dot; pr_err_accum] where p is the inertial position 
%                  of the platform, pdot is the inertial velocity of the 
%                  platform, and pr_err_accum is the integral of the error 
%                  in the relative position of the platform
%           a    = structure containing environmental constants and gain
%                  values
% Outputs:  sdot = time derivative of input state vector
%                = [p_dot; p_ddot; pr_err] where pdot is the inertial 
%                  velocity of the platform, p_ddot is the inertial 
%                  acceleration of the platform,and pr_err is the error in 
%                  the relative position of the platform

% Current states
p = s(1);
p_dot = s(2);
pr_err_accum = s(3);
pm = s(4);
pm_dot = s(5);
pm_ddot = s(6);
p_theta = s(7);

specs = a.specs;

p_thetadot = 1/(pm_dot^2 + 1);



% Inserting measured accel manually ---------------------------------------
sz = 3;
if t > specs.accel_resolution

    measured_a_x = a.measured_accel_horz(t, a.o_d_n_a_c_h, a.biasTempDistAccel, a.accel_drift_horz, a.noiseDistAccel, a.real_accel, a.theta_err);
    measured_a_y = measured_a_x;
    measured_a_z = a.measured_accel_vert(t, a.o_d_n_a_c_v, a.biasTempDistAccel, a.accel_drift_vert, a.noiseDistAccel, a.real_accel, a.theta_err);

    measured_g_x = a.measured_gyro(t, a.o_d_n_g_c, a.biasTempDistGyro, a.gyro_drift, a.noiseDistGyro, a.real_ang_rate);
    measured_g_y = measured_g_x;
    measured_g_z = measured_g_x;

    measuredState = [measured_a_x measured_a_y measured_a_z measured_g_x measured_g_y measured_g_z];
    corrected_a = compensateError(measuredState, specs, t);

    pm_ddot_error = pm_ddot - corrected_a(3);
    p_thetadot_error = p_thetadot - corrected_a(4);

    pm_ddot = pm_ddot + (pm_ddot_error/1);
    p_thetadot = p_thetadot + (p_thetadot_error/1);

    % pm_ddot = corrected_a(3);
    % p_thetadot = corrected_a(4);

    % scatter(a.fi1, t, a_error, "red")
    % hold on
    % 
    % scatter(a.fi1, t, pm_ddot, "blue")
    % hold on
    % 
    % scatter(a.fi1, t, corrected_a(3), "green")
    % hold on

end


% scatter(a.fi1, t, corrected_a(1), sz, "filled", color='blue')
% hold on
% scatter(a.fi1, t, corrected_a(3), sz, color='blue')
% hold on

% 
% error = abs(corrected_a(3) - p_ddot);
% scatter(a.fi1, t, error, sz, "red")
% hold on
% scatter(a.fi3, t, corrected_a(4), sz, "filled", color='blue')
% hold on

% ----------------------------------------------------------------------------------------------

% Error in relative position (distance to center of operation region)
pr = pm-a.d(t);
pr_err = pr-a.pr_d;

% For testing gains without mixing proportions
% ka = a.ka; % Acceleration [kg]
% kv = a.kv; % Velocity     [kg/s]
% ks = a.ks; % Position     [kg*s^-2]
% kp = a.kp; % Proportional [kg*s^-2]
% kd = a.kd; % Derivative   [kg/s]
% ki = a.ki; % Integral     [kg*s^-3]

% Control gain proportions

% I = a.I(pr); % Proportion of inertial stability control to apply
I = 0.5;

ka = a.ka*I; % Acceleration [kg]
kv = a.kv*I; % Velocity     [kg/s]
ks = a.ks*I; % Position     [kg*s^-2]

% k = a.K(pr);     % Proportion of relative position control to apply
k = 0.5;

k_h = a.K_h(pr);
kp = a.kp*k_h;     % Proportional [kg*s^-2]
kd = a.kd*k_h;     % Derivative   [kg/s]
ki = a.ki*k_h;     % Integral     [kg*s^-3]

% Derivative of states
s_dot = zeros(7,1);

% Derivative of position
s_dot(1) = p_dot;  % inertial velocity
s_dot(4) = pm_dot; % measured inertial velocity

% Control Law

% Inertial stability control force
c_i = a.initial_scale(t); % Initial scale of gains
f_i = -(ka*pm_ddot + kv*pm_dot + ks*pm)*c_i;
% Relative position control force
f_pr = -(kp*pr_err + ki*pr_err_accum + kd*(pm_dot-a.d_dot(t)));

% Platform EOM
p_ddot = (f_i+f_pr) / a.m;

% scatter(a.fi1, t, p_ddot, sz, "blue")
% hold on

% Derivative of velocity
s_dot(2) = p_ddot; % Inertial acceleration
s_dot(5) = pm_ddot; % Measured inertial acceleration

% Derivative of measured inertial acceleration
s_dot(6) = a.omega*(p_ddot - pm_ddot);

% Error in relative position
s_dot(3) = pr_err;

% Derivative of angle
s_dot(7) = p_thetadot; % Angular Rate

err_v=abs(p_dot-pm_dot);
err_a=abs(p_ddot-pm_ddot);

t

end

function state = compensateError(measuredState, specs, time)
    k = specs.k;
    dk = specs.dk;
    b_a = specs.b_a;
    b_g = specs.b_g;
    ARW = specs.ARW;
    VRW = specs.VRW;
    accel_noiseDensity = specs.accel_noiseDensity;
    gyro_noiseDensity = specs.gyro_noiseDensity;
    accel_resolution = specs.accel_resolution;
    gyro_resolution = specs.gyro_resolution;
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

    theta_err = b_g/2*time + ARW*sqrt(time);

    a_adjusted_x_y = measured_accel(1) - ACCEL_BIAS(1) - g*sin(theta_err);
    a_adjusted_z = measured_accel(3) - ACCEL_BIAS(3) - g*(1-cos(theta_err));
    a_adjusted = [a_adjusted_x_y; a_adjusted_x_y; a_adjusted_z];

    % Adding calibration error of up to 2 bits of accuracy
    a_adjusted = a_adjusted - 2*accel_resolution;
    
    A_FIX = inv([1+S_x+dS_x  M_xy       M_xz
                M_yx         1+S_y+dS_y M_yz
                M_zx         M_zy       1+S_z+dS_z]);

    corrected_a = A_FIX * a_adjusted;


    G_DEP_BIAS = [B_gx  0     0
                  0     B_gy  0
                  0     0     B_gz];

    g_adjusted = measured_gyro - GYRO_BIAS - G_DEP_BIAS*corrected_a;
    
    % Adding calibration error of up to 2 bits of accuracy
    g_adjusted = g_adjusted - 2*gyro_resolution;

    corrected_g = A_FIX * g_adjusted;

    state(1:3) = corrected_a';
    state(4:6) = corrected_g';
    
end


function s_dot = noError(t, s, a)
% rigidArmControl is the EOM for the 1 DOF model of the inertially
% stabilized platform. It uses inertial acceleration control when the
% platform is close to the center of the operation region, and uses PID
% control on the relative position as the platform goes closer to the
% operational bounderies
%
% Inputs:   t    = current time
%           s    = vector of states
%                = [p; p_dot; pr_err_accum] where p is the inertial position 
%                  of the platform, pdot is the inertial velocity of the 
%                  platform, and pr_err_accum is the integral of the error 
%                  in the relative position of the platform
%           a    = structure containing environmental constants and gain
%                  values
% Outputs:  sdot = time derivative of input state vector
%                = [p_dot; p_ddot; pr_err] where pdot is the inertial 
%                  velocity of the platform, p_ddot is the inertial 
%                  acceleration of the platform,and pr_err is the error in 
%                  the relative position of the platform

% Current states
p = s(1);
p_dot = s(2);
pr_err_accum = s(3);
pm = s(4);
pm_dot = s(5);
pm_ddot = s(6);

% Error in relative position (distance to center of operation region)
pr = pm-a.d(t);
pr_err = pr-a.pr_d;

% For testing gains without mixing proportions
% ka = a.ka; % Acceleration [kg]
% kv = a.kv; % Velocity     [kg/s]
% ks = a.ks; % Position     [kg*s^-2]
% kp = a.kp; % Proportional [kg*s^-2]
% kd = a.kd; % Derivative   [kg/s]
% ki = a.ki; % Integral     [kg*s^-3]

% Control gain proportions

I = a.I(pr); % Proportion of inertial stability control to apply
ka = a.ka*I; % Acceleration [kg]
kv = a.kv*I; % Velocity     [kg/s]
ks = a.ks*I; % Position     [kg*s^-2]

k = a.K(pr);     % Proportion of relative position control to apply
k_h = a.K_h(pr);
kp = a.kp*k_h;     % Proportional [kg*s^-2]
kd = a.kd*k_h;     % Derivative   [kg/s]
ki = a.ki*k_h;     % Integral     [kg*s^-3]

% Derivative of states
s_dot = zeros(6,1);

% Derivative of position
s_dot(1) = p_dot;  % inertial velocity
s_dot(4) = pm_dot; % measured inertial velocity

% Control Law

% Inertial stability control force
c_i = a.initial_scale(t); % Initial scale of gains
f_i = -(ka*pm_ddot + kv*pm_dot + ks*pm)*c_i;
% Relative position control force
f_pr = -(kp*pr_err + ki*pr_err_accum + kd*(pm_dot-a.d_dot(t)));

% Platform EOM
p_ddot = (f_i+f_pr) / a.m;

% Derivative of velocity
s_dot(2) = p_ddot; % Inertial acceleration
s_dot(5) = pm_ddot; % Measured inertial acceleration

% Derivative of measured inertial acceleration
s_dot(6) = a.omega*(p_ddot - pm_ddot) ;

% Error in relative position
s_dot(3) = pr_err;

% s_dot

err_v=abs(p_dot-pm_dot);
err_a=abs(p_ddot-pm_ddot);

end

function states = KalmanFilter(t, signal, noise_std, q)
    n = length(t);
    
    dim = size(signal, 1);

    err_measure = noise_std;
    err_estimate = err_measure;

    states = zeros(size(signal));

    last_estimate = signal(:, 1);

    for i = 1:n
        mea = signal(:, i);
        
        E = err_measure + err_estimate;

        kalman_gain = err_estimate ./ E;
        K = diag(kalman_gain);

        curr_estimate = last_estimate + K * (mea - last_estimate);

        diff = diag(abs(last_estimate - curr_estimate));
        err_estimate = (diag(ones(dim,1)) - K)*err_estimate + diff*q;

        last_estimate = curr_estimate;
        states(:, i) = curr_estimate;
    end
end

function s_fixed = fixedIntegration(s_i_fixed, dt, tspan, type, a)
    startTime = tspan(1);
    finishTime = tspan(2);

    s_fixed = [];

    % Derivative of states
    s_dot = zeros(size(s_i_fixed));

    figure
    hold on

    for t=startTime:dt:finishTime
        t
    
        % Current states
        p = s_i_fixed(1, 1);
        p_dot = s_i_fixed(2, 1);
        pr_err_accum = s_i_fixed(3, 1);
        pm = s_i_fixed(4, 1);
        pm_dot = s_i_fixed(5, 1);
        pm_ddot = s_i_fixed(6, 1);
        
        % Error in relative position (distance to center of operation region)
        pr = pm-a.d(t);
        pr_err = pr-a.pr_d;
        
        % For testing gains without mixing proportions
        % ka = a.ka; % Acceleration [kg]
        % kv = a.kv; % Velocity     [kg/s]
        % ks = a.ks; % Position     [kg*s^-2]
        % kp = a.kp; % Proportional [kg*s^-2]
        % kd = a.kd; % Derivative   [kg/s]
        % ki = a.ki; % Integral     [kg*s^-3]
        
        % Control gain proportions
        
        I = a.I(pr); % Proportion of inertial stability control to apply
        ka = a.ka*I; % Acceleration [kg]
        kv = a.kv*I; % Velocity     [kg/s]
        ks = a.ks*I; % Position     [kg*s^-2]
        
        k = a.K(pr);     % Proportion of relative position control to apply
        k_h = a.K_h(pr);
        kp = a.kp*k_h;     % Proportional [kg*s^-2]
        kd = a.kd*k_h;     % Derivative   [kg/s]
        ki = a.ki*k_h;     % Integral     [kg*s^-3]
        
        % Derivative of position
        s_dot(1, 1) = p_dot;  % inertial velocity
        s_dot(4, 1) = pm_dot; % measured inertial velocity
        
        % Control Law
        
        % Inertial stability control force
        c_i = a.initial_scale(t); % Initial scale of gains
        f_i = -(ka*pm_ddot + kv*pm_dot + ks*pm)*c_i;
        % Relative position control force
        f_pr = -(kp*pr_err + ki*pr_err_accum + kd*(pm_dot-a.d_dot(t)));
        
        % Platform EOM
        p_ddot = (f_i+f_pr) / a.m;
        
        % Derivative of velocity
        s_dot(2, 1) = p_ddot; % Inertial acceleration
        s_dot(5, 1) = pm_ddot; % Measured inertial acceleration
        
        % Derivative of measured inertial acceleration
        s_dot(6, 1) = a.omega*(p_ddot - pm_ddot);
        
        % Error in relative position
        s_dot(3, 1) = pr_err;
    
        % Integrate from derivative
        s_i_fixed(:, 1) = fdm_integrator(s_i_fixed(:, 2), s_dot, dt, type);
        
        s_i_fixed
        s_dot

        % Update y_dot_dot
        s_dot = insert_value(s_dot, s_dot(:,1));

        s_dot

        % U = insert_value(U, fdm_integrator(U, U_dot, dt, 2)); 
        
        s_fixed = [s_fixed s_i_fixed(:, 1)];
        
        scatter(t, s_fixed(1))
        hold on

    end

end
