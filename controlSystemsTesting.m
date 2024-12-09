close all; clear; clc;

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
% specs.accel_noiseDensity = 0;
% specs.gyro_noiseDensity = 0;
% specs.accel_resolution = 0.000001;
% specs.gyro_resolution = 0.000001;
% specs.accel_noiseDensity = specs.accel_noiseDensity / 100;
% specs.gyro_noiseDensity = specs.gyro_noiseDensity / 100;
a.specs = specs;

% Measured signal

a.n_a = @(t) 0.5*specs.VRW*t.^(-0.5);
a.n_g = @(t) 0.5*specs.ARW*t.^(-0.5);
a.theta_err = @(t) specs.b_g*t + specs.ARW*sqrt(t);

a.accel_noise_std = specs.accel_noiseDensity * sqrt(specs.accel_bandwidth); % Noise standard deviation (microg)
a.gyro_noise_std = specs.gyro_noiseDensity * sqrt(specs.gyro_bandwidth); % Noise standard deviation (microg)

% a.accel_quantized = @(t, real_accel) specs.accel_resolution*floor(a.real_accel(t)/specs.accel_resolution);
% a.gyro_quantized = @(t, real_ang_rate) specs.gyro_resolution*floor(a.real_ang_rate(t)/specs.gyro_resolution);
% 
% a.quant_noise_accel = @(t, accel_quantized) a.accel_quantized(t, a.real_accel) + normrnd(0,a.accel_noise_std);
% a.quant_noise_gyro = @(t, gyro_quantized) a.gyro_quantized(t, a.real_ang_rate) + normrnd(0,a.gyro_noise_std);

% a.drift_error_accel_vert = @(t, n_a, quant_noise_accel, theta_err) (1 + specs.k)*a.quant_noise_accel(t, a.accel_quantized) + specs.b_a + a.n_a(t) + a.g*(1-cos(a.theta_err(t)));
% a.drift_error_accel_horz = @(t, n_a, quant_noise_accel, theta_err) (1 + specs.k)*a.quant_noise_accel(t, a.accel_quantized) + specs.b_a + a.n_a(t) + a.g*sin(theta_err(t));
% a.drift_error_gyro = @(t, n_g, quant_noise_gyro, theta_err) (1 + specs.k)*a.quant_noise_gyro(t, a.gyro_quantized) + specs.b_g + a.n_g(t);


a.accel_quantized = @(t, p_ddot) specs.accel_resolution*floor(p_ddot/specs.accel_resolution);
a.gyro_quantized = @(t, p_thetadot) specs.gyro_resolution*floor(p_thetadot/specs.gyro_resolution);

a.quant_noise_accel = @(t, accel_quantized, p_ddot) a.accel_quantized(t, p_ddot) + normrnd(0,a.accel_noise_std);
a.quant_noise_gyro = @(t, gyro_quantized, p_thetadot) a.gyro_quantized(t, p_thetadot) + normrnd(0,a.gyro_noise_std);

a.drift_error_accel_vert = @(t, n_a, quant_noise_accel, theta_err, p_ddot) (1 + specs.k)*a.quant_noise_accel(t, a.accel_quantized, p_ddot) + specs.b_a + a.n_a(t) + a.g*(1-cos(a.theta_err(t)));
a.drift_error_accel_horz = @(t, n_a, quant_noise_accel, theta_err, p_ddot) (1 + specs.k)*a.quant_noise_accel(t, a.accel_quantized, p_ddot) + specs.b_a + a.n_a(t) + a.g*sin(theta_err(t));
a.drift_error_gyro = @(t, n_g, quant_noise_gyro, theta_err, p_thetadot) (1 + specs.k)*a.quant_noise_gyro(t, a.gyro_quantized, p_thetadot) + specs.b_g + a.n_g(t);

% Simulation time
% startTime = T/4;
startTime = 0;
tspan = [startTime 10]; % [s]

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


% Running Simulation
op = odeset('RelTol',1e-8,'AbsTol',1e-8); % Tolerance options
[t_reg, s_reg] = ode23s(@(t,s)noError(t,s,a),tspan,s0,op);

s0 = [p0; p_dot0; pr_err_accum0; pm0; pm_dot; pm_ddot; p_theta0];

% Running Simulation
op = odeset('RelTol',1e-8,'AbsTol',1e-8); % Tolerance options
[t, s] = ode23s(@(t,s)rigidArmControl(t,s,a),tspan,s0,op);

% plot(a.fi1, t,a.real_accel(t), color='black')

% Feeding states back through EOM to calculating inertial acceleration of
% the platfor
p_ddot = zeros(size(t));
p_thetadot = zeros(size(t));
for i = 1:length(t)
    s_dot = rigidArmControl(t(i),s(i,:),a);
    p_ddot(i) = s_dot(2);
    p_thetadot(i) = s_dot(7);
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
title('Platform Angular Rate Before Control')
subplot(1,2,2)
plot(t, p_thetadot)
xlabel('Time (s)')
ylabel('Angular Rate (deg/s)')
title('Platform Angular Rate After Control')

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

sim_t = round(t(end),3);
time_index = find(round(t_reg,3) == sim_t);
time_index = time_index(1);
plat_err = s(:,1);
plat_no_err = s_reg(:,1);
pos_err_worse = zeros(1, time_index);
pos_err_avg = zeros(1, time_index);
pos_err_best = zeros(1, time_index);

figure;
plot(t,s(:,1))
hold on
plot(t_reg(1:time_index), s_reg(1:time_index,1))
title('Relative Position vs Time')
xlabel('Time (s)')
ylabel('Position (m)')
legend('Platform Position with Error', 'Platform Position without Error')

for i=1:time_index
    platform_without_err = round(plat_no_err(i),3);
    integrated_index = find(round(plat_err,3) == platform_without_err);
    platform_with_err = plat_err(integrated_index(1));
    all_err = abs(platform_with_err - plat_no_err(i));
    pos_err_worse(i) = max(all_err);
    pos_err_best(i) = min(all_err);
    pos_err_avg(i) = mean(all_err);
end

%% Plotting Error
sz = 2;
figure
subplot(3,1,1)
scatter(t_reg(1:time_index), pos_err_worse*100, sz, 'filled', displayName="Positional Error")
hold on
plot(t,a.d(t)/50, displayName="Deck Disturbance")
title('Worst Relative Position Error vs Time')
xlabel('Time (s)')
ylabel('Error (cm)')
legend

subplot(3,1,2)
scatter(t_reg(1:time_index), pos_err_worse*100, sz, 'filled', displayName="Positional Error")
hold on
plot(t,a.d(t)/50, displayName="Deck Disturbance")
title('Average Relative Position Error vs Time')
xlabel('Time (s)')
ylabel('Error (cm)')
legend

subplot(3,1,3)
scatter(t_reg(1:time_index), pos_err_worse*100, sz, 'filled', displayName="Positional Error")
hold on
plot(t,a.d(t)/50, displayName="Deck Disturbance")
title('Best Relative Position Error vs Time')
xlabel('Time (s)')
ylabel('Error (cm)')
legend

figure
scatter(t_reg(1:time_index), pos_err_worse*100, sz, 'filled', displayName="Positional Error")
hold on
plot(t,a.d(t)/50, displayName="Deck Disturbance")
title('Worst Case Relative Position Error vs Time')
xlabel('Time (s)')
ylabel('Error (cm)')
legend


%%

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

    measured_a_x = a.drift_error_accel_horz(t, a.n_a, a.quant_noise_accel, a.theta_err, pm_ddot);
    measured_a_y = measured_a_x;
    measured_a_z = a.drift_error_accel_vert(t, a.n_a, a.quant_noise_accel, a.theta_err, pm_ddot);

    measured_g_x = a.drift_error_gyro(t, a.n_g, a.quant_noise_gyro, a.theta_err, p_thetadot);
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
s_dot(6) = a.omega*(p_ddot - pm_ddot);

% Error in relative position
s_dot(3) = pr_err;

err_v=abs(p_dot-pm_dot);
err_a=abs(p_ddot-pm_ddot);

end