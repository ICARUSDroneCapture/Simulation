Aclose all; clear; clc;

sampling_rate = 10; %1 Hz
time_step = 1/sampling_rate;
accel_resolution = 0.4; % mg
accel_resolution = accel_resolution/1000*9.8; % m/s
bias_instability = 19*10^-6; % microg
bias_instability = bias_instability * 9.8;

tspan = 1:time_step:50;

%% Modeling Platform Position Accounting for Error

accel_resolution = 0.122 / 1000 * 9.8; % m/s

figure;

a.error = accel_resolution;
op = odeset('RelTol',1e-12,'AbsTol',1e-12);
[t_err, s_err] = ode45(@(t,s)controlWithError(t,s,a),tspan,s0,op);

% Position
subplot(1,2,1);
plot(t,s(:,1), Color="blue")
hold on
plot(t_err,s_err(:,1))
hold on
plot(t,a.d(t))
title('Position vs Time with Error')
xlabel('Time (s)')
ylabel('Position (m)')
legend('Platform (no error)', 'Platform (considering error)', 'Deck')

% Acceleration
subplot(1,2,2);
% Feeding states back through EOM to calculating inertial acceleration of
% the platfor
p2dot_err = zeros(size(t_err));
for i = 1:length(t_err)
    sdot_err = controlWithError(t_err(i),s_err(i,:),a);
    p2dot_err(i) = sdot_err(2);
end
plot(t_err,p2dot_err)
hold on
plot(t,a.d_ddot(t))
yline(p_ddot_max,'--')
yline(-p_ddot_max,'--')
title('Acceleration vs Time')
xlabel('Time (s)')
ylabel('Acceleration (m/s^2)')
legend('Platform (considering error)', 'Deck')

%% Calculating Max Error as a Function of Resolution

figure(6)

error_dev = 100;

deck_pos_error = zeros(1, error_dev);
deck_acc_error = zeros(1, error_dev);

min_exp_acc = 0.357; % mg
min_exp_acc = 0.357 * 9.8 / 1000; % m/s
bias_instability = 19*10^-6 * 9.8; % m/s

min_res = min_exp_acc / 2;
max_res = min_exp_acc / 10;

range_res = linspace(min_res, max_res, error_dev);

for err = 1:error_dev
    a.error = range_res(err) + bias_instability;
    [t_err, s_err] = ode45(@(t,s)controlWithError(t,s,a),tspan,s0,op);
    p2dot_err = zeros(size(t));
    for i = 1:length(t)
        sdot_err = controlWithError(t(i),s(i,:),a);
        p2dot_err(i) = sdot_err(2);
    end
    deck_pos_error(err) = abs(max(s(:,1)) - max(s_err(:,1)));
    deck_acc_error(err) = abs(max(p2dot_err) - max(p_ddot));
end

subplot(1, 2, 1)
plot(range_res/9.8*1000, deck_pos_error*100)
yline(0.5, '--')
ylim([0 0.75])
title('Position Deviation vs Time')
xlabel('Resolution (mg)')
ylabel('Position Error (cm)')

subplot(1, 2, 2)
plot(range_res/9.8*1000, deck_acc_error)
title('Acceleration Deviation vs Time')
xlabel('Resolution (mg)')
ylabel('Acceleration Error (m/s)')


%% Comparing without error, with error, with more error


figure(7)

error_dev = 100;

deck_pos_error = zeros(1, error_dev);
deck_acc_error = zeros(1, error_dev);

min_exp_acc = 0.357; % mg
min_exp_acc = 0.357 * 9.8 / 1000; % m/s

min_res = min_exp_acc / 2;
max_res = min_exp_acc / 10;

range_res = linspace(min_res, max_res, error_dev);
bias_instability = 19*10^-6; % microg
bias_instability = bias_instability * 9.8;

for err = 1:error_dev
    a.error = 0;
    [t_err, s_err] = ode45(@(t,s)controlWithError(t,s,a),tspan,s0,op);
    p2dot_err = zeros(size(t));
    for i = 1:length(t)
        sdot_err = controlWithError(t(i),s(i,:),a);
        p2dot_err(i) = sdot_err(2);
    end
    deck_pos_error(err) = abs(max(s(:,1)) - max(s_err(:,1)));
    deck_acc_error(err) = abs(max(p2dot_err) - max(p_ddot));
end

plot(range_res/9.8*1000, abs(deck_pos_error)*100, DisplayName="No Error Accounted")
hold on


for err = 1:error_dev
    a.error = range_res(err);
    [t_err, s_err] = ode45(@(t,s)controlWithError(t,s,a),tspan,s0,op);
    p2dot_err = zeros(size(t));
    for i = 1:length(t)
        sdot_err = controlWithError(t(i),s(i,:),a);
        p2dot_err(i) = sdot_err(2);
    end
    deck_pos_error(err) = abs(max(s(:,1)) - max(s_err(:,1)));
    deck_acc_error(err) = abs(max(p2dot_err) - max(p_ddot));
end

plot(range_res/9.8*1000, abs(deck_pos_error)*100, DisplayName="Resolution Error")
hold on

for err = 1:error_dev
    a.error = range_res(err) + bias_instability;
    [t_err, s_err] = ode45(@(t,s)controlWithError(t,s,a),tspan,s0,op);
    p2dot_err = zeros(size(t));
    for i = 1:length(t)
        sdot_err = controlWithError(t(i),s(i,:),a);
        p2dot_err(i) = sdot_err(2);
    end
    deck_pos_error(err) = abs(max(s(:,1)) - max(s_err(:,1)));
    deck_acc_error(err) = abs(max(p2dot_err) - max(p_ddot));
end

plot(range_res/9.8*1000, abs(deck_pos_error)*100, DisplayName="Resolution and Bias Instability Error")
hold on


for err = 1:error_dev
    a.error = bias_instability;
    [t_err, s_err] = ode45(@(t,s)controlWithError(t,s,a),tspan,s0,op);
    p2dot_err = zeros(size(t));
    for i = 1:length(t)
        sdot_err = controlWithError(t(i),s(i,:),a);
        p2dot_err(i) = sdot_err(2);
    end
    deck_pos_error(err) = abs(max(s(:,1)) - max(s_err(:,1)));
    deck_acc_error(err) = abs(max(p2dot_err) - max(p_ddot));
end

plot(range_res/9.8*1000, abs(deck_pos_error)*100, DisplayName="Bias Instability Error")
hold off

yline(0.5, '--', "Positional Error Threshold")
ylim([0 0.75])
title('Position Deviation vs Time')
xlabel('Resolution (mg)')
ylabel('Position Error (cm)')
legend

%% Functions

function sdot = controlWithError(t, s, a)    

    % Current states
    p = s(1);
    pdot = s(2);
    pr_err_accum = s(3);
    
    % Error in relative position (distance to center of operation region)
    pr_err = s(1)-a.d(t)-a.pr_ref;
    
    % Magnitude of relative position PID control (largest near bounds of
    % operation region, smallest in center of operation region)
    
    % Piecewise-linear gain proportion
    
    % Radius/distace from center for zero relative position control
    r_k = a.r_k;
    h_k = a.h_k; % Baseline proportion of relative position control used
    
    k = (((h_k-1)/(r_k+a.pr_ref-1))*(abs(pr_err)-r_k)+h_k)*(abs(pr_err) > r_k) ...
            + h_k*(abs(pr_err) <= r_k);
    kp = a.kp*k; % Proportional [kg*s^-2]
    kd = a.kd*k; % Derivative   [kg/s]
    ki = a.ki*k; % Integral     [kg*s^-3]
    
    % For testing acceleration and relative position control seperately
    % kp = a.kp; % Proportional [kg*s^-2]
    % kd = a.kd; % Derivative   [kg/s]
    % ki = a.ki; % Integral     [kg*s^-3]
    
    % Magnitude of acceleration control (largest in center of operation region,
    % smallest near boudnaries of operation region)
    
    % Piecewise-linear control law
    
    % Radius/distace from center for full inertial control
    r_g = a.r_g; 
    
    c = ((-1/(a.pr_ref-r_g))*(abs(pr_err)-r_g)+1)*(abs(pr_err) > r_g) ...
                + 1*(abs(pr_err) <= r_g);
    G = a.G*c; % Acceleration gain [kg]
    H = a.H*c; % Velocity gain     [kg/s]
    
    % For testing acceleration and relative position control seperately
    % G = a.G; % Acceleration gain [kg]
    % H = a.H; % Velocity gain     [kg/s]
    
    % Derivative of states
    sdot = zeros(3,1);
    % Inertial Velocity
    sdot(1) = pdot;
    % Relative position control force
    f_pr = -(kp*(p-a.d(t)-a.pr_ref) + ki*pr_err_accum + kd*(pdot-a.ddot(t)));
    % Inertial Acceleration
    sdot(2) = (-H*pdot + f_pr - a.m*a.g - G*a.error) / (a.m + G);
    % Error in relative position
    sdot(3) = pr_err;

end

