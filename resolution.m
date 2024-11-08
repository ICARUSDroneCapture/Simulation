

figure(2)

plot(t,a.d2dot(t), 'Color', [0.2 0.5 0.9 0.2])
hold on


sampling_rate = 10; %1 Hz
time_step = 1/sampling_rate;
accel_resolution = 0.005; % mg

tspan = 0:time_step:10;

lower_bound = a.d2dot(tspan) - accel_resolution;
upper_bound = a.d2dot(tspan) + accel_resolution;

error = 1 * accel_resolution;

y_prev = 0;

for i = 0:(length(tspan)-1)
    x1 = tspan(i+1);
    x2 = x1 + time_step;
    y = upper_bound(i+1);
    plot([x1; x2], [y; y], color='blue')
    if i ~= 0
        plot([x1; x1], [y; y_prev], color='blue')
    end
    hold on

    y_prev = y;
end

y_prev = 0;

for i = 0:(length(tspan)-1)
    x1 = tspan(i+1);
    x2 = x1 + time_step;
    y = lower_bound(i+1);
    plot([x1; x2], [y; y], color='blue')
    if i ~= 0
        plot([x1; x1], [y; y_prev], color='blue')
    end
    hold on

    y_prev = y;
end

title('Acceleration vs Time')
xlabel('Time (s)')
ylabel('Acceleration (m/s^2)')
legend('Deck Acceleration')


figure(3)

a.error = 2 * accel_resolution;
op = odeset('RelTol',1e-12,'AbsTol',1e-12);
[t_err, s_err] = ode45(@(t,s)controlWithError(t,s,a),tspan,s0,op);

% Maximum acceleration metric
p2dot_max = 0.005*beta^2;

% Position
subplot(1,2,1);
plot(t_err,s_err(:,1))
hold on
plot(t,a.d(t))
title('Position vs Time with Error')
xlabel('Time (s)')
ylabel('Position (m)')
legend('Platform (considering error)', 'Deck')

% Acceleration
subplot(1,2,2);
% Feeding states back through EOM to calculating inertial acceleration of
% the platfor
p2dot_err = zeros(size(t));
for i = 1:length(t)
    sdot_err = controlWithError(t(i),s(i,:),a);
    p2dot_err(i) = sdot_err(2);
end
plot(t,p2dot_err)
hold on
plot(t,a.d2dot(t))
yline(p2dot_max,'--')
yline(-p2dot_max,'--')
title('Acceleration vs Time')
xlabel('Time (s)')
ylabel('Acceleration (m/s^2)')
legend('Platform (considering error)', 'Deck')
disp(max(abs(p2dot))<p2dot_max)


figure(4)


error_dev = 100;

deck_pos_error = zeros(1, error_dev);
deck_acc_error = zeros(1, error_dev);

range_res = linspace(0.0000357, 0.0001785, error_dev);

for err = 1:error_dev
    a.error = range_res(err);
    [t_err, s_err] = ode45(@(t,s)controlWithError(t,s,a),tspan,s0,op);
    p2dot_err = zeros(size(t));
    for i = 1:length(t)
        sdot_err = controlWithError(t(i),s(i,:),a);
        p2dot_err(i) = sdot_err(2);
    end
    deck_pos_error(err) = max(s(:,1)) - max(s_err(:,1));
    deck_acc_error(err) = max(p2dot) - max(p2dot_err);
end

subplot(1, 2, 1)
plot(range_res*1000, abs(deck_pos_error)*100)
yline(0.5, '--')
title('Position Deviation vs Time')
xlabel('Resolution (mg)')
ylabel('Position Error (cm)')

subplot(1, 2, 2)
plot(range_res*1000, abs(deck_acc_error))
title('Acceleration Deviation vs Time')
xlabel('Resolution (mg)')
ylabel('Acceleration Error (m/s)')

function sdot = controlWithError(t, s, a)


% GRAB nearest timestep

% Current states
p = s(1);
pdot = s(2);
pr_err_accum = s(3);

% Error in relative position (distance to center of operation region)
pr_err = s(1)-a.d(t)-a.pr_ref;

% For testing acceleration and relative position control seperately
kp = a.kp; % Proportional
kd = a.kd; % Derivative
ki = a.ki; % Integral

% For testing acceleration and relative position control seperately
G = a.G;

% Derivative of states
sdot = zeros(3,1);
% Inertial Velocity
sdot(1) = pdot;



% Inertial Acceleration
sdot(2) = -(a.m*a.g + kd*(pdot-a.ddot(t)) + kp*(p-a.d(t)-a.pr_ref) ...
    + ki*pr_err_accum - G*a.error) / (a.m + G);
% Error in relative position
sdot(3) = pr_err;

end
