close all; clear; clc;

% Rigid arm system with control
a.m = 1;
a.g = 9.81;
tspan = [0 10];

% Disturbance equations
alpha = 0.5; % wave amplitdue (m)
hdeck = 1; % inertial reference deck hight (m)
k = 2;
beta = k*(2*pi/7.5); % wave frequency (rad/s)

% Maximum acceleration metric
p2dot_max = 0.005*beta^2;

a.d = @(t) alpha*sin(beta*t) + hdeck;
a.ddot = @(t) beta*alpha*cos(beta*t);
a.d2dot = @(t) -beta^2*alpha*sin(beta*t);

% Control constants
a.G = 1000;
a.b = 100;
a.k = 1000;
a.ki = 0;
% a.G = 100;
% a.b = 40;
% a.k = 100;
% a.ki = 15;

a.pr_ref = 0.5;

s0 = [a.d(tspan(1))+a.pr_ref; a.ddot(tspan(1)); 0];

op = odeset('RelTol',1e-12,'AbsTol',1e-12);
[t, s] = ode45(@(t,s)rigidArmControl(t,s,a),tspan,s0,op);

figure;

% Position
subplot(1,2,1);
plot(t,s(:,1))
hold on
plot(t,a.d(t))
title('Position vs Time')
xlabel('Time (s)')
ylabel('Position (m)')
legend('Platform', 'Deck')

% Acceleration
subplot(1,2,2);
p2dot = zeros(size(t));
for i = 1:length(t)
    sdot = rigidArmControl(t(i),s(i,:),a);
    p2dot(i) = sdot(2);
end
plot(t,p2dot)
hold on
plot(t,a.d2dot(t))
yline(p2dot_max,'--')
yline(-p2dot_max,'--')
title('Acceleration vs Time')
xlabel('Time (s)')
ylabel('Acceleration (m/s^2)')
legend('Platform', 'Deck')


pr = s(:,1) - a.d(t);
figure;
plot(t, pr)


function sdot = rigidArmControl(t, s, a)

pr_err = s(1)-a.d(t)-a.pr_ref;

kp = a.k*abs(pr_err) / 0.5;
kd = a.b*abs(pr_err) / 0.5;
ki = a.ki*abs(pr_err) / 0.5;
G = 0;
if (abs(pr_err) < 0.5)
    G = a.G * (1 - abs(pr_err) / 0.5);
end

sdot = zeros(3,1);
sdot(1) = s(2);
sdot(2) = -(a.m*a.g + kd*(s(2)-a.ddot(t)) + kp*(s(1)-a.d(t)-a.pr_ref) + ki*(s(1)-a.d(t)-a.pr_ref)) / (a.m + G);
sdot(3) = s(1)-a.d(t)-a.pr_ref;

end