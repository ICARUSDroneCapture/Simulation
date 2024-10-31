close all; clear; clc;

% Rigid arm system with control
a.m = 1; % Mass (kg)
a.g = 9.81; % Acceleration of gravity (m/s^2)
tspan = [0 10]; % Simulation time (s)


% Disturbance equations

% Parameters
alpha = 0.5; % wave amplitdue (m)
Tmax = 7.5; % Maximum period
hdeck = 1; % inertial reference deck hight (m) (arbitrary)

% Values
k = 2;
T = Tmax / k; % Period of deck disturbance (s)
beta = (2*pi/T); % wave frequency (rad/s)

% Deck motion functions
a.d = @(t) alpha*sin(beta*t) + hdeck;
a.ddot = @(t) beta*alpha*cos(beta*t);
a.d2dot = @(t) -beta^2*alpha*sin(beta*t);


% Gains, desired position, and initial state\

% Control constants
a.G = 700; % Need at least 700
a.kd = 0;
a.kp = 0;
a.ki = 0;

% Desired deck position
a.pr_ref = 0.5; % desired relative position of platform (m)

% Initial State
s0 = [a.d(tspan(1))+a.pr_ref; a.ddot(tspan(1)); 0];

op = odeset('RelTol',1e-12,'AbsTol',1e-12);
[t, s] = ode45(@(t,s)rigidArmControl(t,s,a),tspan,s0,op);

% Maximum acceleration metric
p2dot_max = 0.005*beta^2;

% Plotting Position vs Time and Acceleration vs Time
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
% Feeding states back through EOM to calculating inertial acceleration of
% the platfor
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
close all;
disp(max(abs(p2dot))<p2dot_max)


function sdot = rigidArmControl(t, s, a)
% rigidArmControl is the EOM for the 1 DOF model of the inertially
% stabilized platform. It uses inertial acceleration control when the
% platform is close to the center of the operation region, and uses PID
% control on the relative position as the platform goes closer to the
% operational bounderies
%
% Inputs:   t    = current time
%           s    = vector of states
%                = [p; pdot; pr_err_accum] where p is the inertial position 
%                  of the platform, pdot is the inertial velocity of the 
%                  platform, and pr_err_accum is the integral of the error 
%                  in the relative position of the platform
%           a    = structure containing environmental constants and gain
%                  values
% Outputs:  sdot = time derivative of input state vector
%                = [pdot; p2dot; pr_err] where pdot is the inertial 
%                  velocity of the platform, p2dot is the inertial 
%                  acceleration of the platform,and pr_err is the error in 
%                  the relative position of the platform

% Current states
p = s(1);
pdot = s(2);
pr_err_accum = s(3);

% Error in relative position (distance to center of operation region)
pr_err = s(1)-a.d(t)-a.pr_ref;

% Magnitude of relative position PID control (largest near bounds of
% operation region, smallest in center of operation region)
% kp = a.kp*abs(pr_err) / 0.5; % Proportional
% kd = a.kd*abs(pr_err) / 0.5; % Derivative
% ki = a.ki*abs(pr_err) / 0.5; % Integral

% For testing acceleration and relative position control seperately
kp = a.kp; % Proportional
kd = a.kd; % Derivative
ki = a.ki; % Integral

% Magnitude of acceleration control (largest in center of operation region,
% smallest near boudnaries of operation region)
% G = 0;
% if (abs(pr_err) < 0.5)
%     G = a.G * (1 - abs(pr_err) / 0.5);
% end

% For testing acceleration and relative position control seperately
G = a.G;

% Derivative of states
sdot = zeros(3,1);
% Inertial Velocity
sdot(1) = pdot;
% Inertial Acceleration
sdot(2) = -(a.m*a.g + kd*(pdot-a.ddot(t)) + kp*(p-a.d(t)-a.pr_ref) ...
    + ki*pr_err_accum) / (a.m + G);
% Error in relative position
sdot(3) = pr_err;

end