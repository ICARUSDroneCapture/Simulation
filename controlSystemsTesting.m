close all; clear; clc;

% Rigid arm system with control
a.m = 1; % Mass (kg)
a.g = 9.81; % Acceleration of gravity (m/s^2)
tspan = [0 50]; % Simulation time (s)


% Disturbance equations

% Parameters
alpha = 0.4; % wave amplitdue (m)
Tmax = 7.5; % Maximum period
hdeck = 1; % inertial reference deck hight (m) (arbitrary)

% Wave frequency
k = 1;
T = Tmax / k; % Period of deck disturbance (s)
beta = (2*pi/T); % wave frequency (rad/s)

% Deck motion functions
a.d = @(t) alpha*sin(beta*t) + hdeck; % [m]
a.ddot = @(t) beta*alpha*cos(beta*t); % [m/s]
a.d2dot = @(t) -beta^2*alpha*sin(beta*t); % [m*s^-2]


% Control constants

% Inertial Stabilization Control
a.G = 700;  % Need at least 2800; Acceleration Control [kg]
a.H = 5000;  % Need at least  (look at G/H to get 2s settling time)  ; Velocity Control [kg/s]
% Relative Position Control
a.kp = 3000; % Need at least 500; Proportional [kg*s^-2]
a.kd = 500; % Need at least 200; Derivative [kg/s]
a.ki = 200; % Need at least 200; Integral [kg*s^-3]

% Radius from center of operation region for full inertial control
a.r_g = 0.4; % [m]
% Radius from center of operation region for base-line proportion of
% relative position control
a.r_k = 0.4; % [m]
a.h_k = 0.05; % Base-line proportion

% Desired deck position
a.pr_ref = 0.5; % desired relative position of platform (m)

% Initial States
p0 =  a.d(tspan(1))+a.pr_ref; % Initial inertial platform position
pdot0 = a.ddot(tspan(1));  % Initial inertial platform velocity
pr_err_accum0 = 0;
s0 = [p0; pdot0; pr_err_accum0];

op = odeset('RelTol',1e-12,'AbsTol',1e-12);
[t, s] = ode45(@(t,s)rigidArmControl(t,s,a),tspan,s0,op);

% Maximum acceleration metric
p2dot_max = 0.005*beta^2;

% Plotting Position vs Time and Acceleration vs Time
figure;
sgtitle('Non-Zero Relative Positon Control during Inertial Control')

% Position
subplot(1,2,1);
plot(t,s(:,1))
hold on
plot(t,a.d(t))
title('Inertial Position vs Time')
xlabel('Time (s)')
ylabel('Position (m)')
legend('Platform', 'Deck','Location','southeast')

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
title('Inertial Acceleration vs Time')
xlabel('Time (s)')
ylabel('Acceleration (m/s^2)')
legend('Platform', 'Deck','Location','southeast')


% Plotting relative position
figure;
plot(t,s(:,1)-a.d(t))
hold on
% Plotting inertial control region
x = [tspan, flip(tspan)];
yf = [a.pr_ref-a.r_g, a.pr_ref-a.r_g, a.pr_ref+a.r_g, a.pr_ref+a.r_g];
yta = [a.pr_ref+a.r_g, a.pr_ref+a.r_g, 1, 1];
ytb = [0, 0, a.pr_ref-a.r_g, a.pr_ref-a.r_g];
fill(x,yf,'y','FaceAlpha',0.2,'EdgeColor','none')
fill(x,yta,'g','FaceAlpha',0.2,'EdgeColor','none')
fill(x,ytb,'g','FaceAlpha',0.2,'EdgeColor','none')
% Plotting Relative position control region
x = [tspan, flip(tspan)];
yta = [a.pr_ref+a.r_k, a.pr_ref+a.r_k, 1, 1];
ytb = [0, 0, a.pr_ref-a.r_k, a.pr_ref-a.r_k];
fill(x,yta,'b','FaceAlpha',0.2,'EdgeColor','none')
fill(x,ytb,'b','FaceAlpha',0.2,'EdgeColor','none')
title('Relative Position vs Time')
xlabel('Time (s)')
ylabel('Position (m)')
legend('','Full Inertial','Transitional Inertial','','Transitional Relative Position','')


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
sdot(2) = (-H*pdot + f_pr - a.m*a.g) / (a.m + G);
% Error in relative position
sdot(3) = pr_err;

end