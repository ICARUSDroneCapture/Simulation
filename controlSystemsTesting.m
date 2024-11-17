close all; clear; clc;

simulationParameters;

% Simulation time
tspan = [0 20]; % [s]

% Initial States
p0 =  a.d(tspan(1))+a.pr_d;   % Platform position [m]
p_dot0 = a.d_dot(tspan(1));   % Platform velocity [m/s]
pr_err_accum0 = 0;            % Integral of relative position error [m*s]
pm0 = p0;                     % Platform inetegrated position [m]
pm_dot = p_dot0;              % Platform integrated velocity [m/s]
pm_ddot = a.d_ddot(tspan(1)); % Platform measured acceleration [m*s^-2]

s0 = [p0; p_dot0; pr_err_accum0; pm0; pm_dot; pm_ddot];

% Running Simulation
op = odeset('RelTol',1e-3,'AbsTol',1e-6); % Tolerance options
[t, s] = ode45(@(t,s)rigidArmControl(t,s,a),tspan,s0,op);

% Plotting Position vs Time and Acceleration vs Time
figure;
sgtitle('Non-Zero Relative Positon Control during Inertial Control')

% Position
subplot(1,3,1);
plot(t,s(:,1))
hold on
plot(t,a.d(t))
title('Inertial Position vs Time')
xlabel('Time (s)')
ylabel('Position (m)')
legend('Platform', 'Deck','Location','southeast')

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
% Feeding states back through EOM to calculating inertial acceleration of
% the platfor
p_ddot = zeros(size(t));
for i = 1:length(t)
    s_dot = rigidArmControl(t(i),s(i,:),a);
    p_ddot(i) = s_dot(2);
end
plot(t,p_ddot)
hold on
plot(t,a.d_ddot(t))
yline(p_ddot_max,'--')
yline(-p_ddot_max,'--')
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

% Error in relative position (distance to center of operation region)
pr = pm-a.d(t);
pr_err = pr-a.pr_d;

% For testing gains without mixing proportions
ka = a.ka; % Acceleration [kg]
kv = a.kv; % Velocity     [kg/s]
ks = a.ks; % Position     [kg*s^-2]
kp = a.kp; % Proportional [kg*s^-2]
kd = a.kd; % Derivative   [kg/s]
ki = a.ki; % Integral     [kg*s^-3]

% Control gain proportions

% I = a.I(pr); % Proportion of inertial stability control to apply
% ka = a.ka*I; % Acceleration [kg]
% kv = a.kv*I; % Velocity     [kg/s]
% ks = a.ks*I; % Position     [kg*s^-2]
% 
% k = a.K(pr); % Proportion of relative position control to apply
% kp = a.kp*k; % Proportional [kg*s^-2]
% kd = a.kd*k; % Derivative   [kg/s]
% ki = a.ki*k; % Integral     [kg*s^-3]

% Derivative of states
s_dot = zeros(6,1);

% Derivative of position
s_dot(1) = p_dot;  % inertial velocity
s_dot(4) = pm_dot; % measured inertial velocity

% Control Law

% Inertial stability control force
f_i = -(ka*pm_ddot + kv*pm_dot + ks*pm);
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

t;

end