close all; clear; clc;

simulationParameters;

for i = 1:10

a.m = i; % Mass [kg]
f_comp = m0*a.g; % Gravity compensation force mass [N]

% Simulation time
tspan = [0 30]; % [s]

% Initial States
p0 =  a.d(tspan(1))+a.pr_d; % Platform position [m]
p_dot0 = a.d_dot(tspan(1)); % Platform velocity [m/s]
pr_err_accum0 = 0;          % Integral of error in relative position [m*s]
s0 = [p0; p_dot0; pr_err_accum0];

% Running Simulation
op = odeset('RelTol',1e-12,'AbsTol',1e-12); % Tolerance options
[t, s] = ode45(@(t,s)rigidArmControl(t,s,a),tspan,s0,op);

pf = s((t>25),1);
error = (max(pf)-min(pf))*1e2

% Plotting Position vs Time and Acceleration vs Time

% Position
figure(2);
hold on
plot(t,s(:,1))
hold off

% Acceleration
figure(3);
hold on
% Feeding states back through EOM to calculating inertial acceleration of
% the platfor
p_ddot = zeros(size(t));
f_comp = m0*a.g; % Gravity compensation force mass [N]
for i = 1:length(t)
    s_dot = rigidArmControl(t(i),s(i,:),a);
    p_ddot(i) = s_dot(2);
end
plot(t,p_ddot)
hold off

end

figure(2);
hold on
plot(t,a.d(t))
title('Inertial Position vs Time')
xlabel('Time (s)')
ylabel('Position (m)')
% legend('Platform', 'Deck','Location','southeast')
hold off

figure(3);
hold on
plot(t,a.d_ddot(t))
yline(p_ddot_max,'--')
yline(-p_ddot_max,'--')
% ylim([min(a.d_ddot(t))*1.25 max(a.d_ddot(t))*1.25])
title('Inertial Acceleration vs Time')
xlabel('Time (s)')
ylabel('Acceleration (m/s^2)')
% legend('Platform', 'Deck','Location','southeast')
hold off



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

% Compensation force to account for interfering control forces
% (When control forces cause non-zero steady state velocity with zero
% acceleration)
global f_comp

% Current states
p = s(1);
p_dot = s(2);
pr_err_accum = s(3);

% Error in relative position (distance to center of operation region)
pr = p-a.d(t);
pr_err = pr-a.pr_d;

% For testing gains without mixing proportions
% I = 1;
% ka = a.ka; % Acceleration [kg]
% kv = a.kv; % Velocity     [kg/s]
% ks = a.ks; % Position     [kg*s^-2]
% kp = a.kp; % Proportional [kg*s^-2]
% kd = a.kd; % Derivative   [kg/s]
% ki = a.ki; % Integral     [kg*s^-3]

% Control gain proportions
c_i = a.int_scale_i(t); % Initial scale of intertial stability gains
c_k = a.int_scale_k(t); % Initial scale of relative position gains

I = a.I(pr)*c_i; % Proportion of inertial stability control to apply
ka = a.ka*I; % Acceleration [kg]
kv = a.kv*I; % Velocity     [kg/s]
ks = a.ks*I; % Position     [kg*s^-2]

k = max(a.K(pr),c_k); % Proportion of relative position control to apply
kp = a.kp*k; % Proportional [kg*s^-2]
kd = a.kd*k; % Derivative   [kg/s]
ki = a.ki*k; % Integral     [kg*s^-3]

% Derivative of states
s_dot = zeros(3,1);

% Inertial Velocity
s_dot(1) = p_dot;

% Relative position control force
f_pr = -(kp*pr_err + ki*pr_err_accum + kd*(p_dot-a.d_dot(t)));

% Inertial Acceleration
s_dot(2) = (-kv*p_dot - ks*p + f_pr - a.m*a.g + f_comp) / (a.m + ka);
% abs(s_dot(2)) < 1e-3 && p_dot > 1e-3
if ( (abs(s_dot(2)) < 1e-7) && (abs(p_dot) > 1e-8) && (I == 1))
    f_comp = -kv*p_dot + f_comp;
end

% Error in relative position
s_dot(3) = pr_err;


end