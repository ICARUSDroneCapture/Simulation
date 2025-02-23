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

C = a.C(pr)*c_i; % Proportion of inertial stability control to apply
ka = a.ka*C; % Acceleration [kg]
kv = a.kv*C; % Velocity     [kg/s]
ks = a.ks*C; % Position     [kg*s^-2]

% Proportion of relative position control
k = max(a.K(pr),c_k); 
kp = a.kp_c*C + a.kp_b*k; % Proportional [kg*s^-2]
kd = a.kd_c*C + a.kd_b*k; % Derivative   [kg/s]
ki = a.ki_c*C + a.ki_b*k; % Integral     [kg*s^-3]

% Derivative of states
s_dot = zeros(3,1);

% Inertial Velocity
s_dot(1) = p_dot;

% Relative position control force
f_pr = -(kp*pr_err + ki*pr_err_accum + kd*(p_dot-a.d_dot(t)));

% Inertial Acceleration
s_dot(2) = (-kv*p_dot - ks*p + f_pr - a.m*a.g + f_comp) / (a.m + ka);
% abs(s_dot(2)) < 1e-3 && p_dot > 1e-3
% if ( (abs(s_dot(2)) < 1e-7) && (abs(p_dot) > 1e-8) && (I == 1))
%     f_comp = -kv*p_dot + f_comp;
% end

% Error in relative position
s_dot(3) = pr_err;


end