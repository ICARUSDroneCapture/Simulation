function dS = EOM_V3_5(t, S, a)
%EOM_V3.5 Summary of this function goes here
%   Detailed explanation goes here

q1 = S(1); % Joint angle
q1_dot = S(2); % Joint velocity
q1_err_accum = S(3); % Error in desired angle integrationpm_ddot = S(4:6); % Measured acceleration

Iyy = a.Itot1(2,2);

%%% Control Law %%%

C = a.C(q1);
ka = C*a.ka;

% Proportion of relative position control
B = a.B(q1);
kp = a.kp_c*C + a.kp_b*B; % Proportional 
kd = a.kd_c*C + a.kd_b*B; % Derivative   
ki = a.ki_c*C + a.ki_b*B; % Integral

q1_err = q1 - a.q1_ref;
tau_r = -(kp*q1_err + ki*q1_err_accum + kd*q1_dot);

% friction

tau_f = -a.MU*sign(q1_dot) - a.BF*q1_dot; 

%%% EOM %%%

tau_c = tau_r + tau_f;

% q1_ddot = -(a.m1*a.thetad_ddot(t)*a.r1^2 ...
%     - a.g*a.m1*cos(q1 + a.thetad(t))*a.r1 + a.N*a.BF*q1_dot - a.N*tau ...
%     + a.N*a.MU*sign(q1_dot) + Iyy*a.thetad_ddot(t))/(a.m1*a.r1^2 + Iyy);

% q1_ddot = (a.N*tau - (a.m1*a.r1^2 + Iyy)*a.thetad_ddot(t) ...
%     + a.g*a.m1*a.r1*cos(q1 + a.thetad(t)))/(a.m1*a.r1^2 + Iyy);

q1_ddot = (2*a.N*tau_c - 2*Iyy*a.thetad_ddot(t) - ...
    2*a.m1*a.r1^2*a.thetad_ddot(t) - a.N*ka*a.l1^2*a.thetad_ddot(t) ...
    + 2*a.g*a.m1*a.r1*cos(q1 + a.thetad(t)) + a.N*ka*a.l1^2*sin(2*q1 ...
    + 2*a.thetad(t))*q1_dot^2 + a.N*ka*a.l1^2*sin(2*q1 ...
    + 2*a.thetad(t))*a.thetad_dot(t)^2 - a.N*ka*a.l1^2*cos(2*q1 ...
    + 2*a.thetad(t))*a.thetad_ddot(t) + 2*a.N*ka*a.l1^2*sin(2*q1 ...
    + 2*a.thetad(t))*q1_dot*a.thetad_dot(t))/(2*Iyy + 2*a.m1*a.r1^2 ...
    + a.N*ka*a.l1^2 + a.N*ka*a.l1^2*cos(2*q1 + 2*a.thetad(t)));

dS = [q1_dot;
      q1_ddot;
      q1_err];

end

