function dS = EOM(t, S, a)
%EOM Summary of this function goes here
%   Detailed explanation goes here

q1 = S(1); % Joint angle
q1_dot = S(2); % Joint velocity
q1_err_accum = S(3); % Error in desired angle integration
pm_ddot = S(4:6); % Measured acceleration

% Equation of motion for bar
dJdt = 0;

M = a.m1*a.r1^2*cos(q1 + a.thetad(t))^2 ...
    + a.m1*a.r1^2*sin(q1 + a.thetad(t))^2 + a.Itot1;

M_dstar = [0, a.m1*a.r1^2*cos(q1 + a.thetad(t))^2 ...
               + a.m1*a.r1^2*sin(q1 + a.thetad(t))^2 + a.Itot1, 0];
C = 0;
C_dstar = 0;
j_q = 0;
 
C_star = [0, 0, 0];
j_star_q = [0, 0, 0];

G = -a.g*a.m1*a.r1*cos(q1 + a.thetad(t));

%%% Control Law %%%

q1_ref = -pi/4;

% gains of inertial stability control (cartesian)
Ka = [0 0 0; 
      0 0 0; 
      0 0 300*1];
Kv = [0 0 0;
      0 0 0;
      0 0 375*1];

% gains of motor 1 (joint)
Kp1 = 4*0;
Ki1 = 10*0;
Kd1 = 15*0;

F_I = -R_I_B*(Ka*pm_ddot + Kv*p_dot);

tau_I = J.'*F_I;

tau_r = -(Kp1*(q1 - q1_ref) + Ki1*(q1_err_accum) + Kd1*(q1_dot));

tau = tau_I + tau_r;

%%% EOM %%%

q1_ddot = M^(-1) * (-dJdt - M_dstar*a.thetad_ddot - (C+C_dstar-j_q)*q1_dot ...
    - (C_star-j_star_q)*a.thetad_dot - a.MU*sign(q1_dot) - a.B*q1_dot - G ...
    + N*tau);

% Calculating resulting platform acceleration
p_ddot = [-a.l1*sin(q1 + a.thetad(t))*(q1_ddot + a.thetad_ddot(t)) ...
          -a.l1*cos(q1 + a.thetad(t))*(q1_dot + a.thetad_dot(t))^2;
           0;
           a.l1*sin(q1 + a.thetad(t))*(q1_dot + a.thetad_dot(t))^2 ...
          -a.l1*cos(q1 + a.thetad(t)) *(q1_ddot + a.thetad_ddot(t))];

% Derivative states
dS(1) = q1_dot;
dS(2) = q1_ddot;
dS(3) = q1_err;
dS(4:6) = a.omega*(p_ddot - pm_ddot);

end

