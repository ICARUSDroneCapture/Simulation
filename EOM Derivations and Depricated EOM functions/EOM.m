function dS = EOM(t, S, a)
%EOM Summary of this function goes here
%   Detailed explanation goes here

q1 = S(1); % Joint angle
q1_dot = S(2); % Joint velocity
q1_err_accum = S(3); % Error in desired angle integration
pm_dot = S(4:6); % Measured velocity
pm_ddot = S(7:9); % Measured acceleration

% Equation of motion for bar
dJdt = 0;

Ixy = a.Itot1(1,2);
Iyy = a.Itot1(2,2);
Iyz = a.Itot1(3,2);

M = a.m1*a.r1^2 + Iyy;

M_dstar = [Ixy*cos(q1 + a.thetad(t)) + Iyz*sin(q1 + a.thetad(t)), ...
           a.m1*a.r1^2 + Iyy, ...
           Iyz*cos(q1 + a.thetad(t)) - Ixy*sin(q1 + a.thetad(t))];

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
      0 0 300*1]*0;
Kv = [0 0 0;
      0 0 0;
      0 0 375*1]*0;

% gains of motor 1 (joint)
Kp1 = 4*0;
Ki1 = 10*0;
Kd1 = 15*0;

R_I_B = [cos(a.thetad(t)) 0 -sin(a.thetad(t));
        0 1 0;
        sin(a.thetad(t)) 0 cos(a.thetad(t))];

F_I = -R_I_B*(Ka*pm_ddot + Kv*pm_dot);

% the jocabain
J = [-a.l1*sin(q1); 0;- a.l1*cos(q1)];

tau_I = J.'*F_I;

q1_err = q1 - q1_ref;
tau_r = -(Kp1*q1_err + Ki1*(q1_err_accum) + Kd1*(q1_dot));

tau = tau_I + tau_r;

%%% EOM %%%

thetad_ddot = [0; a.thetad_ddot(t); 0];
thetad_dot = [0; a.thetad_dot(t); 0];

q1_ddot = M^(-1) * (-dJdt - M_dstar*thetad_ddot ...
    - (C+C_dstar-j_q)*q1_dot - (C_star-j_star_q)*thetad_dot ...
    - a.MU*sign(q1_dot) - a.B*q1_dot - G + a.N*tau);

% Calculating resulting platform acceleration
p_ddot = [(-a.l1*sin(q1 + a.thetad(t))*(q1_ddot + a.thetad_ddot(t)) ...
          -a.l1*cos(q1 + a.thetad(t))*(q1_dot + a.thetad_dot(t))^2);
           0;
           (a.l1*sin(q1 + a.thetad(t))*(q1_dot + a.thetad_dot(t))^2 ...
          -a.l1*cos(q1 + a.thetad(t)) *(q1_ddot + a.thetad_ddot(t)))];

% Derivative states
% dS(1) = q1_dot;
% dS(2) = q1_ddot;
% dS(3) = q1_err;
% dS(4:6) = pm_ddot;
% dS(7:9) = a.omega*(p_ddot - pm_ddot);

pm_dddot = a.omega*(p_ddot - pm_ddot);

dS = [q1_dot;
      q1_ddot;
      q1_err;
      pm_ddot;
      pm_dddot];

end

