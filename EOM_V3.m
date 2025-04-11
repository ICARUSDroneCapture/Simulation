function dS = EOM_V3(t, S, a)
%EOM Summary of this function goes here
%   Detailed explanation goes here

q1 = S(1); % Joint angle
q1_dot = S(2); % Joint velocity
q1_err_accum = S(3); % Error in desired angle integration
pm_ddot = S(4:6); % Measured acceleration

Iyy = a.Itot1(2,2);
p_dot = [-a.l1*sin(q1 + a.thetad(t))*(q1_dot + a.thetad_dot(t));
          0;
         -a.l1*cos(q1 + a.thetad(t))*(q1_dot + a.thetad_dot(t))];

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

R_I_B = [cos(a.thetad(t)) 0 -sin(a.thetad(t));
        0 1 0;
        sin(a.thetad(t)) 0 cos(a.thetad(t))];

F_I = -R_I_B*(Ka*pm_ddot + Kv*p_dot);

% the jocabain
J = [-a.l1*sin(q1); 0;- a.l1*cos(q1)];

tau_I = J.'*F_I;

q1_err = q1 - q1_ref;
tau_r = -(Kp1*q1_err + Ki1*(q1_err_accum) + Kd1*(q1_dot));

tau = tau_I + tau_r;

%%% EOM %%%

q1_ddot = -(a.m1*a.thetad_ddot(t)*a.r1^2 ...
    - a.g*a.m1*cos(q1 + a.thetad(t))*a.r1 + a.B*q1_dot - a.N*tau ...
    + a.MU*sign(q1_dot) + Iyy*a.thetad_ddot(t))/(a.m1*a.r1^2 + Iyy);

% Calculating resulting platform acceleration
p_ddot = [(-a.l1*sin(q1 + a.thetad(t))*(q1_ddot + a.thetad_ddot(t)) ...
          -a.l1*cos(q1 + a.thetad(t))*(q1_dot + a.thetad_dot(t))^2);
           0;
           (a.l1*sin(q1 + a.thetad(t))*(q1_dot + a.thetad_dot(t))^2 ...
          -a.l1*cos(q1 + a.thetad(t)) *(q1_ddot + a.thetad_ddot(t)))];

pm_dddot = a.omega*(p_ddot-pm_ddot);

dS = [q1_dot;
      q1_ddot;
      q1_err;
      pm_dddot];

end

