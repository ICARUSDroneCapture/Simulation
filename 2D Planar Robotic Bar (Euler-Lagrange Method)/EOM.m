function dS = EOM(t, S, a)
%EOM Summary of this function goes here
%   Detailed explanation goes here

q1 = S(1);
q1_dot = S(2);

p_ddot = [-l1*sin(q1(t) + theta2(t))*(diff(q1(t), t, t) ...
          + diff(theta2(t), t, t)) - l1*cos(q1(t) ...
          + theta2(t))*(diff(q1(t), t) + diff(theta2(t), t))^2;
           0;
           l1*sin(q1(t) + theta2(t))*(diff(q1(t), t) ...
         + diff(theta2(t), t))^2 - l1*cos(q1(t) + theta2(t)) ...
         *(diff(q1(t), t, t) + diff(theta2(t), t, t))];

tau = 0;

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

theta_ddot = M^(-1) * (-M_dstar*a.thetad_ddot - (C+C_dstar-j_q)*q1_dot ...
    - (C_star-j_star_q)*a.thetad_dot - a.MU*sign(q1_dot) - a.B*q1_dot + N*tau);

EOM = diff(j,'t')+M*diff(q1(t),'t',2)+M_dstar*diff(theta_D,'t',2) ...
    + (C+C_dstar-j_q)*diff(q1(t),'t')+(C_star-j_star_q)*diff(theta_D,'t')...
    + MU*sign(diff(q1(t),'t'))+B*diff(q1(t),'t')+G == N*tau1(t); 

dS(1) = q1_dot;
dS(2) = q1_ddot;

end

