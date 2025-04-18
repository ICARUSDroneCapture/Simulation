function p_ddot = realAccel(t, q1, q1_dot, q1_ddot, a)
%REALACCEL Summary of this function goes here
%   Detailed explanation goes here

% Calculating resulting platform acceleration
p_ddot = [(-a.l1*sin(q1 + a.thetad(t))*(q1_ddot + a.thetad_ddot(t)) ...
          -a.l1*cos(q1 + a.thetad(t))*(q1_dot + a.thetad_dot(t))^2);
           0;
           (a.l1*sin(q1 + a.thetad(t))*(q1_dot + a.thetad_dot(t))^2 ...
          -a.l1*cos(q1 + a.thetad(t)) *(q1_ddot + a.thetad_ddot(t)))];

end

