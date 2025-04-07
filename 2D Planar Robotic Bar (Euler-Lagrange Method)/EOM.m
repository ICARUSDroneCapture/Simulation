function dS = EOM(t, S, a)
%EOM Summary of this function goes here
%   Detailed explanation goes here

theta = S(1);
theta_dot = S(2);

tau = 0;

theta_ddot = 0;

dS(1) = theta_dot;
dS(2) = theta_ddot;

end

