% Contributors: Mohammed Al Alawi
% Course number: ASEN 4018
% File name: forwardKinematics.m
% Created: 10/19/2024


function F = forwardKinematicsFun(theta, desiredPosition)
% This function states the 3D nonlinear forward kinematics of the a 3D arm.
% Inputs:

% Outputs:
%calling constants
constants;
F(1) = desiredPosition(1) - (platform.l1*cos(theta(1))+platform.l2*cos(theta(1)+theta(2)))*cos(theta(3));
F(2) = desiredPosition(2) - (platform.l1*cos(theta(1))+platform.l2*cos(theta(1)+theta(2)))*sin(theta(3));
F(3) = desiredPosition(3) + platform.l1*sin(theta(1))+platform.l2*sin(theta(1)+theta(2));
end

%solving the inverse kinematic problem
desiredPosition = [0.5,0,0.9]; %[m]
fun = @(theta)forwardKinematicsFun(theta,desiredPosition);
theta0 = [0,0,0];
theta = fsolve(fun,theta0)

