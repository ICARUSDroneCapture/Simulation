% Contributors: Mohammed Al Alawi, Kirin Kawamoto
% Course number: ASEN 4018
% Created: 4/3/2025

% This struct file is to store the system paramenters that are constants.

% =================== %
a.g = 9.81; %[m/s^2]

% a.m1 = 0.39; %[kg]
% 
% a.l1 = 0.639; %[m]

a.m1 = 3.3015; %[kg]
a.l1 = 0.905; %[m]

Ixx1 = 0; %[kg.m^2]
Ixy1 = 0.0000905; %[kg.m^2]
Ixz1 = 0; %[kg.m^2]
Iyx1 = Ixy1; %[kg.m^2]
% Iyy1 = (1/12)*a.m1*(a.l1)^2; %[kg.m^2]
Iyy1 = 0.2925; %[kg.m^2]
Iyz1 = 0; %[kg.m^2]
Izx1 = Ixz1; %[kg.m^2]
Izy1 = Iyz1; %[kg.m^2]
% Izz1 = (1/12)*a.m1*(a.l1)^2; %[kg.m^2]
Izz1 = 0.2925; %[kg.m^2]

a.r1 = a.l1/2; %[m]

a.Irotor1 = 0.15; %[kg.m2]

N1 = 50; %Gearhead ration of motor 1

eta1 = 0.85; %Motor 1 efficiency

mu1 = 1; %Static friction torque of joint 1 (about cg) [N.m]

b1 = 0.15; %Kinetic friction coefficient of joint 1 (damping ratio)

% moments of inertia matrices
I1 = [Ixx1 Ixy1 Ixz1;
      Iyx1 Iyy1 Iyz1; 
      Izx1 Izy1 Izz1]; %link 1

% total moment of inertia
Irotor1 = [0 0 0;
           0 a.Irotor1 0;
           0 0 0];

a.Itot1 = I1+(N1^2)*Irotor1;

% Low-Pass filter on measured acceleration
a.f_c = 50; % Cutoff frequency [Hz]
a.omega = 2*pi*a.f_c; % Angular frequency [rad/s]

% important matrices
a.MU = mu1; %static friction matrix
a.B = b1; %kinetic friction matrix
a.N = N1*eta1; %gearhead ratio/efficiency matrix

a.thetad = @(t) (pi*sin((pi*t)/15).^2)/2;
a.thetad_dot = @(t) (pi^2*cos((pi*t)/15).*sin((pi*t)/15))/15;
a.thetad_ddot = @(t) (pi^3*cos((pi*t)/15).^2)/225 - (pi^3*sin((pi*t)/15).^2)/225;