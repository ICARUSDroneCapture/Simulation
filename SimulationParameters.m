% Contributors: Mohammed Al Alawi, Kirin Kawamoto
% Course number: ASEN 4018
% Created: 4/3/2025

% This struct file is to store the system paramenters that are constants.

% =================== %
a.g = 9.81; %[m/s^2]

a.m1 = 0.39; %[kg]

a.l1 = 0.639; %[m]

a.Ixx1 = 0; %[kg.m^2]
a.Ixy1 = 0; %[kg.m^2]
a.Ixz1 = 0; %[kg.m^2]
a.Iyx1 = a.Ixy1; %[kg.m^2]
a.Iyy1 = (1/12)*a.m1*(a.l1)^2; %[kg.m^2]
a.Iyz1 = 0; %[kg.m^2]
a.Izx1 = a.Ixz1; %[kg.m^2]
a.Izy1 = a.Iyz1; %[kg.m^2]
a.Izz1 = (1/12)*a.m1*(a.l1)^2; %[kg.m^2]

a.r1 = a.l1/2; %[m]

a.Irotor1 = 0.15; %[kg.m2]

a.N1 = 50; %Gearhead ration of motor 1

a.eta1 = 0.85; %Motor 1 efficiency

a.mu1 = 1; %Static friction torque of joint 1 (about cg) [N.m]

a.b1 = 0.15; %Kinetic friction coefficient of joint 1 (damping ratio)

% moments of inertia matrices
I1 = [platform.Ixx1 platform.Ixy1 platform.Ixz1;
      platform.Iyx1 platform.Iyy1 platform.Iyz1; 
      platform.Izx1 platform.Izy1 platform.Izz1]; %link 1

% total moment of inertia
Irotor1 = [0 0 0;
           0 platform.Irotor1 0;
           0 0 0];

Itot1 = I1+(platform.N1^2)*Irotor1;