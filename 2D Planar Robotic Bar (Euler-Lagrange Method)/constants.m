% Contributors: Mohammed Al Alawi
% Course number: ASEN 4018
% File name: constants
% Created: 2/28/2025

% This struct file is to store the system paramenters that are constants.

% =================== %
platform.g = 9.81; %[m/s^2]

platform.m1 = 3.3015; %[kg]

platform.Ixx1 = 0.0105; %[kg.m^2]
platform.Ixy1 = 0.0000905; %[kg.m^2]
platform.Ixz1 = 0; %[kg.m^2]
platform.Iyx1 = platform.Ixy1; %[kg.m^2]
platform.Iyy1 = 0.2925; %[kg.m^2]
platform.Iyz1 = 0; %[kg.m^2]
platform.Izx1 = platform.Ixz1; %[kg.m^2]
platform.Izy1 = platform.Iyz1; %[kg.m^2]
platform.Izz1 = 0.2925; %[kg.m^2]

platform.l1 = 0.905; %[m]

platform.r1 = platform.l1/2; %[m]

platform.Irotor1 = 0.15; %[kg.m2]

platform.N1 = 50; %Gearhead ration of motor 1

platform.eta1 = 0.85; %Motor 1 efficiency

platform.mu1 = 1; %Static friction torque of joint 1 (about cg) [N.m]

platform.b1 = 0.15; %Kinetic friction coefficient of joint 1 (damping ratio)
