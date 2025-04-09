% Contributors: Mohammed Al Alawi
% Course number: ASEN 4018
% File name: constants
% Created: 10/13/2024

% This struct file is to store the system paramenters that are constants.

% =================== %
platform.g = 9.81; %[m/s^2]

platform.m1 = 3.3015; %[kg]
platform.m2 = 3.3015; %[kg]

platform.Ixx1 = 0.0105; %[kg.m^2]
platform.Ixy1 = 0.0000905; %[kg.m^2]
platform.Ixz1 = 0; %[kg.m^2]
platform.Iyx1 = platform.Ixy1; %[kg.m^2]
platform.Iyy1 = 0.2925; %[kg.m^2]
platform.Iyz1 = 0; %[kg.m^2]
platform.Izx1 = platform.Ixz1; %[kg.m^2]
platform.Izy1 = platform.Iyz1; %[kg.m^2]
platform.Izz1 = 0.2925; %[kg.m^2]

platform.Ixx2 = 0.0105; %[kg.m^2]
platform.Ixy2 = 0.0000905; %[kg.m^2]
platform.Ixz2 = 0; %[kg.m^2]
platform.Iyx2 = platform.Ixy2; %[kg.m^2]
platform.Iyy2 = 0.2925; %[kg.m^2]
platform.Iyz2 = 0; %[kg.m^2]
platform.Izx2 = platform.Ixz2; %[kg.m^2]
platform.Izy2 = platform.Iyz2; %[kg.m^2]
platform.Izz2 = 0.2925; %[kg.m^2]


platform.l1 = 0.905; %[m]
platform.l2 = 0.905; %[m]

platform.r1 = platform.l1/2; %[m]
platform.r2 = platform.l2/2; %[m]


platform.Irotor1 = 0.00015; %[kg.m2]
platform.Irotor2 = 0.00015; %[kg.m2]
platform.Irotor3 = 0.00015; %[kg.m2]

platform.N1 = 50; %Gearhead ration of motor 1
platform.N2 = 50; %Gearhead ration of motor 2
platform.N3 = 50; %Gearhead ration of motor 2

platform.eta1 = 1; %Motor 1 efficiency
platform.eta2 = 1; %Motor 2 efficiency
platform.eta3 = 1; %Motor 3 efficiency

platform.mu1 = 0.08*platform.N1; %Static friction torque of joint 1 (about cg) [N.m]
platform.mu2 = 0.08*platform.N2; %Static friction torque of joint 2 (about cg) [N.m]
platform.mu3 = 0.08*platform.N3; %Static friction torque of joint 3 (about cg) [N.m]

platform.b1 = 0.1897*platform.N1; %Kinetic friction coefficient of joint 1 (damping ratio)
platform.b2 = 0.1897*platform.N2; %Kinetic friction coefficient of joint 2 (damping ratio)
platform.b3 = 0.1897*platform.N3; %Kinetic friction coefficient of joint 3 (damping ratio)