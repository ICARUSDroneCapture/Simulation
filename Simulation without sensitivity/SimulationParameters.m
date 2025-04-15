% Contributors: Mohammed Al Alawi, Kirin Kawamoto
% Course number: ASEN 4018
% Created: 4/3/2025

% This struct file is to store the system paramenters that are constants.

% =================== %
a.g = 9.81; %[m/s^2]

a.m1 = 3.08; %[kg]

a.l1 = 0.6; %[m]

% a.m1 = 3.3015; %[kg]
% a.l1 = 0.905; %[m]

Ixx1 = 0; %[kg.m^2]
Ixy1 = 0.0000905; %[kg.m^2]
Ixz1 = 0; %[kg.m^2]
Iyx1 = Ixy1; %[kg.m^2]
Iyy1 = (1/12)*a.m1*(a.l1)^2; %[kg.m^2]
% Iyy1 = 0.2925; %[kg.m^2]
Iyz1 = 0; %[kg.m^2]
Izx1 = Ixz1; %[kg.m^2]
Izy1 = Iyz1; %[kg.m^2]
Izz1 = (1/12)*a.m1*(a.l1)^2; %[kg.m^2]
% Izz1 = 0.2925; %[kg.m^2]

a.r1 = a.l1/2; %[m]

a.Irotor1 = 0.00015; %[kg.m2]

N1 = 50; %Gearhead ration of motor 1

eta1 = 1; %Motor 1 efficiency

mu1 = 0.08; %Static friction torque of joint 1 (about cg) [N.m]

b1 = 0.1897; %Kinetic friction coefficient of joint 1 (damping ratio)

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
a.BF = b1; %kinetic friction matrix
a.N = N1*eta1; %gearhead ratio/efficiency matrix

a.thetad = @(t) (pi*sin((pi*t)/15).^2)/2;
a.thetad_dot = @(t) (pi^2*cos((pi*t)/15).*sin((pi*t)/15))/15;
a.thetad_ddot = @(t) (pi^3*cos((pi*t)/15).^2)/225 - (pi^3*sin((pi*t)/15).^2)/225;

% Inertial Stabilization Control
a.ka = 0.5;  % Acceleration Control [kg]

% Relative Position Control at center
scale = 1;
a.kp_c = scale*0.05;  % Proportional [kg*s^-2]
a.kd_c = scale*5;  % Derivative [kg/s]    
a.ki_c = scale*0.5;  % Integral [kg*s^-3] 

% Relative Position Control at boundaries
a.kp_b = 50;  % Proportional [kg*s^-2]
a.kd_b = 10;  % Derivative [kg/s]    
a.ki_b = 0.5;  % Integral [kg*s^-3]

% a.kp_b = 0;  % Proportional [kg*s^-2]
% a.kd_b = 0;  % Derivative [kg/s]    
% a.ki_b = 0;  % Integral [kg*s^-3]

a.q1_ref = -pi/4;

% -------------------------- Parameters -------------------------------- %
a.w = pi; % Range of inputs
d = a.q1_ref; % Center of input region

% Piecewise radii
c = 0.8;
r_c = 45*c; % Full isolation control radius
b = 0.8;
r_b = 45*b; % Zero relative position control radius

% Polynomial order
n = 1;


% --------------------------- Mixing Functions ------------------------- %

% Center gain scale (inertial isolation)
a_c = -1 / abs(r_c - a.w/2)^n;
k_c = 1;
a.C = @(x) 1 .* (abs(x - d) <= r_c) ...
       + (a_c*abs(x-(d+sign(x-d)*r_c)).^n + k_c) .* (abs(x - d) > r_c & ...
       abs(x - d) <= a.w/2);

% Boundary gain scale (relative position control)
a_b = 1 / abs(r_b - a.w/2)^n;
k_b = 0;
a.B = @(x) (a_b*abs(x-(d+sign(x-d)*r_b)).^n + k_b) .* (abs(x - d) > r_b & ...
       abs(x - d) <= a.w/2) ...
       + 1 .* (abs(x - d) > a.w/2);
