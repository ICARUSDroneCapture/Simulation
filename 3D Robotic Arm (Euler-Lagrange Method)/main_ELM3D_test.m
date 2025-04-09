% Contributors: Mohammed Al Alawi
% Course number: ASEN 4028
% File name: main_ELM3D_test
% Created: 2/19/2025

% housekeeping
clear; clc; close all

% calling constants
constants;

% the symbolic variables of the system
syms q1(t) q2(t) q3(t) d1(t) d2(t) d3(t) tau1(t) tau2(t) tau3(t)

% deck movement
dx = 0.2*t; %[m]
dy = 0*t; %[m]
period = 7.5; %[s]
dz = 0.5*sin((2*pi/period)*t); %[m]

d =[dx;dy;dz];
d_dot = diff(d,'t');

% the deck rotations
angle1 = 0*t; % deck rotation about its x-axis [rad]
angle2 = 0*t; %10*(pi/180)*sin((2*pi/(period*2))*t); %deck rotation about its y-axis [rad]
angle3 = 0*t; %deck rotation about its z-axis [rad]

theta_D = [angle1;angle2;angle3];

% Initial Conditions (IC)
q1_init = -0.7248;
Dq1_init = 0;
q2_init = 0.2228;
Dq2_init = 0;
q3_init = (1/18)*pi;
Dq3_init = 0;

% moments of inertia
I1 = [platform.Ixx1 platform.Ixy1 platform.Ixz1;platform.Iyx1 platform.Iyy1 platform.Iyz1; platform.Izx1 platform.Izy1 platform.Izz1]; %link 1
I2 = [platform.Ixx2 platform.Ixy2 platform.Ixz2;platform.Iyx2 platform.Iyy2 platform.Iyz2; platform.Izx2 platform.Izy2 platform.Izz2]; %link 2

Irotor1 = [0 0 0;0 platform.Irotor1 0;0 0 0];
Irotor2 = [0 0 0;0 platform.Irotor2 0;0 0 0];
Irotor3 = [0 0 0;0 0 0;0 0 platform.Irotor3];

Itot1 = I1+(platform.N1^2)*Irotor1+(platform.N3^2)*Irotor3;
Itot2 = I2+(platform.N2^2)*Irotor2;

% rotation matrices
rot1 = [cos(d3) -sin(d3) 0;sin(d3) cos(d3) 0;0 0 1]; %rotation about deck's z axis
rot2 = [cos(d2) 0 sin(d2);0 1 0;-sin(d2) 0 cos(d2)]; %rotation about deck's y axis
rot3 = [1 0 0;0 cos(d1) -sin(d1);0 sin(d1) cos(d1)]; %rotation about deck's x axis

rot_1_0 = [cos(q1)*cos(q3) -sin(q3) sin(q1)*cos(q3);cos(q1)*sin(q3) cos(q3) sin(q1)*sin(q3);-sin(q1) 0 cos(q1)]; %rotation of link 1 with respect to the base
rot_2_0 = [cos(q1+q2)*cos(q3) -sin(q3) sin(q1+q2)*cos(q3);cos(q1+q2)*sin(q3) cos(q3) sin(q1+q2)*sin(q3);-sin(q1+q2) 0 cos(q1+q2)]; % rotation of link2 with respect to the base

rot_0_I = rot1*rot2*rot3; % deck's 321 rotation matrix with respect to inertial frame

R1 = rot_0_I*rot_1_0; % link 1 rotation matrix with respect to the inertial frame
R2 = rot_0_I*rot_2_0; % link 2 rotation matrix with respect to the inertial frame

% controls
%-------------------------------------------------------------------------%
% JOINT-BASED CONTROL (Relative Position Control)
% reference angles
ref_q1 = -(2/6)*pi;
ref_q2 = (2/9)*pi;
ref_q3 = pi/3;

% gains of motor 1
Kp1 = 20;
Ki1 = 10;
Kd1 = 6;

% gains of motor 2
Kp2 = 6;
Ki2 = 2;
Kd2= 6;

% gains of motor 3
Kp3 = 10;
Ki3 = 6;
Kd3 = 10;

tau1(t) =  Kp1*(ref_q1-q1(t))+Ki1*int(ref_q1-q1(t))+Kd1*diff(ref_q1-q1(t),'t');
tau2(t) =  Kp2*(ref_q2-q2(t))+Ki2*int(ref_q2-q2(t))+Kd2*diff(ref_q2-q2(t),'t');
tau3(t) =  Kp3*(ref_q3-q3(t))+Ki3*int(ref_q3-q3(t))+Kd3*diff(ref_q3-q3(t),'t');

% END of JOINT-BASED CONTROL
%-------------------------------------------------------------------------%
%-------------------------------------------------------------------------%
% CARTESIAN-BASED CONTROL (Relative Position & Inertial Stability Controls)
% gains of inertial stability
Ka = 600;
Kv = 750;

% forward kinematics
r_0 = [platform.l2*cos(q1+q2)*cos(q3)+platform.l1*cos(q1)*cos(q3);platform.l2*cos(q1+q2)*sin(q3)+platform.l1*cos(q1)*sin(q3);-platform.l2*sin(q1+q2)-platform.l1*sin(q1)];
r_I = rot_0_I*r_0;

% inertial acceleration
p_ddot = diff(d+r_I,'t',2);

% Jocabian matrix
J1_1 = (cos(d1)*sin(d3)-cos(d3)*sin(d1)*sin(d2))*(platform.l1*sin(q1)*sin(q3)+platform.l2*sin(q3)*sin(q1+q2)) ...
    -(sin(d1)*sin(d3)+cos(d1)*cos(d3)*sin(d2))*(platform.l1*cos(q1)+platform.l2*cos(q1+q2)) ...
    -cos(d2)*cos(d3)*(platform.l2*cos(q3)*sin(q1+q2)+platform.l1*cos(q3)*sin(q1));

J1_2 = (cos(d1)*sin(d3)-cos(d3)*sin(d1)*sin(d2))*(platform.l2*sin(q3)*sin(q1+q2)) ...
    -(sin(d1)*sin(d3)+cos(d1)*cos(d3)*sin(d2))*(platform.l2*cos(q1+q2)) ...
    -cos(d2)*cos(d3)*(platform.l2*cos(q3)*sin(q1+q2));

J1_3 = (cos(d1)*sin(d3)-cos(d3)*sin(d1)*sin(d2))*(-platform.l2*cos(q3)*cos(q1+q2)-platform.l1*cos(q1)*cos(q3)) ...
    -cos(d2)*cos(d3)*(platform.l2*sin(q3)*cos(q1+q2)+platform.l1*cos(q1)*sin(q3));

J2_1 = -(cos(d1)*cos(d3)+sin(d1)*sin(d2)*sin(d3))*(platform.l1*sin(q1)*sin(q3)+platform.l2*sin(q3)*sin(q1+q2)) ...
    +(cos(d3)*sin(d1)-cos(d1)*sin(d2)*sin(d3))*(platform.l1*cos(q1)+platform.l2*cos(q1+q2)) ...
    -cos(d2)*sin(d3)*(platform.l2*cos(q3)*sin(q1+q2)+platform.l1*cos(q3)*sin(q1));

J2_2 = -(cos(d1)*cos(d3)+sin(d1)*sin(d2)*sin(d3))*(platform.l2*sin(q3)*sin(q1+q2)) ...
    +(cos(d3)*sin(d1)-cos(d1)*sin(d2)*sin(d3))*(platform.l2*cos(q1+q2)) ...
    -cos(d2)*sin(d3)*(platform.l2*cos(q3)*sin(q1+q2));

J2_3 = -(cos(d1)*cos(d3)+sin(d1)*sin(d2)*sin(d3))*(-platform.l2*cos(q3)*cos(q1+q2)-platform.l1*cos(q1)*cos(q3)) ...
    -cos(d2)*sin(d3)*(platform.l2*sin(q3)*cos(q1+q2)+platform.l1*cos(q1)*sin(q3));

J3_1 = sin(d2)*(platform.l2*cos(q3)*sin(q1+q2)+platform.l1*cos(q3)*sin(q1)) ...
    -cos(d2)*sin(d1)*(platform.l1*sin(q1)*sin(q3)+platform.l2*sin(q3)*sin(q1+q2)) ...
    -cos(d1)*cos(d2)*(platform.l1*cos(q1)+platform.l2*cos(q1+q2));

J3_2 = sin(d2)*(platform.l2*cos(q3)*sin(q1+q2)) ...
    -cos(d2)*sin(d1)*(platform.l2*sin(q3)*sin(q1+q2)) ...
    -cos(d1)*cos(d2)*(platform.l2*cos(q1+q2));

J3_3 = sin(d2)*(platform.l2*sin(q3)*cos(q1+q2)+platform.l1*cos(q1)*sin(q3)) ...
    -cos(d2)*sin(d1)*(-platform.l2*cos(q3)*cos(q1+q2)-platform.l1*cos(q1)*cos(q3));

J = formula([J1_1, J1_2, J1_3;J2_1, J2_2, J2_3;J3_1, J3_2, J3_3]);

F_control = -Ka*p_ddot+Kv*int(-p_ddot);

tau = formula((J.')*F_control);

tau1(t) = tau(1);
tau2(t) = tau(2);
tau3(t) = tau(3);

% End of CARTESIAN-BASED CONTROL
%-------------------------------------------------------------------------%

% jocabian matrices
% Jw1 & Jw2
Jw1 = [0 0 0;1 0 0;0 0 1];
Jw2 = [0 0 0;1 1 0;0 0 1];

% Jv1 & Jv1s
Jv11_1 = (cos(d1)*sin(d3)-cos(d3)*sin(d1)*sin(d2))*(platform.r1*sin(q1)*sin(q3)) ...
    -(sin(d1)*sin(d3)+cos(d1)*cos(d3)*sin(d2))*(platform.r1*cos(q1)) ...
    -cos(d2)*cos(d3)*(platform.r1*cos(q3)*sin(q1));

Jv11_2 = 0;

Jv11_3 = (cos(d1)*sin(d3)-cos(d3)*sin(d1)*sin(d2))*(-platform.r1*cos(q1)*cos(q3)) ...
    -cos(d2)*cos(d3)*(platform.r1*cos(q1)*sin(q3));

Jv12_1 = -(cos(d1)*cos(d3)+sin(d1)*sin(d2)*sin(d3))*(platform.r1*sin(q1)*sin(q3)) ...
    +(cos(d3)*sin(d1)-cos(d1)*sin(d2)*sin(d3))*(platform.r1*cos(q1)) ...
    -cos(d2)*sin(d3)*(platform.r1*cos(q3)*sin(q1));

Jv12_2 = 0;

Jv12_3 = -(cos(d1)*cos(d3)+sin(d1)*sin(d2)*sin(d3))*(-platform.r1*cos(q1)*cos(q3)) ...
    -cos(d2)*sin(d3)*(platform.r1*cos(q1)*sin(q3));

Jv13_1 = sin(d2)*(platform.r1*cos(q3)*sin(q1)) ...
    -cos(d2)*sin(d1)*(platform.r1*sin(q1)*sin(q3)) ...
    -cos(d1)*cos(d2)*(platform.r1*cos(q1));

Jv13_2 = 0;

Jv13_3 = sin(d2)*(platform.r1*cos(q1)*sin(q3)) ...
    -cos(d2)*sin(d1)*(-platform.r1*cos(q1)*cos(q3));

Jv1s1_1 = -(platform.r1*sin(q1))*(cos(d1)*sin(d3)-cos(d3)*sin(d1)*sin(d2)) ...
    -sin(q3)*(platform.r1*cos(q1))*(sin(d1)*sin(d3)+cos(d1)*cos(d3)*sin(d2));

Jv1s1_2 = -(platform.r1*sin(q1))*(cos(d1)*cos(d2)*cos(d3)) ...
    -sin(q3)*(platform.r1*cos(q1))*(cos(d2)*cos(d3)*sin(d1)) ...
    -cos(q3)*cos(d3)*sin(d2)*(platform.r1*cos(q1));

Jv1s1_3 = -(platform.r1*sin(q1))*(cos(d3)*sin(d1)-cos(d1)*sin(d2)*sin(d3)) ...
    -sin(q3)*(platform.r1*cos(q1))*(-cos(d1)*cos(d2)-sin(d1)*sin(d2)*sin(d3)) ...
    -cos(q3)*cos(d2)*sin(d3)*(platform.r1*cos(q1));

Jv1s2_1 = -(platform.r1*sin(q1))*(-cos(d1)*cos(d3)-sin(d1)*sin(d2)*sin(d3)) ...
    +sin(q3)*(platform.r1*cos(q1))*(cos(d1)*sin(d2)*sin(d3)-cos(d3)*sin(d1));

Jv1s2_2 = -(platform.r1*sin(q1))*(cos(d1)*cos(d2)*sin(d3)) ...
    +sin(q3)*(platform.r1*cos(q1))*(cos(d2)*sin(d1)*sin(d3)) ...
    -cos(q3)*sin(d2)*sin(d3)*(platform.r1*cos(q1));

Jv1s2_3 = -(platform.r1*sin(q1))*(sin(d1)*sin(d3)+cos(d1)*cos(d3)*sin(d2)) ...
    +sin(q3)*(platform.r1*cos(q1))*(-cos(d1)*sin(d3)+cos(d3)*sin(d1)*sin(d2)) ...
    +cos(q3)*cos(d2)*cos(d3)*(platform.r1*cos(q1));

Jv1s3_1 = cos(d2)*sin(d1)*(platform.r1*sin(q1)) ...
    +cos(d1)*cos(d2)*sin(q3)*(platform.r1*cos(q1));

Jv1s3_2 = -cos(q3)*cos(d2)*(platform.r1*cos(q1)) ...
    +cos(d1)*sin(d2)*(platform.r1*sin(q1)) ...
    -sin(q3)*sin(d1)*sin(d2)*(platform.r1*cos(q1));

Jv1s3_3 = 0;

Jv1 = [Jv11_1, Jv11_2, Jv11_3;Jv12_1, Jv12_2, Jv12_3;Jv13_1, Jv13_2, Jv13_3];
Jv1s = [Jv1s1_1, Jv1s1_2, Jv1s1_3;Jv1s2_1, Jv1s2_2, Jv1s2_3;Jv1s3_1, Jv1s3_2, Jv1s3_3];

% Jv2 & Jv2s
Jv21_1 = (cos(d1)*sin(d3)-cos(d3)*sin(d1)*sin(d2))*(platform.l1*sin(q1)*sin(q3)+platform.r2*sin(q3)*sin(q1+q2)) ...
    -(sin(d1)*sin(d3)+cos(d1)*cos(d3)*sin(d2))*(platform.l1*cos(q1)+platform.r2*cos(q1+q2)) ...
    -cos(d2)*cos(d3)*(platform.r2*cos(q3)*sin(q1+q2)+platform.l1*cos(q3)*sin(q1));

Jv21_2 = (cos(d1)*sin(d3)-cos(d3)*sin(d1)*sin(d2))*(platform.r2*sin(q3)*sin(q1+q2)) ...
    -(sin(d1)*sin(d3)+cos(d1)*cos(d3)*sin(d2))*(platform.r2*cos(q1+q2)) ...
    -cos(d2)*cos(d3)*(platform.r2*cos(q3)*sin(q1+q2));

Jv21_3 = (cos(d1)*sin(d3)-cos(d3)*sin(d1)*sin(d2))*(-platform.r2*cos(q3)*cos(q1+q2)-platform.l1*cos(q1)*cos(q3)) ...
    -cos(d2)*cos(d3)*(platform.r2*sin(q3)*cos(q1+q2)+platform.l1*cos(q1)*sin(q3));

Jv22_1 = -(cos(d1)*cos(d3)+sin(d1)*sin(d2)*sin(d3))*(platform.l1*sin(q1)*sin(q3)+platform.r2*sin(q3)*sin(q1+q2)) ...
    +(cos(d3)*sin(d1)-cos(d1)*sin(d2)*sin(d3))*(platform.l1*cos(q1)+platform.r2*cos(q1+q2)) ...
    -cos(d2)*sin(d3)*(platform.r2*cos(q3)*sin(q1+q2)+platform.l1*cos(q3)*sin(q1));

Jv22_2 = -(cos(d1)*cos(d3)+sin(d1)*sin(d2)*sin(d3))*(platform.r2*sin(q3)*sin(q1+q2)) ...
    +(cos(d3)*sin(d1)-cos(d1)*sin(d2)*sin(d3))*(platform.r2*cos(q1+q2)) ...
    -cos(d2)*sin(d3)*(platform.r2*cos(q3)*sin(q1+q2));

Jv22_3 = -(cos(d1)*cos(d3)+sin(d1)*sin(d2)*sin(d3))*(-platform.r2*cos(q3)*cos(q1+q2)-platform.l1*cos(q1)*cos(q3)) ...
    -cos(d2)*sin(d3)*(platform.r2*sin(q3)*cos(q1+q2)+platform.l1*cos(q1)*sin(q3));

Jv23_1 = sin(d2)*(platform.r2*cos(q3)*sin(q1+q2)+platform.l1*cos(q3)*sin(q1)) ...
    -cos(d2)*sin(d1)*(platform.l1*sin(q1)*sin(q3)+platform.r2*sin(q3)*sin(q1+q2)) ...
    -cos(d1)*cos(d2)*(platform.l1*cos(q1)+platform.r2*cos(q1+q2));

Jv23_2 = sin(d2)*(platform.r2*cos(q3)*sin(q1+q2)) ...
    -cos(d2)*sin(d1)*(platform.r2*sin(q3)*sin(q1+q2)) ...
    -cos(d1)*cos(d2)*(platform.r2*cos(q1+q2));

Jv23_3 = sin(d2)*(platform.r2*sin(q3)*cos(q1+q2)+platform.l1*cos(q1)*sin(q3)) ...
    -cos(d2)*sin(d1)*(-platform.r2*cos(q3)*cos(q1+q2)-platform.l1*cos(q1)*cos(q3));

Jv2s1_1 = -(platform.l1*sin(q1)+platform.r2*sin(q1+q2))*(cos(d1)*sin(d3)-cos(d3)*sin(d1)*sin(d2)) ...
    -sin(q3)*(platform.l1*cos(q1)+platform.r2*cos(q1+q2))*(sin(d1)*sin(d3)+cos(d1)*cos(d3)*sin(d2));

Jv2s1_2 = -(platform.l1*sin(q1)+platform.r2*sin(q1+q2))*(cos(d1)*cos(d2)*cos(d3)) ...
    -sin(q3)*(platform.l1*cos(q1)+platform.r2*cos(q1+q2))*(cos(d2)*cos(d3)*sin(d1)) ...
    -cos(q3)*cos(d3)*sin(d2)*(platform.l1*cos(q1)+platform.r2*cos(q1+q2));

Jv2s1_3 = -(platform.l1*sin(q1)+platform.r2*sin(q1+q2))*(cos(d3)*sin(d1)-cos(d1)*sin(d2)*sin(d3)) ...
    -sin(q3)*(platform.l1*cos(q1)+platform.r2*cos(q1+q2))*(-cos(d1)*cos(d2)-sin(d1)*sin(d2)*sin(d3)) ...
    -cos(q3)*cos(d2)*sin(d3)*(platform.l1*cos(q1)+platform.r2*cos(q1+q2));

Jv2s2_1 = -(platform.l1*sin(q1)+platform.r2*sin(q1+q2))*(-cos(d1)*cos(d3)-sin(d1)*sin(d2)*sin(d3)) ...
    +sin(q3)*(platform.l1*cos(q1)+platform.r2*cos(q1+q2))*(cos(d1)*sin(d2)*sin(d3)-cos(d3)*sin(d1));

Jv2s2_2 = -(platform.l1*sin(q1)+platform.r2*sin(q1+q2))*(cos(d1)*cos(d2)*sin(d3)) ...
    +sin(q3)*(platform.l1*cos(q1)+platform.r2*cos(q1+q2))*(cos(d2)*sin(d1)*sin(d3)) ...
    -cos(q3)*sin(d2)*sin(d3)*(platform.l1*cos(q1)+platform.r2*cos(q1+q2));

Jv2s2_3 = -(platform.l1*sin(q1)+platform.r2*sin(q1+q2))*(sin(d1)*sin(d3)+cos(d1)*cos(d3)*sin(d2)) ...
    +sin(q3)*(platform.l1*cos(q1)+platform.r2*cos(q1+q2))*(-cos(d1)*sin(d3)+cos(d3)*sin(d1)*sin(d2)) ...
    +cos(q3)*cos(d2)*cos(d3)*(platform.l1*cos(q1)+platform.r2*cos(q1+q2));

Jv2s3_1 = cos(d2)*sin(d1)*(platform.l1*sin(q1)+platform.r2*sin(q1+q2)) ...
    +cos(d1)*cos(d2)*sin(q3)*(platform.l1*cos(q1)+platform.r2*cos(q1+q2));

Jv2s3_2 = -cos(q3)*cos(d2)*(platform.l1*cos(q1)+platform.r2*cos(q1+q2)) ...
    +cos(d1)*sin(d2)*(platform.l1*sin(q1)+platform.r2*sin(q1+q2)) ...
    -sin(q3)*sin(d1)*sin(d2)*(platform.l1*cos(q1)+platform.r2*cos(q1+q2));

Jv2s3_3 = 0;

Jv2 = [Jv21_1, Jv21_2, Jv21_3;Jv22_1, Jv22_2, Jv22_3;Jv23_1, Jv23_2, Jv23_3];
Jv2s = [Jv2s1_1, Jv2s1_2, Jv2s1_3;Jv2s2_1, Jv2s2_2, Jv2s2_3;Jv2s3_1, Jv2s3_2, Jv2s3_3];

% constructing the inertia matrices M, M_star, M_dstar
M_matrix = platform.m1*(Jv1.')*Jv1+(Jw1.')*R1*Itot1*(R1.')*Jw1+ platform.m2*(Jv2.')*Jv2+(Jw2.')*R2*Itot2*(R2.')*Jw2;
M = formula(M_matrix);

M_star_matrix =  platform.m1*(Jv1s.')*Jv1s+R1*Itot1*(R1.')+ platform.m2*(Jv2s.')*Jv2s+R2*Itot2*(R2.');
M_star = formula(M_star_matrix);

M_dstar_matrix =  platform.m1*(Jv1.')*Jv1s+(Jw1.')*R1*Itot1*(R1.')+ platform.m2*(Jv2.')*Jv2s+(Jw2.')*R2*Itot2*(R2.');
M_dstar = formula(M_dstar_matrix);

% decomposing the inertia matrices
m11 = M(1,1);
m12 = M(1,2);
m13 = M(1,3);
m21 = M(2,1);
m22 = M(2,2);
m23 = M(2,3);
m31 = M(3,1);
m32 = M(3,2);
m33 = M(3,3);

ms11 = M_star(1,1);
ms12 = M_star(1,2);
ms13 = M_star(1,3);
ms21 = M_star(2,1);
ms22 = M_star(2,2);
ms23 = M_star(2,3);
ms31 = M_star(3,1);
ms32 = M_star(3,2);
ms33 = M_star(3,3);

mds11 = M_dstar(1,1);
mds12 = M_dstar(1,2);
mds13 = M_dstar(1,3);
mds21 = M_dstar(2,1);
mds22 = M_dstar(2,2);
mds23 = M_dstar(2,3);
mds31 = M_dstar(3,1);
mds32 = M_dstar(3,2);
mds33 = M_dstar(3,3);

% constructing the C, C_star, C_dstar matrices
% C matrix
c111 = formula(0.5*(diff(m11,q1)+diff(m11,q1)-diff(m11,q1)));
c112 = formula(0.5*(diff(m21,q1)+diff(m21,q1)-diff(m11,q2)));
c113 = formula(0.5*(diff(m31,q1)+diff(m31,q1)-diff(m11,q3)));
c121 = formula(0.5*(diff(m12,q1)+diff(m11,q2)-diff(m12,q1)));
c211 = c121;
c122 = formula(0.5*(diff(m22,q1)+diff(m21,q2)-diff(m12,q2)));
c212 = c122;
c123 = formula(0.5*(diff(m32,q1)+diff(m31,q2)-diff(m12,q3)));
c213 = c123;
c131 = formula(0.5*(diff(m13,q1)+diff(m11,q3)-diff(m13,q1)));
c311 = c131;
c132 = formula(0.5*(diff(m23,q1)+diff(m21,q3)-diff(m13,q2)));
c312 = c132;
c133 = formula(0.5*(diff(m33,q1)+diff(m31,q3)-diff(m13,q3)));
c313 = c133;
c221 = formula(0.5*(diff(m12,q2)+diff(m12,q2)-diff(m22,q1)));
c222 = formula(0.5*(diff(m22,q2)+diff(m22,q2)-diff(m22,q2)));
c223 = formula(0.5*(diff(m32,q2)+diff(m32,q2)-diff(m22,q3)));
c231 = formula(0.5*(diff(m13,q2)+diff(m12,q3)-diff(m23,q1)));
c321 = c231;
c232 = formula(0.5*(diff(m23,q2)+diff(m22,q3)-diff(m23,q2)));
c322 = c232;
c233 = formula(0.5*(diff(m33,q2)+diff(m32,q3)-diff(m23,q3)));
c323 = c233;
c331 = formula(0.5*(diff(m13,q3)+diff(m13,q3)-diff(m33,q1)));
c332 = formula(0.5*(diff(m23,q3)+diff(m23,q3)-diff(m33,q2)));
c333 = formula(0.5*(diff(m33,q3)+diff(m33,q3)-diff(m33,q3)));

q1_dot = diff(q1(t),'t');
q2_dot = diff(q2(t),'t');
q3_dot = diff(q3(t),'t');

C = [c111*q1_dot+c121*q2_dot+c131*q3_dot c211*q1_dot+c221*q2_dot+c231*q3_dot c311*q1_dot+c321*q2_dot+c331*q3_dot;
    c112*q1_dot+c122*q2_dot+c132*q3_dot c212*q1_dot+c222*q2_dot+c232*q3_dot c312*q1_dot+c322*q2_dot+c332*q3_dot;
    c113*q1_dot+c123*q2_dot+c133*q3_dot c213*q1_dot+c223*q2_dot+c233*q3_dot c313*q1_dot+c323*q2_dot+c333*q3_dot];

% C_star matrix
cs111 = formula(0.5*(diff(mds11,d1)+diff(mds11,d1)-diff(ms11,q1)));
cs112 = formula(0.5*(diff(mds21,d1)+diff(mds21,d1)-diff(ms11,q2)));
cs113 = formula(0.5*(diff(mds31,d1)+diff(mds31,d1)-diff(ms11,q3)));
cs121 = formula(0.5*(diff(mds12,d1)+diff(mds11,d2)-diff(ms12,q1)));
cs211 = cs121;
cs122 = formula(0.5*(diff(mds22,d1)+diff(mds21,d2)-diff(ms12,q2)));
cs212 = cs122;
cs123 = formula(0.5*(diff(mds32,d1)+diff(mds31,d2)-diff(ms12,q3)));
cs213 = cs123;
cs131 = formula(0.5*(diff(mds13,d1)+diff(mds11,d3)-diff(ms13,q1)));
cs311 = cs131;
cs132 = formula(0.5*(diff(mds23,d1)+diff(mds21,d3)-diff(ms13,q2)));
cs312 = cs132;
cs133 = formula(0.5*(diff(mds33,d1)+diff(mds31,d3)-diff(ms13,q3)));
cs313 = cs133;
cs221 = formula(0.5*(diff(mds12,d2)+diff(mds12,d2)-diff(ms22,q1)));
cs222 = formula(0.5*(diff(mds22,d2)+diff(mds22,d2)-diff(ms22,q2)));
cs223 = formula(0.5*(diff(mds32,d2)+diff(mds32,d2)-diff(ms22,q3)));
cs231 = formula(0.5*(diff(mds13,d2)+diff(mds12,d3)-diff(ms23,q1)));
cs321 = cs231;
cs232 = formula(0.5*(diff(mds23,d2)+diff(mds22,d3)-diff(ms23,q2)));
cs322 = cs232;
cs233 = formula(0.5*(diff(mds33,d2)+diff(mds32,d3)-diff(ms23,q3)));
cs323 = cs233;
cs331 = formula(0.5*(diff(mds13,d3)+diff(mds13,d3)-diff(ms33,q1)));
cs332 = formula(0.5*(diff(mds23,d3)+diff(mds23,d3)-diff(ms33,q2)));
cs333 = formula(0.5*(diff(mds33,d3)+diff(mds33,d3)-diff(ms33,q3)));

d1_dot = diff(d1(t),'t');
d2_dot = diff(d2(t),'t');
d3_dot = diff(d3(t),'t');

C_star = [cs111*d1_dot+cs121*d2_dot+cs131*d3_dot cs211*d1_dot+cs221*d2_dot+cs231*d3_dot cs311*d1_dot+cs321*d2_dot+cs331*d3_dot;
    cs112*d1_dot+cs122*d2_dot+cs132*d3_dot cs212*d1_dot+cs222*d2_dot+cs232*d3_dot cs312*d1_dot+cs322*d2_dot+cs332*d3_dot;
    cs113*d1_dot+cs123*d2_dot+cs133*d3_dot cs213*d1_dot+cs223*d2_dot+cs233*d3_dot cs313*d1_dot+cs323*d2_dot+cs333*d3_dot];

% C_dstar matrix
cds111 = formula((diff(mds11,q1)+diff(m11,d1)-diff(mds11,q1)));
cds112 = formula((diff(mds21,q1)+diff(m21,d1)-diff(mds11,q2)));
cds113 = formula((diff(mds31,q1)+diff(m31,d1)-diff(mds11,q3)));
cds121 = formula((diff(mds12,q1)+diff(m11,d2)-diff(mds12,q1)));
cds211 = formula((diff(mds11,q2)+diff(m12,d1)-diff(mds21,q1)));
cds122 = formula((diff(mds22,q1)+diff(m21,d2)-diff(mds12,q2)));
cds212 = formula((diff(mds21,q2)+diff(m22,d1)-diff(mds21,q2)));
cds123 = formula((diff(mds32,q1)+diff(m31,q2)-diff(mds12,q3)));
cds213 = formula((diff(mds31,q2)+diff(m32,q1)-diff(mds21,q3)));
cds131 = formula((diff(mds13,q1)+diff(m11,d3)-diff(mds13,q1)));
cds311 = formula((diff(mds11,q3)+diff(m13,d1)-diff(mds31,q1)));
cds132 = formula((diff(mds23,q1)+diff(m21,d3)-diff(mds13,q2)));
cds312 = formula((diff(mds21,q3)+diff(m23,d1)-diff(mds31,q2)));
cds133 = formula((diff(mds33,q1)+diff(m31,d3)-diff(mds13,q3)));
cds313 = formula((diff(mds31,q3)+diff(m33,d1)-diff(mds31,q3)));
cds221 = formula((diff(mds12,q2)+diff(m12,d2)-diff(mds22,q1)));
cds222 = formula((diff(mds22,q2)+diff(m22,d2)-diff(mds22,q2)));
cds223 = formula((diff(mds32,q2)+diff(m32,d2)-diff(mds22,q3)));
cds231 = formula((diff(mds13,q2)+diff(m12,d3)-diff(mds23,q1)));
cds321 = formula((diff(mds12,q3)+diff(m13,d2)-diff(mds32,q1)));
cds232 = formula((diff(mds23,q2)+diff(m22,d3)-diff(mds23,q2)));
cds322 = formula((diff(mds22,q3)+diff(m23,d2)-diff(mds32,q2)));
cds233 = formula((diff(mds33,q2)+diff(m32,d3)-diff(mds23,q3)));
cds323 = formula((diff(mds32,q3)+diff(m33,d2)-diff(mds32,q3)));
cds331 = formula((diff(mds13,q3)+diff(m13,d3)-diff(mds33,q1)));
cds332 = formula((diff(mds23,q3)+diff(m23,d3)-diff(mds33,q2)));
cds333 = formula((diff(mds33,q3)+diff(m33,d3)-diff(mds33,q3)));

C_dstar = [cds111*d1_dot+cds121*d2_dot+cds131*d3_dot cds211*d1_dot+cds221*d2_dot+cds231*d3_dot cds311*d1_dot+cds321*d2_dot+cds331*d3_dot;
    cds112*d1_dot+cds122*d2_dot+cds132*d3_dot cds212*d1_dot+cds222*d2_dot+cds232*d3_dot cds312*d1_dot+cds322*d2_dot+cds332*d3_dot;
    cds113*d1_dot+cds123*d2_dot+cds133*d3_dot cds213*d1_dot+cds223*d2_dot+cds233*d3_dot cds313*d1_dot+cds323*d2_dot+cds333*d3_dot];

% constructing the G vector
% potential energy of the system
z1_I = cos(d2(t))*sin(d1(t))*(platform.r1*cos(q1(t))*sin(q3(t))) - sin(d2(t))*(platform.r1*cos(q1(t))*cos(q3(t))) - cos(d1(t))*cos(d2(t))*(platform.r1*sin(q1(t)));
z2_I = cos(d2(t))*sin(d1(t))*(platform.l1*cos(q1(t))*sin(q3(t)) + platform.r2*sin(q3(t))*cos(q1(t) + q2(t))) - sin(d2(t))*(platform.l1*cos(q1(t))*cos(q3(t)) + platform.r2*cos(q3(t))*cos(q1(t) + q2(t))) - cos(d1(t))*cos(d2(t))*(platform.l1*sin(q1(t)) + platform.r2*sin(q1(t) + q2(t)));

P = platform.m1*platform.g*z1_I+platform.m2*platform.g*z2_I;

g1 = diff(P,q1);
g2 = diff(P,q2);
g3 = diff(P,q3);

G = formula([g1;g2;g3]);

% constructing the j and j_star vectors
j = formula(platform.m1*(Jv1.')*d_dot+platform.m2*(Jv2.')*d_dot);
j_star =  formula(platform.m1*(Jv1s.')*d_dot+platform.m2*(Jv2s.')*d_dot);

% constructing the j_q and j_q_star matrices
j_q = formula([diff(j(1),q1), diff(j(2),q1), diff(j(3),q1);diff(j(1),q2), diff(j(2),q2), diff(j(3),q2);diff(j(1),q3), diff(j(2),q3), diff(j(3),q3)]);
j_q_star = formula([diff(j_star(1),q1), diff(j_star(2),q1), diff(j_star(3),q1);diff(j_star(1),q2), diff(j_star(2),q2), diff(j_star(3),q2);diff(j_star(1),q3), diff(j_star(2),q3), diff(j_star(3),q3)]);

% important matrices
MU = [platform.mu1 0 0;0 platform.mu2 0;0 0 platform.mu3]; %static friction matrix
B = [platform.b1 0 0;0 platform.b2 0;0 0 platform.b3]; %kinetic friction matrix
N = [platform.N1*platform.eta1 0 0;0 platform.N2*platform.eta2 0;0 0 platform.N3*platform.eta3]; %gearhead ratio/efficiency matrix

% constructing EOM
EOM = diff(j,'t')+M_dstar*diff([d1(t);d2(t);d3(t)],'t',2)+(C_star-j_q_star)*diff([d1(t);d2(t);d3(t)],'t')+M*diff([q1(t);q2(t);q3(t)],'t',2)+(C+C_dstar+B-j_q)*diff([q1(t);q2(t);q3(t)],'t')+MU*sign(diff([q1(t);q2(t);q3(t)],'t'))+G == N*[tau1(t);tau2(t);tau3(t)];

% substituting deck rotations in EOM
EOM = subs(EOM,d1,angle1);
EOM = subs(EOM,d2,angle2);
EOM = subs(EOM,d3,angle3);

% making the system compatible to be solved by any MATLAB solver
[V,S] = odeToVectorField(EOM(1,1),EOM(2,1),EOM(3,1));
MFun = matlabFunction(V,'vars',{'t','Y'});

%solving the system
freq = 160; %Hz
time_interval = [0 30]; %seconds

initial_conditions = [q2_init Dq2_init q1_init Dq1_init q3_init Dq3_init]; %[q2 Dq2 q1 Dq1 q3 Dq3]

[sol.x, sol.y]= rk4_solver(MFun,time_interval,initial_conditions,1/freq);
sol.y = sol.y';

% plot q3, q3_dot, q1, q1_dot, q2, q2_dot over time
plot2_test(sol.x,sol.y)

%animate
animate3D_test(sol.x,sol.y,d, theta_D)