% Contributors: Mohammed Al Alawi, Kirin Kawamoto
% Course number: ASEN 4028
% File name: main
% Created: 4/3/2025

% housekeeping
clear; clc; close all

% calling constants
constants;

% controls
%-------------------------------------------------------------------------%
% JOINT-BASED & CARTESIAN-BASED CONTROL
% reference angle
ref_q1 = -pi/4;

% reference relative position (inertial)
r_I_ref = [0.5;0;0.2];

% forward kinematics
r_I = [platform.l1*cos(q1(t)+theta2(t));0;-platform.l1*sin(q1(t)+theta2(t))];

% the jocabain
J = [-platform.l1*sin(q1(t)+theta2(t));
     0;
     -platform.l1*cos(q1(t)+theta2(t))];

% inertial acceleration
p = d+r_I;
p_ddot = diff(p,'t',2);

% gains of inertial stability control (cartesian)
Ka = 300*1;
Kv = 375*1;

% gains of relative position control (cartesian)
Kp = 15*0;
Ki = 2*0;
Kd = 10*0;

% gains of motor 1 (joint)
Kp1 = 4*0;
Ki1 = 10*0;
Kd1 = 15*0;

F_control = Kp*(r_I_ref-r_I)+Ki*int(r_I_ref-r_I)+Kd*diff(r_I_ref-r_I)-Ka*p_ddot+Kv*int(-p_ddot);

tau = J.'*F_control ;