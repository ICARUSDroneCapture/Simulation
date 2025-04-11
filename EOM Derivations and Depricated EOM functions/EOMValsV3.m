% Contributors: Mohammed Al Alawi
% Course number: ASEN 4028
% File name: main_ELM1D
% Created: 2/28/2025


% housekeeping
clear; clc; close all

% moments of inertia matrices
syms Ixx Iyy Izz Ixy Ixz Iyz

Itot1 = [Ixx Ixy Ixz;
         Ixy Iyy Iyz;
         Ixz Iyz Izz];

% defining the symbolic math we use for angles
syms q1(t) theta1(t) theta2(t) theta3(t) tau1(t) 

% the base movement in the inertial frame 
dx = 0*t; %[m] (DO NOT CHANGE)
dy = 0*t; %[m] (DO NOT CHANGE)
dz = 0*t; %[m] (DO NOT CHANGE)

d =[dx;dy;dz];
d_dot = diff(d,'t');

% the deck rotations
period = 7.5*4; %s
angle1 = 0*t; % deck rotation about its x-axis [rad] (DO NOT CHANGE)
angle2 = 90*(pi/180)*(sin((2*pi/(period))*t))^2; %deck rotation about its y-axis [rad]
angle3 = 0*t; %deck rotation about its z-axis [rad] (DO NOT CHANGE)

% Initial Conditions
q1_init = -0.1;
Dq1_init = 0;

% controls
%-------------------------------------------------------------------------%
% JOINT-BASED & CARTESIAN-BASED CONTROL
% reference angle
% ref_q1 = -pi/4;
% 
syms l1
% % forward kinematics
r_I = [l1*cos(q1(t)+theta2(t));0;-l1*sin(q1(t)+theta2(t))];
% 
% % the jocabain
% J = [-l1*sin(q1(t));0;-l1*cos(q1(t))];
% 
% % inertial acceleration
p = d+r_I;
p_ddot = diff(p,'t',2);
p_ddot = simplify(p_ddot)
% 
% R_I_B = [cos(theta2(t)) 0 -sin(theta2(t));
%         0 1 0;
%         sin(theta2(t)) 0 cos(theta2(t))];
% % p_ddot_B = R_I_B*p_ddot;
% 
% syms ka tau_c
% % gains of inertial stability control (cartesian)
% Ka = [0 0 0;
%       0 0 0; 
%       0 0 ka];
% Kv = [0 0 0; 
%       0 0 0; 
%       0 0 375*1];
% 
% % gains of motor 1 (joint)
% Kp1 = 4*0;
% Ki1 = 10*0;
% Kd1 = 15*0;
% 
% F_control = - R_I_B*(Ka*p_ddot + Kv*int(p_ddot));
% 
% F_control = - R_I_B*Ka*p_ddot;
% 
% tau = J.'*F_control ;
% 
% tau1(t) =  tau_c+tau(1);
%-------------------------------------------------------------------------%

% defining the rotation matrices
R1 = [cos(q1(t)+theta2(t)) 0 sin(q1(t)+theta2(t));0 1 0;-sin(q1(t)+theta2(t)) 0 cos(q1(t)+theta2(t))];

syms r1
% defining the jocabian matrices
Jw1 = [0;1;0];
Jv1 = [-r1*sin(q1(t)+theta2(t));0;-r1*cos(q1(t)+theta2(t))];
Jv1_star = [0 -r1*sin(q1(t)+theta2(t)) 0;0 0 0;0 -r1*cos(q1(t)+theta2(t)) 0];

syms m1
% constructing the inertia matrix M
M = (m1*(Jv1.')*Jv1+(Jw1.')*R1*Itot1*(R1.')*Jw1);

% constructing the inertia matrix M_star
M_star = m1*(Jv1_star.')*Jv1_star+R1*Itot1*(R1.');

% constructing the inertia matrix M_dstar
M_dstar = m1*(Jv1.')*Jv1_star+(Jw1.')*R1*Itot1*(R1.');

% constructing the j vector
j = m1*(Jv1.')*d_dot;

%constructing the j_q matrix
j_q = diff(j(1),q1);
j_q = j_q(t);

% constructing the j_star vector
j_star = m1*(Jv1_star.')*d_dot;

%constructing the j_star_q matrix
j_star_q = [diff(j_star(1),q1) diff(j_star(2),q1) diff(j_star(3),q1)];
j_star_q = j_star_q(t);

% constructing the centripetal/corilios matrix C
m11 = M(1,1);

c111 = 0.5*(diff(m11,q1)+diff(m11,q1)-diff(m11,q1));
c111 = c111(t);

C = c111*diff(q1(t),'t');

% constructing the centripetal/corilios matrix C_star
m_star11 = M_star(1,1);
m_star12 = M_star(1,2);
m_star13 = M_star(1,3);
m_star21 = M_star(2,1);
m_star22 = M_star(2,2);
m_star23 = M_star(2,3);
m_star31 = M_star(3,1);
m_star32 = M_star(3,2);
m_star33 = M_star(3,3);


m_dstar11 = M_dstar(1,1);
m_dstar12 = M_dstar(1,2);
m_dstar13 = M_dstar(1,3);

c_star111 = 0.5*(diff(m_dstar11,theta1)+diff(m_dstar11,theta1)-diff(m_star11,q1));
c_star111 = c_star111(t);

c_star211 = 0.5*(diff(m_dstar11,theta2)+diff(m_dstar12,theta1)-diff(m_star21,q1));
c_star211 = c_star211(t);
c_star121 = c_star211;

c_star311 = 0.5*(diff(m_dstar11,theta3)+diff(m_dstar13,theta1)-diff(m_star31,q1));
c_star311 = c_star311(t);
c_star131 = c_star311;

c_star231 = 0.5*(diff(m_dstar13,theta2)+diff(m_dstar12,theta3)-diff(m_star23,q1));
c_star231 = c_star231(t);
c_star321 = c_star231;

c_star221 = 0.5*(diff(m_dstar12,theta2)+diff(m_dstar12,theta2)-diff(m_star22,q1));
c_star221 = c_star221(t);

c_star331 = 0.5*(diff(m_dstar13,theta3)+diff(m_dstar13,theta3)-diff(m_star33,q1));
c_star331 = c_star331(t);

C_star =[c_star111*diff(theta1(t),'t')+c_star121*diff(theta2(t),'t')+c_star131*diff(theta3(t),'t'),c_star211*diff(theta1(t),'t')+c_star221*diff(theta2(t),'t')+c_star231*diff(theta3(t),'t'),c_star311*diff(theta1(t),'t')+c_star321*diff(theta2(t),'t')+c_star331*diff(theta3(t),'t')];

% constructing the centripetal/corilios matrix C_dstar
c_dstar111 = (diff(m_dstar11,q1)+diff(m11,theta1)-diff(m_dstar11,q1));
c_dstar111 = c_dstar111(t);

c_dstar121 = (diff(m_dstar12,q1)+diff(m11,theta2)-diff(m_dstar12,q1));
c_dstar121 = c_dstar121(t);

c_dstar131 = (diff(m_dstar13,q1)+diff(m11,theta3)-diff(m_dstar13,q1));
c_dstar131 = c_dstar131(t);


C_dstar = c_dstar111*diff(theta1(t),'t')+c_dstar121*diff(theta2(t),'t')+c_dstar131*diff(theta3(t),'t');

syms g
% constructing the G matrix (I call it g vector in my notes)
G = g*(-m1*r1*cos(q1(t)+theta2(t)));

% constructing EOM
theta_D = [theta1(t);theta2(t);theta3(t)];

syms MU B N
% important matrices

% EOM
EOM = diff(j,'t')+M*diff(q1(t),'t',2)+M_dstar*diff(theta_D,'t',2)+(C+C_dstar-j_q)*diff(q1(t),'t')+(C_star-j_star_q)*diff(theta_D,'t')+MU*sign(diff(q1(t),'t'))+B*diff(q1(t),'t')+G == N*tau1(t);

% substituting deck rotations in EOM
EOM = subs(EOM,theta1,angle1);
% EOM = subs(EOM,theta2,angle2);
EOM = subs(EOM,theta3,angle3);

EOM = simplify(EOM);
syms QDOT
EOM_V2 = simplify(subs(EOM, diff(q1(t), t, t), QDOT));
[QDOT, params, conditions] = solve(EOM_V2, QDOT,"ReturnConditions",true);
QDOT = simplify(QDOT)