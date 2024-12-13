% Contributors: Mohammed Al Alawi
% Course number: ASEN 4018
% File name: main2
% Created: 10/19/2024


%HouseKeeping
clear; clc; close all

% calling constants
constants;

syms phi1(t) alpha(t) 

%Deck Input
dx = 0.2*t;
peroid = 7.5; %s
dz = 0.4*sin((2*pi/peroid)*t);

%define design choice
r_a1_star_x = -2;
r_a1_star_z = 1;

r_1_star_x = platform.r1*sin(phi1);
r_1_star_z = platform.r1*(cos(phi1)-1);

p1_x = dx+r_a1_star_x+r_1_star_x; %inertial x position of the cg of link 1
p1_z = dz+r_a1_star_z+r_1_star_z; %inertial z position of the cg of link 1

l1_x = platform.l1*sin(phi1);
l1_z = platform.l1*cos(phi1);
l2_x = platform.l2*sin(phi1+alpha);
l2_z = platform.l2*cos(phi1+alpha);
r1_x = platform.r1*sin(phi1);
r1_z = platform.r1*cos(phi1);


p_ee_x = p1_x + l1_x - r1_x + l2_x; %inertial x position of the end effector
p_ee_z = p1_z + l1_z - r1_z + l2_z; %inertial z position of the end effector

p_b_x = p1_x - r1_x; %inertial x position of the base
p_b_z = p1_z - r1_z; %inertial z position of the base

%lets construct the transformation matrix T_B_I
T_B_I =[1 0 p_b_x; 0 1 p_b_z;0 0 1];
%T_I_B = inv(T_B_I)
T_I_B = [1, 0, 2 - t/5;0, 1, - (2*sin((4*pi*t)/15))/5 - 3/5; 0, 0, 1];

%p_ee_base = T_I_B*[p_ee_x;p_ee_z;1]
p_ee_base_x = (4*sin(phi1(t)))/5 + sin(alpha(t) + phi1(t))/2;
p_ee_base_z = (4*cos(phi1(t)))/5 + cos(alpha(t) + phi1(t))/2;

%The Jacobian matrix for the 2-DoF Arm
J11 = platform.l1*cos(phi1)+platform.l2*cos(phi1+alpha);
J12 = platform.l2*cos(phi1+alpha);
J21 = -platform.l1*sin(phi1)-platform.l2*sin(phi1+alpha);
J22 = -platform.l2*sin(phi1+alpha);

%Inverse Jacobian Matrix
det_J = J11*J22-J12*J21;
J11_inv = (1/det_J)*J22;
J12_inv = -(1/det_J)*J12;
J21_inv = -(1/det_J)*J21;
J22_inv = (1/det_J)*J11;

phi1_dot = J11_inv*diff(p_ee_base_x,1)+J12_inv*diff(p_ee_base_z,1);
alpha_dot = J21_inv*diff(p_ee_base_x,1)+J22_inv*diff(p_ee_base_z,1);

J_dot11 = diff(J11,phi1)*phi1_dot + diff(J11,alpha)*alpha_dot;
J_dot12 = diff(J12,phi1)*phi1_dot + diff(J12,alpha)*alpha_dot;
J_dot21 = diff(J21,phi1)*phi1_dot + diff(J21,alpha)*alpha_dot;
J_dot22 = diff(J22,phi1)*phi1_dot + diff(J22,alpha)*alpha_dot;

%(diff(p_ee_x,2)-J_dot11*phi1_dot+J_dot12*alpha_dot)
%(diff(p_ee_z,2)-J_dot21*phi1_dot+J_dot22*alpha_dot)

phi1_double_dot = J11_inv*(diff(p_ee_base_x,2)-J_dot11*phi1_dot+J_dot12*alpha_dot)+J12_inv*(diff(p_ee_base_z,2)-J_dot21*phi1_dot+J_dot22*alpha_dot);
alpha_double_dot = J21_inv*(diff(p_ee_base_x,2)-J_dot11*phi1_dot+J_dot12*alpha_dot)+J22_inv*(diff(p_ee_base_z,2)-J_dot21*phi1_dot+J_dot22*alpha_dot);

%Reference paramenters
ref_phi1 = 1.2126; %rad
ref_alpha = -2.3288; %rad
ref_phi1_dot = 0; %rad/s
ref_alpha_dot = 0; %rad/s
ref_phi1_double_dot = 0; %rad/s^2
ref_alpha_double_dot = 0; %rad/s^2
ref_phi1_int = 1.2126; %rad*s
ref_alpha_int = -2.3288; %rad*s

%Aceleration Control Gains
Kv_phi1 = 60;
Ka_phi1 = 100;
Kv_alpha = 30;
Ka_alpha = 60;

%PID Gains
Kp_phi1 = 5;
Kd_phi1 = 30;
Ki_phi1 = 20;
Kp_alpha = 5;
Kd_alpha = 5;
Ki_alpha = 5;

%defining two opposite trapezoidal functions

%defining two opposite square functions
sqr = 2*floor((1/2.1)*t)-floor((1/1.05)*t)+1;
sqr2 = -2*floor((1/2.1)*t)+floor((1/1.05)*t);

tau_c1 = sqr*(Kv_phi1*(ref_phi1_dot-phi1_dot)+Ka_phi1*(ref_phi1_double_dot-phi1_double_dot))+sqr2*(Kp_phi1*(ref_phi1-phi1)+Kd_phi1*diff(ref_phi1_dot-phi1)+Ki_phi1*int((ref_phi1_int-phi1),[0 t]));
tau_c2 = sqr*(Kv_alpha*(ref_alpha_dot-alpha_dot)+Ka_alpha*(ref_alpha_double_dot-alpha_double_dot))+sqr2*(Kp_alpha*(ref_alpha-alpha)+Kd_alpha*diff(ref_alpha_dot-alpha)+Ki_alpha*int((ref_alpha_int-alpha),[0 t]));

F_a2_x = -platform.m2*(diff(dx,2)+platform.l1*(cos(phi1)*diff(phi1,2)-sin(phi1)*(diff(phi1))^2)-platform.r2*(cos(phi1+alpha)*diff(phi1+alpha,2)-sin(phi1+alpha)*(diff(phi1+alpha))^2));
F_a2_z = -platform.m2*(diff(dz,2)-platform.l1*(cos(phi1)*(diff(phi1))^2 + sin(phi1)*diff(phi1,2))-platform.r2*(cos(phi1+alpha)*(diff(phi1+alpha))^2 + sin(phi1+alpha)*diff(phi1+alpha,2))+platform.g);

diffeq1 = tau_c2+platform.r2*(F_a2_x*cos(phi1+alpha)-F_a2_z*sin(phi1+alpha)) == platform.Iy2*diff(phi1+alpha,2);

F_d_x = platform.m1*(diff(dx,2)+platform.r1*(cos(phi1)*diff(phi1,2)-sin(phi1)*(diff(phi1))^2)) - F_a2_x;
F_d_z = platform.m1*(diff(dz,2)-platform.r1*(cos(phi1)*(diff(phi1))^2 + sin(phi1)*diff(phi1,2)))-F_a2_z+platform.m1*platform.g;

diffeq2 = tau_c1-tau_c2+platform.r1*(F_d_z*sin(phi1)-F_d_x*cos(phi1))+(platform.l1-platform.r1)*(F_a2_x*cos(phi1)-F_a2_z*sin(phi1)) == platform.Iy1*diff(phi1,2);

[V,S] = odeToVectorField(diffeq1,diffeq2);

M = matlabFunction(V,'vars',{'t','Y'});

%solving the system
freq = 160; %Hz
time_interval = [0 20]; %seconds
initial_conditions = [ref_phi1 0 ref_alpha 0]; %[phi1 phi1_dot alpha alpha_dot]

[sol.x sol.y]= rk4_solver(M,time_interval,initial_conditions,1/freq);
sol.y = sol.y';

figure()
subplot(2,2,1)
plot(sol.x,sol.y(1,:))
xlabel("Time [s]")
ylabel("\theta_1 [rad]")
title("\theta_1 Over Time")
yline(pi/2,'r')
yline(-pi/2,'r')
yline(ref_phi1, 'b--')

subplot(2,2,2)
plot(sol.x,sol.y(2,:))
xlabel("Time [s]")
ylabel("\theta_1^{dot} [rad/s]")
title("\theta_1^{dot} Over Time")

subplot(2,2,3)
plot(sol.x,sol.y(3,:))
xlabel("Time [s]")
ylabel("\theta_2 [rad]")
title("\theta_2 Over Time")
yline(pi,'r')
yline(-pi,'r')
yline(ref_alpha,'b--')

subplot(2,2,4)
plot(sol.x,sol.y(4,:))
xlabel("Time [s]")
ylabel("\theta_2^{dot} [rad/s]")
title("\theta_2^{dot} Over Time")

figure()
subplot(2,1,1)
plot(sol.x, double(subs(dx,sol.x)))
xlabel("Time [s]")
ylabel("X-axis [m]")
title("Deck X-Movement Over Time")

subplot(2,1,2)
plot(sol.x,double(subs(dz,sol.x)))
xlabel("Time [s]")
ylabel("Z-axis [m]")
title("Deck Z-Movement Over Time")


%define design choice
dxt = double(subs(dx,sol.x))';
dzt = double(subs(dz,sol.x))';
d = [dxt,dzt];


r_a1_star_x = ones(length(sol.x),1)*-2;
r_a1_star_z = ones(length(sol.x),1)*1;
r_a1_star = [r_a1_star_x,r_a1_star_z]; %design parameter

r_1_star = [platform.r1*sin(sol.y(1,:))', platform.r1*(cos(sol.y(1,:))-1)'];

p1 = d+r_a1_star+r_1_star; %inertial position of the cg of link 1

l1 = [platform.l1*sin(sol.y(1,:))',platform.l1*cos(sol.y(1,:))'];
l2 = [platform.l2*sin(sol.y(1,:)+sol.y(3,:))',platform.l2*cos(sol.y(1,:)+sol.y(3,:))'];
r1 = [platform.r1*sin(sol.y(1,:))',platform.r1*cos(sol.y(1,:))'];
r2 = [platform.r2*sin(sol.y(1,:)+sol.y(3,:))',platform.r2*cos(sol.y(1,:)+sol.y(3,:))'];

p2 = p1 + l1 - r1 + r2; %inertial position of the cg of link 2

p_ee = p1 + l1 - r1 + l2; %inertial position of the end effector
p_ee_relative = r_a1_star+r_1_star+l1-r1+l2; %relative position of the end effector wrt the deck cg

figure()
subplot(2,1,1)
plot(sol.x,p_ee(:,1))
xlabel("Time [s]")
ylabel("X-axis [m]")
title("End Effector Movement In The Inertial X Direction")

subplot(2,1,2)
plot(sol.x,p_ee(:,2))
xlabel("Time [s]")
ylabel("Z-axis [m]")
title("End Effector Movement In The Inertial Z Direction")

figure()
subplot(2,1,1)
plot(sol.x,p_ee_relative(:,1))
xlabel("Time [s]")
ylabel("X-axis [m]")
title("End Effector Relative Movement In The Inertial X Direction")

subplot(2,1,2)
plot(sol.x,p_ee_relative(:,2))
xlabel("Time [s]")
ylabel("Z-axis [m]")
title("End Effector Relative Movement In The Inertial Z Direction")

%calculating the inertial acceleration of the end effector
p_ee_x_dot = gradient(p_ee(:,1),sol.x);
p_ee_z_dot = gradient(p_ee(:,2),sol.x);
p_ee_x_double_dot = gradient(p_ee_x_dot,sol.x);
p_ee_z_double_dot = gradient(p_ee_z_dot,sol.x);

figure()
subplot(2,1,1)
plot(sol.x,p_ee_x_double_dot);
xlabel("Time [s]")
ylabel("Acceleration [m/s^2]")
title("The X Inertial Acceleration of The End Effector.")
subplot(2,1,2)
plot(sol.x,p_ee_z_double_dot);
xlabel("Time [s]")
ylabel("Acceleration [m/s^2]")
title("The Z Inertial Acceleration of The End Effector.")

%plotting two opposite square functions
sqr = 2*floor((1/2.1)*sol.x)-floor((1/1.05)*sol.x)+1;
sqr2 = -2*floor((1/2.1)*sol.x)+floor((1/1.05)*sol.x);

figure()
plot(sol.x,sqr)
hold on
plot(sol.x,sqr2)
xlabel("Time [s]")
ylabel("Amplitude [m/m]")
legend("Acceleration Control", "Relative Position Control")
title("Controls Implementation For Inertial Stability")

%plotting the control torques
phi1_double_dot = gradient(sol.y(1,:),sol.x);
alpha_double_dot = gradient(sol.y(3,:),sol.x);

tau_c1 = sqr.*(Kv_phi1*(ref_phi1_dot*ones(size(sol.y(1,:)))-sol.y(2,:))+Ka_phi1*(ref_phi1_double_dot*ones(size(sol.y(1,:)))-phi1_double_dot))+sqr2.*(Kp_phi1*(ref_phi1*ones(size(sol.y(1,:)))-sol.y(1,:))+Kd_phi1*(ones(size(sol.y(2,:)))-sol.y(2,:))+Ki_phi1*cumtrapz(ref_phi1_int*ones(size(sol.y(1,:)))-sol.y(1,:),sol.x));
tau_c2 = sqr.*(Kv_alpha*(ref_alpha_dot*ones(size(sol.y(3,:)))-sol.y(4,:))+Ka_alpha*(ref_alpha_double_dot*ones(size(sol.y(3,:)))-alpha_double_dot))+sqr2.*(Kp_alpha*(ref_alpha*ones(size(sol.y(3,:)))-sol.y(3,:))+Kd_alpha*(ones(size(sol.y(4,:)))-sol.y(4,:))+Ki_alpha*cumtrapz(ref_phi1_int*ones(size(sol.y(3,:)))-sol.y(3,:),sol.x));

figure()
subplot(2,1,1)
plot(sol.x,tau_c1)
xlabel("Time [s]")
ylabel("Torque [N.m]")
title("Control Torque of Joint 1 Over Time")

subplot(2,1,2)
plot(sol.x,tau_c2)
xlabel("Time [s]")
ylabel("Torque [N.m]")
title("Control Torque of Joint 2 Over Time")