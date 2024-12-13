% Contributors: Mohammed Al Alawi
% Course number: ASEN 4018
% File name: main2
% Created: 10/13/2024


%HouseKeeping
clear; clc; close all

% calling constants
constants;

%solving the inverse kinematic problem
% desiredPosition = [0.3,0.5]; %[m]
% fun = @(theta)forwardKinematics(theta,desiredPosition);
% theta0 = [0,0];
% theta = fsolve(fun,theta0)

syms phi1(t) alpha(t)

dx = 0.2*t;
peroid = 7.5; %s
dz = 0.4*sin((2*pi/peroid)*t);

%Reference paramenters
ref_phi1 = 1.2126; %rad
ref_alpha = -2.3288; %rad
ref_phi1_dot = 0; %rad/s
ref_alpha_dot = 0; %rad/s
ref_phi1_int = 1.2126; %rad*s
ref_alpha_int = -2.3288; %rad*s

%PID Gains
Kp_phi1 = 5;
Kd_phi1 = 30;
Ki_phi1 = 20;
Kp_alpha = 5;
Kd_alpha = 5;
Ki_alpha = 5;

tau_c1 = Kp_phi1*(ref_phi1-phi1)+Kd_phi1*diff(ref_phi1_dot-phi1)+Ki_phi1*int((ref_phi1_int-phi1),[0 t]);
tau_c2 = Kp_alpha*(ref_alpha-alpha)+Kd_alpha*diff(ref_alpha_dot-alpha)+Ki_alpha*int((ref_alpha_int-alpha),[0 t]);

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
time_interval = [0 10]; %seconds

initial_conditions = [-(4/9)*pi 0 (8/9)*pi 0]; %[phi1 phi1_dot alpha alpha_dot]

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

%Lets discover how jocobian is doing 
%The jocobian matrix for the 2-DoF Arm
J = [platform.l1*cos(phi1)+platform.l2*cos(phi1+alpha),platform.l2*cos(phi1+alpha);-platform.l1*sin(phi1)-platform.l2*sin(phi1+alpha), -platform.l2*sin(phi1+alpha)];

%need to evaluate the deck forces in order to be able to evaluate tau = J'*F
tau_c1 = Kp_phi1*(ref_phi1*ones(size(sol.y(1,:)))-sol.y(1,:))+Kd_phi1*(ones(size(sol.y(2,:)))-sol.y(2,:))+Ki_phi1*cumtrapz(ref_phi1_int*ones(size(sol.y(1,:)))-sol.y(1,:),sol.x);
tau_c2 = Kp_alpha*(ref_alpha*ones(size(sol.y(3,:)))-sol.y(3,:))+Kd_alpha*(ones(size(sol.y(4,:)))-sol.y(4,:))+Ki_alpha*cumtrapz(ref_phi1_int*ones(size(sol.y(3,:)))-sol.y(3,:),sol.x);

phi1_double_dot = gradient(sol.y(1,:),sol.x);
alpha_double_dot = gradient(sol.y(3,:),sol.x);

F_a2_x = -platform.m2*(double(subs(diff(dx,2),sol.x))+platform.l1*(cos(sol.y(1,:)).*phi1_double_dot-sin(sol.y(1,:)).*(sol.y(2,:)).^2)-platform.r2*(cos(sol.y(1,:)+sol.y(3,:)).*(phi1_double_dot+alpha_double_dot)-sin(sol.y(1,:)+sol.y(3,:)).*(sol.y(2,:)+sol.y(4,:)).^2));
F_a2_z = -platform.m2*(double(subs(diff(dz,2),sol.x))-platform.l1*(cos(sol.y(1,:)).*(sol.y(2,:)).^2 + sin(sol.y(1,:)).*phi1_double_dot)-platform.r2*(cos(sol.y(1,:)+sol.y(3,:)).*(sol.y(2,:)+sol.y(4,:)).^2 + sin(sol.y(1,:)+sol.y(3,:)).*(phi1_double_dot+alpha_double_dot))+platform.g);
 
F_d_x = platform.m1*(double(subs(diff(dx,2),sol.x))+platform.r1*(cos(sol.y(1,:)).*sol.y(2,:)-sin(sol.y(1,:)).*(sol.y(2,:)).^2)) - F_a2_x;
F_d_z = platform.m1*(double(subs(diff(dz,2),sol.x))-platform.r1*(cos(sol.y(1,:)).*(sol.y(2,:)).^2 + sin(sol.y(1,:)).*phi1_double_dot))-F_a2_z+(platform.m1*platform.g)*ones(size(sol.x));

figure()
subplot(2,1,1)
plot(sol.x,F_d_x)
xlabel("Time [s]")
ylabel("Force [N]")
title("Deck X-Forces on Arm Over Time")

subplot(2,1,2)
plot(sol.x,F_d_z)
xlabel("Time [s]")
ylabel("Force [N]")
title("Deck Z-Forces on Arm Over Time")

%Now lets test the jocobian of the reference condition and its effects on the
%required control torques to hold joints stationary

% J_ref = [platform.l1*cos(ref_phi1)+platform.l2*cos(ref_phi1+ref_alpha),platform.l2*cos(ref_phi1+ref_alpha);-platform.l1*sin(ref_phi1)-platform.l2*sin(ref_phi1+ref_alpha), -platform.l2*sin(ref_phi1+ref_alpha)];
% 
% torque_of_joints = J_ref'*[F_d_x;F_d_z-ones(size(F_d_z)).*((platform.m1+platform.m2)*platform.g)];

figure()
subplot(2,1,1)
plot(sol.x,tau_c1)
hold on
%plot(sol.x,torque_of_joints(1,:))
xlabel("Time [s]")
ylabel("Torque [N.m]")
title("Control Torque of Joint 1 Over Time")
%legend("Actual Dynamics", "Jacobian Approach")
hold off

subplot(2,1,2)
plot(sol.x,tau_c2)
hold on
%plot(sol.x,torque_of_joints(2,:))
xlabel("Time [s]")
ylabel("Torque [N.m]")
title("Control Torque of Joint 2 Over Time")
%legend("Actual Dynamics", "Jacobian Approach")
hold off
