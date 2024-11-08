% Contributors: Mohammed Al Alawi
% Course number: ASEN 4018
% File name: main
% Created: 10/13/2024


%HouseKeeping
clear; clc; close all

% calling constants
constants;

syms phi1(t) phi2(t)

dx = 0.2*t;
dz = sin(0.25*t);

%Reference paramenters
ref_phi1 = pi/3; %rad
ref_phi2 = pi/12; %rad
ref_phi1_dot = 0; %rad/s
ref_phi2_dot = 0; %rad/s
ref_phi1_double_dot = 0; %rad/s^2
ref_phi2_double_dot = 0; %rad/s^2

%PID Gains
Kp_phi1 = 50;
Kd_phi1 = 10;
Ki_phi1 = 0;
Kp_phi2 = 10;
Kd_phi2 = 5;
Ki_phi2 = 0.1;

tau_c1 = Kp_phi1*(ref_phi1-phi1)+Kd_phi1*diff(ref_phi1_dot-phi1)+Ki_phi1*int(ref_phi1_double_dot-phi1,[0 t]);
tau_c2 = Kp_phi2*(ref_phi2-phi2)+Kd_phi2*diff(ref_phi2_dot-phi2)+Ki_phi2*int(ref_phi2_double_dot-phi2,[0 t]);

F_a2_x = -platform.m2*(diff(dx,2)+platform.l1*(cos(phi1)*diff(phi1,2)-sin(phi1)*(diff(phi1))^2)-platform.r2*(cos(phi2)*diff(phi2,2)-sin(phi2)*(diff(phi2))^2));
F_a2_z = -platform.m2*(diff(dz,2)-platform.l1*(cos(phi1)*(diff(phi1))^2 + sin(phi1)*diff(phi1,2))-platform.r2*(cos(phi2)*(diff(phi2))^2 + sin(phi2)*diff(phi2,2))+platform.g);

diffeq1 = tau_c2+platform.r2*(F_a2_x*cos(phi2)-F_a2_z*sin(phi2)) == platform.Iy2*diff(phi2,2);

F_d_x = platform.m1*(diff(dx,2)+platform.r1*(cos(phi1)*diff(phi1,2)-sin(phi1)*(diff(phi1))^2)) - F_a2_x;
F_d_z = platform.m1*(diff(dz,2)-platform.r1*(cos(phi1)*(diff(phi1))^2 + sin(phi1)*diff(phi1,2)))-F_a2_z+platform.m1*platform.g;

diffeq2 = tau_c1-tau_c2+platform.r1*(F_d_z*sin(phi1)-F_d_x*cos(phi1))+(platform.l1-platform.r1)*(F_a2_x*cos(phi1)-F_a2_z*sin(phi1)) == platform.Iy1*diff(phi1,2);

[V,S] = odeToVectorField(diffeq1,diffeq2);

M = matlabFunction(V,'vars',{'t','Y'});

%solving the system
time_interval = [0 10];
initial_conditions = [pi/8 0 pi/5 0]; %[phi1 phi1_dot phi2 phi2_dot]

sol = ode45(M,time_interval,initial_conditions);

figure()
subplot(2,2,1)
plot(sol.x,sol.y(1,:))
xlabel("Time [s]")
ylabel("\phi_1 [rad]")
title("\phi_1 Over Time")
yline(pi/2,'r')
yline(-pi/2,'r')
yline(ref_phi1, 'b--')

subplot(2,2,2)
plot(sol.x,sol.y(2,:))
xlabel("Time [s]")
ylabel("\phi_1^{dot} [rad/s]")
title("\phi_1^{dot} Over Time")

subplot(2,2,3)
plot(sol.x,sol.y(3,:))
xlabel("Time [s]")
ylabel("\phi_2 [rad]")
title("\phi_2 Over Time")
yline(pi,'r')
yline(-pi,'r')
yline(ref_phi2,'b--')

subplot(2,2,4)
plot(sol.x,sol.y(4,:))
xlabel("Time [s]")
ylabel("\phi_2^{dot} [rad/s]")
title("\phi_2^{dot} Over Time")

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
l2 = [platform.l2*sin(sol.y(3,:))',platform.l2*cos(sol.y(3,:))'];
r1 = [platform.r1*sin(sol.y(1,:))',platform.r1*cos(sol.y(1,:))'];
r2 = [platform.r2*sin(sol.y(3,:))',platform.r2*cos(sol.y(3,:))'];

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
