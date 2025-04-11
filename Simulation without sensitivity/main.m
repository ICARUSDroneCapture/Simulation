% Contributors: Mohammed Al Alawi, Kirin Kawamoto
% Course number: ASEN 4028
% File name: main
% Created: 4/3/2025

% housekeeping
clear; clc; close all

% calling constants
SimulationParameters;

% solving the system
dt = 0.0001; %[d]
time_interval = [0 30]; %seconds

initial_conditions = [-0.1; 0; 0; 
                       0; 0; 0]; %[q1; Dq1; int_q1_err; 
                                  % pm_ddot];

MFun = @(t, y)EOM_V3(t, y, a);
% tic
% for i = 1:20
[sol.x, sol.y]= rk4_solver(MFun,time_interval,initial_conditions,dt);
% end
% toc
% sol.y = sol.y';

save('simResults.mat', 'sol')

%% Plotting and analysis
close all; clear; clc;

SimulationParameters;

load('simResults.mat')

m = 10;
x = sol.x(1:m:end);
y = sol.y(1:m:end,:)';

% plot q1,q1_dot,q2, q2_dot over time
q1_ref = -pi/4;
plot1(x,y,q1_ref)

syms t

% the base movement in the inertial frame 
dx = 0*t; %[m] (DO NOT CHANGE)
dy = 0*t; %[m] (DO NOT CHANGE)
dz = 0*t; %[m] (DO NOT CHANGE)

d =[dx;dy;dz];

% the deck rotations
period = 7.5*4; %[s]
angle2 = 90*(pi/180)*(sin((2*pi/(period))*t))^2; %deck rotation about its y-axis [rad]

% plot the inertial position of the deck and the platform
d_fun = matlabFunction(d, "Vars",{t});
B = d_fun(x);

theta2_fun = matlabFunction(angle2);
theta2_eval = theta2_fun(x);

xEE = B(1) + a.l1*cos(y(1,:)+theta2_eval);
zEE = B(3) - a.l1*sin(y(1,:)+theta2_eval);

% figure()
% plot(0,0,'rx',"LineWidth",2)
% hold on
% plot(xEE,zEE,'b-')
% xlabel("X-axis [m]")
% ylabel("Z-axis [m]")
% axis([-a.l1-0.1 a.l1+0.1 -a.l1-0.1 a.l1+0.1]);
% title("The System in The Inertial Frame")
% legend("Hinge Axis","End Effector")
% hold off

figure;
plot(x, zEE)
xlabel('Time (s)')
ylabel('Vertical Position (m)')
title('Inertial Position vs Time')

% calculating the performance of the inertial stability control
% lets find the index of x where the time is greater than 10 seconds
% timeIndex = find(x >= 10);
% max_zEE = max(zEE(timeIndex));
% min_zEE = min(zEE(timeIndex));
% 
% 
% IsolationPercent = ((max_zEE-min_zEE)/a.l1)*100;
% 
% fprintf('Inertial Stability Isolation Percent = %3.2f%% \n',IsolationPercent)
% 
% % animate
% animate(x,y, d, angle2, a)