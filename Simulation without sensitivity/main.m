% Contributors: Mohammed Al Alawi, Kirin Kawamoto
% Course number: ASEN 4028
% File name: main
% Created: 4/3/2025

% housekeeping
clear; clc; close all

% calling constants7
SimulationParameters;

% solving the system
dt = 0.001; %[d]
time_interval = [0 120]; %seconds

initial_conditions = [0; 0; 0]; %[q1; Dq1; int_q1_err; 
                                  % pm_ddot];

MFun = @(t, y)EOM_V3_5(t, y, a);
% tic
% for i = 1:20
[sol.x, sol.y]= rk4_solver(MFun,time_interval,initial_conditions, dt);
% end
% toc
% sol.y = sol.y';

save('simResults.mat', 'sol')

%% Plotting and analysis
close all; clear;

SimulationParameters;

load('simResults.mat')

m = 1;
x = sol.x(1:m:end);
y = sol.y(1:m:end,:)';

% plot q1,q1_dot,q2, q2_dot over time
plot1(x,y,a)

syms t

% the base movement in the inertial frame 
dx = 0*t; %[m] (DO NOT CHANGE)
dy = 0*t; %[m] (DO NOT CHANGE)
dz = 0*t; %[m] (DO NOT CHANGE)

d =[dx;dy;dz];

% plot the inertial position of the deck and the platform
d_fun = matlabFunction(d, "Vars",{t});
B = d_fun(x);

theta2_eval = a.thetad(x);

xEE = B(1) + a.l1*cos(y(1,:)+theta2_eval);
zEE = B(3) - a.l1*sin(y(1,:)+theta2_eval);

xNotIso = B(1) + a.l1*cos(a.q1_ref+theta2_eval);
zNotIso = B(3) - a.l1*sin(a.q1_ref+theta2_eval);
t_eval = x(end)*0.5;
t_idx = (x > t_eval);

isolation = calculateIsolationFFT(zEE(t_idx), zNotIso(t_idx));
fprintf('Average isolation: %0.2f%%\n', isolation*100);

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

subplot(2,1,1)
plot(x, zEE)
hold on
plot(x, zNotIso)
xlabel('Time (s)')
ylabel('Vertical Position (m)')
title('Inertial Position vs Time')
legend('Isolated', 'Not Isolated')

subplot(2,1,2)
plot(x, y(1, :) - a.q1_ref)
hold on
plot(x, y(3, :))
plot(x, y(2, :))
hold off
xlabel('Time (s)')
ylabel('Control States')
title('Control State Values vs Time')
legend( 'q_{err}', 'q_{err\_int}', 'q_{dot}')

% torques = zeros(5, length(x));
% for i = 1:length(x)
%     taus = torqueValues(x(i), y(:, i), a);
%     torques(:, i) = taus;
% end
% 
% subplot(3,1,3)
% hold on
% for i = 1:5
%     plot(x, torques(i, :))
% end
% hold off
% xlabel('Time (s)')
% ylabel('Toprque (Nm)')
% title('Torques vs Time')
% legend('K_a', 'K_p', 'K_i', 'K_d', 'Friction')


% Level of inertial and boundary control
% intert = a.C(y(1,:));
% bound = a.B(y(1,:));
% 
% figure;
% plot(x, intert)
% hold on
% plot(x, bound)
% xlabel('Time (s)')
% ylabel('Gain Proportion (m)')
% title('Inertial and Boundary Gain Proportions')
% legend('Inertial', 'Boundary')

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
% animate(x,y, d, a.thetad, a)