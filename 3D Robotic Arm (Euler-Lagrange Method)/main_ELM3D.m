% Contributors: Mohammed Al Alawi
% Course number: ASEN 4028
% File name: main_ELM3D
% Created: 1/22/2025


% housekeeping
clear; clc; close all

% constants
constants;

% deck movement
syms t
dx = 0*t; %[m]
dy = 0*t; %[m]
period = 7.5; %[s]
dz = 0.5*cos((2*pi/period)*t); %[m]

d =[dx;dy;dz];
d_ddot = diff(d,'t',2);

% the deck rotations
angle1 = 0*t; % deck rotation about its x-axis [rad]
angle2 = 0*t;%10*(pi/180)*sin((2*pi/(period*2))*t); %deck rotation about its y-axis [rad]
angle3 = 0*t; %deck rotation about its z-axis [rad]

theta_D = [angle1;angle2;angle3];

% initializing integral error (for any integral controller used)
integral_error = [0;0;0];

% solving the system
freq = 160; %Hz
time_interval = [0 30]; %seconds

%---------------------------------------------
% Reference parameters (for control)
%---------------------------------------------
% reference angles
ref_q1 = (-0.097971-0.3)*1;
ref_q2 = (-1.9315-0.5)*1;
ref_q3 = 0;

% relative position of end effector in the base frame
r_B_ref = [0.5;0;1.2];
%---------------------------------------------
% End of Reference parameters (for control)
%---------------------------------------------

%--------------
% Gains
%--------------
n1 = 0;
n2 = 0;
% gains of inertial stability control
Ka = 600*n1;
Kv = 750*n2;

n = 1;
% gains of relative position control
Kp = 1.23*n;
Ki = 2*n;
Kd = 3*n;

n = 0;
% gains of motor 1
Kp1 = 6*n;
Ki1 = 3*n;
Kd1 = 6*n;

% gains of motor 2
Kp2 = 3*n;
Ki2 = 1*n;
Kd2 = 3*n;

% gains of motor 3
Kp3 = 5*n;
Ki3 = 3*n;
Kd3 = 3*n;
%--------------
% End of Gains
%--------------

References = [ref_q1;ref_q2;ref_q3;r_B_ref(1);r_B_ref(2);r_B_ref(3)];
Gains =[Ka;Kv;Kp;Ki;Kd;Kp1;Ki1;Kd1;Kp2;Ki2;Kd2;Kp3;Ki3;Kd3];

initial_conditions = [-0.097971-0.3 -1.9315-0.5 0 0 0 0 0 0 0 0 0 0]; %[q1 q2 q3 Dq1 Dq2 Dq3 e_rpX e_rpY e_rpZ e_j1 e_j2 e_j3]

[sol.x, sol.y]= rk4_solver(@(t,x)EOM3D(t,x,Gains,References),time_interval,initial_conditions,1/freq);
sol.y = sol.y';

% plot q3, q3_dot, q1, q1_dot, q2, q2_dot over time
plot2(sol.x,sol.y)

% plotting the inertial position of the deck and the platform
% base translations
xB = double(subs(d(1),t,sol.x));
yB = double(subs(d(2),t,sol.x));
zB = double(subs(d(3),t,sol.x));

% base rotations
d1 = double(subs(theta_D(1),t,sol.x));
d2 = double(subs(theta_D(2),t,sol.x));
d3 = double(subs(theta_D(3),t,sol.x));

% forward kinematics
xEE = xB + cos(d2).*cos(d3).*(platform.l1*cos(sol.y(1,:)).*cos(sol.y(3,:)) + platform.l2*cos(sol.y(3,:)).*cos(sol.y(1,:) + sol.y(2,:))) - (sin(d1).*sin(d3) + cos(d1).*cos(d3).*sin(d2)).*(platform.l1*sin(sol.y(1,:)) + platform.l2*sin(sol.y(1,:) + sol.y(2,:))) - (cos(d1).*sin(d3) - cos(d3).*sin(d1).*sin(d2)).*(platform.l1*cos(sol.y(1,:)).*sin(sol.y(3,:)) + platform.l2*sin(sol.y(3,:)).*cos(sol.y(1,:) + sol.y(2,:)));
yEE = yB + (cos(d1).*cos(d3) + sin(d1).*sin(d2).*sin(d3)).*(platform.l1*cos(sol.y(1,:)).*sin(sol.y(3,:)) + platform.l2*sin(sol.y(3,:)).*cos(sol.y(1,:) + sol.y(2,:))) + (cos(d3).*sin(d1) - cos(d1).*sin(d2).*sin(d3)).*(platform.l1*sin(sol.y(1,:)) + platform.l2*sin(sol.y(1,:) + sol.y(2,:))) + cos(d2).*sin(d3).*(platform.l1*cos(sol.y(1,:)).*cos(sol.y(3,:)) + platform.l2*cos(sol.y(3,:)).*cos(sol.y(1,:) + sol.y(2,:)));
zEE = zB + cos(d2).*sin(d1).*(platform.l1*cos(sol.y(1,:)).*sin(sol.y(3,:)) + platform.l2*sin(sol.y(3,:)).*cos(sol.y(1,:) + sol.y(2,:))) - sin(d2).*(platform.l1*cos(sol.y(1,:)).*cos(sol.y(3,:)) + platform.l2*cos(sol.y(3,:)).*cos(sol.y(1,:) + sol.y(2,:))) - cos(d1).*cos(d2).*(platform.l1*sin(sol.y(1,:)) + platform.l2*sin(sol.y(1,:) + sol.y(2,:)));

% plotting the inertial positions of deck and end-effector
figure()
axis equal
plot3(xB, yB, zB);
xlabel("X-axis [m]")
ylabel("Y-axis [m]")
zlabel("Z-axis [m]")
title("The system in The Inertial Frame")
axis([-platform.l1-platform.l2-0.02-6 platform.l1+platform.l2+0.02+6 -platform.l1-platform.l2-0.02 platform.l1+platform.l2+0.02 -platform.l1-platform.l2-0.02-0.5 platform.l1+platform.l2+0.02+0.5]);
hold on
plot3(xEE, yEE, zEE)
legend("Deck","End-Effector")
grid on
view(-55,20)
hold off

% plotting the inertial acceleration of the deck
d_ddot = eval(subs(d_ddot,'t',sol.x));
figure()
subplot(3,1,1)
plot(sol.x,d_ddot(1,:));
xlabel("Time [s]")
ylabel("Inertial X Acc. [m/s^{s}]")
title("Deck Acceleration")
subplot(3,1,2)
plot(sol.x,d_ddot(2,:));
xlabel("Time [s]")
ylabel("Inertial Y Acc. [m/s^{s}]")
subplot(3,1,3)
plot(sol.x,d_ddot(3,:));
xlabel("Time [s]")
ylabel("Inertial Z Acc. [m/s^{s}]")


% torque values
tau = zeros(3,length(sol.x));
for i=1:length(sol.x)
[~,tau(:,i)] = EOM3D(sol.x(i),sol.y(:,i),Gains,References);
end

% plotting the torques of each joint
figure()
subplot(3,1,1)
plot(sol.x,tau(1,:))
xlabel('Time [s]')
ylabel('Input Torque [N.m]')
title('Input Torque of Joint 1')
subplot(3,1,2)
plot(sol.x,tau(2,:))
xlabel('Time [s]')
ylabel('Input Torque [N.m]')
title('Input Torque of Joint 2')
subplot(3,1,3)
plot(sol.x,tau(3,:))
xlabel('Time [s]')
ylabel('Input Torque [N.m]')
title('Input Torque of Joint 3')

% calculating the performance of the inertial stability control
% lets find the index of sol.x where the time is greater than 10 seconds
timeIndex = find(sol.x >= 10);
max_zEE = max(zEE(timeIndex));
min_zEE = min(zEE(timeIndex));


IsolationPercent = ((max_zEE-min_zEE)/0.5)*100;

fprintf('Inertial Stability Isolation Percent = %3.2f%% \n',IsolationPercent)

% animate
animate3D(sol.x,sol.y,d, theta_D)