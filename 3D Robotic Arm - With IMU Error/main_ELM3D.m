% Contributors: Mohammed Al Alawi
% Course number: ASEN 4028
% File name: main_ELM3D
% Created: 1/22/2025


% housekeeping
clear; clc; close all

% constants
constants;

dir = "3DOFSensorError";

% ----------------------- Get IMU Data With Error -------------------------

DEMO_IMU
close all

% Cutoff first minute of acceleration data
imu_start = length(accel_m_controlled_z) / 4 + 1;
a.imu_data = [accel_m_controlled_x(imu_start:end)'; accel_m_controlled_y(imu_start:end)'; accel_m_controlled_z(imu_start:end)'];

% -------------------------------------------------------------------------

% deck movement (Should be the same as those in EOM3D.m)
syms t
period = 7.5; %[s]
xamplitude = 0; %wave amplitude [m]
yamplitude = 0; %wave amplitude [m]
zamplitude = 0.1; %wave amplitude [m]
dx = xamplitude*cos((2*pi/period)*t); %[m]
dy = yamplitude*cos((2*pi/period)*t); %[m]
dz = zamplitude*cos((2*pi/period)*t); %[m]

d =[dx;dy;dz];
d_ddot = diff(d,'t',2);

% the deck rotations (Should be the same as those in EOM3D.m)
angle1 = 0*t; % deck rotation about its x-axis [rad]
angle2 = 0*t; %10*(pi/180)*sin((2*pi/(period))*t); %deck rotation about its y-axis [rad]
angle3 = 0*t; %deck rotation about its z-axis [rad]

theta_D = [angle1;angle2;angle3];

% solving the system
freq = 160; %Hz
time_interval = [0 60*3]; %seconds

%---------------------------------------------
% Initial Conditions (for simulation)
%---------------------------------------------
q1_0 = -2.7925; %[rad]
q2_0 = 2.6180; %[rad]
q3_0 = 0; %[rad]
Dq1_0 = 0; %[rad/s]
Dq2_0 = 0; %[rad/s]
Dq3_0 = 0; %[rad/s]

%---------------------------------------------
% End of Initial Conditions (for simulation)
%---------------------------------------------

%---------------------------------------------
% Reference parameters (for control)
%---------------------------------------------
% reference angles
ref_q1 = -1.9199;
ref_q2 = 1.7453;
ref_q3 = 0;

% relative position of end effector in the base frame
r_B_ref = [0.4;0;0.85];
%---------------------------------------------
% End of Reference parameters (for control)
%---------------------------------------------

%--------------
% Gains
%--------------
n1 = 1;
n2 = 0;
% gains of inertial stability control
Ka = 2000*n1;
Kv = 750*n2;

n = 1;
% gains of relative position control
Kp = 12.6*n; 
Ki = 0.0001*n;
Kd = 85*n;

n = 0;
% gains of motor 1
Kp1 = 3*n;
Ki1 = 1*n;
Kd1 = 3*n;

% gains of motor 2
Kp2 = 1*n;
Ki2 = 0.5*n;
Kd2 = 1*n;

% gains of motor 3
Kp3 = 1*n;
Ki3 = 0.5*n;
Kd3 = 1*n;
%--------------
% End of Gains
%--------------

References = [ref_q1;ref_q2;ref_q3;r_B_ref(1);r_B_ref(2);r_B_ref(3)];
Gains =[Ka;Kv;Kp;Ki;Kd;Kp1;Ki1;Kd1;Kp2;Ki2;Kd2;Kp3;Ki3;Kd3];

% initial_conditions = [q1_0 q2_0 q3_0 Dq1_0 Dq2_0 Dq3_0 0 0 0 0 0 0]; %[q1 q2 q3 Dq1 Dq2 Dq3 e_rpX e_rpY e_rpZ e_j1 e_j2 e_j3]
initial_conditions = [q1_0 q2_0 q3_0 Dq1_0 Dq2_0 Dq3_0 0 0 0 0 0 0 0 0 0]; %[q1 q2 q3 Dq1 Dq2 Dq3 e_rpX e_rpY e_rpZ e_j1 e_j2 e_j3]

[sol.x, sol.y]= rk4_solver(@(t,x)EOM3D(t,x,platform,Gains,References, a),time_interval,initial_conditions,1/freq);
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

% calculating the performance of the inertial stability control
% times up to which isolation is measured (vector)
timeIsolation = period:period:time_interval(2)-period;

% initiaze the IsolationPercent variable
IsolationPercentx = zeros(size(timeIsolation));
IsolationPercenty = zeros(size(timeIsolation));
IsolationPercentz = zeros(size(timeIsolation));

for i = 1:length(timeIsolation)
% lets find the index of sol.x where the time is greater than timeIsolation(i) seconds
% 0.5 m: sea state 2 amplitude
timeIndex = find(sol.x >= timeIsolation(i));
max_xEE = max(xEE(timeIndex));
min_xEE = min(xEE(timeIndex));

IsolationPercentx(i) = ((max_xEE-min_xEE)/0.5)*100;

max_yEE = max(yEE(timeIndex));
min_yEE = min(yEE(timeIndex));

IsolationPercenty(i) = ((max_yEE-min_yEE)/0.5)*100;

max_zEE = max(zEE(timeIndex));
min_zEE = min(zEE(timeIndex));

IsolationPercentz(i) = ((max_zEE-min_zEE)/0.5)*100;
end

% plotting the inertial positions of deck and end-effector
fig = figure;
axis equal
plot3(xB, yB, zB,"LineWidth",1.2);
xlabel("X-axis [m]",'FontWeight','bold')
ylabel("Y-axis [m]",'FontWeight','bold')
zlabel("Z-axis [m]",'FontWeight','bold')
title("The system in The Inertial Frame")
axis([-platform.l1-platform.l2-0.02-6 platform.l1+platform.l2+0.02+6 -platform.l1-platform.l2-0.02 platform.l1+platform.l2+0.02 -platform.l1-platform.l2-0.02-0.5 platform.l1+platform.l2+0.02+0.5]);
hold on
plot3(xEE, yEE, zEE,"LineWidth",1.2)
legend("Deck","End-Effector")
grid on
view(-55,20)
hold off

saveas(fig, "../figures/" + dir + "/end_effector_inertial_pos.png")

% plotting the inertial acceleration of the deck
d_ddot = eval(subs(d_ddot,'t',sol.x));
fig = figure;
subplot(3,1,1)
plot(sol.x,d_ddot(1,:),"LineWidth",1.2);
xlabel("Time [s]",'FontWeight','bold')
ylabel("Inertial X Acc. [m/s^{s}]",'FontWeight','bold')
title("Deck Acceleration")
subplot(3,1,2)
plot(sol.x,d_ddot(2,:),"LineWidth",1.2);
xlabel("Time [s]",'FontWeight','bold')
ylabel("Inertial Y Acc. [m/s^{s}]",'FontWeight','bold')
subplot(3,1,3)
plot(sol.x,d_ddot(3,:),"LineWidth",1.2);
xlabel("Time [s]",'FontWeight','bold')
ylabel("Inertial Z Acc. [m/s^{s}]",'FontWeight','bold')

saveas(fig, "../figures/" + dir + "/end_effector_accel.png")

% plotting the X inertial position of the deck and end effector as well as
% the isolation percentage
fig = figure;
subplot(2,1,1)
plot(sol.x,xEE,"LineWidth",1.2)
xlabel("Time [s]",'FontWeight','bold')
ylabel("X-axis [m]",'FontWeight','bold')
title("Inertial X Position of the Deck & End-Effector")
hold on 
plot(sol.x,xB,"LineWidth",1.2)
hold off
legend("End Effector","Deck")

subplot(2,1,2)
plot(timeIsolation,IsolationPercentx,"x",'MarkerSize',6,"LineWidth",2)
xlabel("Time [seconds]",'FontWeight','bold')
ylabel("Isolation Percentage [%]",'FontWeight','bold')
hold on
yline(15,'r--',"LineWidth",2)
hold off
legend("Isolation Points","Required Percentage")
title("Isolation Percentage (in X-Axis) at Various Time Points to Finish Time")

saveas(fig, "../figures/" + dir + "/x_isolation.png")

% plotting the Y inertial position of the deck and end effector as well as
% the isolation percentage
fig = figure;
subplot(2,1,1)
plot(sol.x,yEE,"LineWidth",1.2)
xlabel("Time [s]",'FontWeight','bold')
ylabel("Y-axis [m]",'FontWeight','bold')
title("Inertial Y Position of the Deck & End-Effector")
hold on 
plot(sol.x,yB,"LineWidth",1.2)
hold off
legend("End Effector","Deck")

subplot(2,1,2)
plot(timeIsolation,IsolationPercenty,"x",'MarkerSize',6,"LineWidth",2)
xlabel("Time [seconds]",'FontWeight','bold')
ylabel("Isolation Percentage [%]",'FontWeight','bold')
hold on
yline(15,'r--',"LineWidth",2)
hold off
legend("Isolation Points","Required Percentage")
title("Isolation Percentage (in Y-Axis) at Various Time Points to Finish Time")

saveas(fig, "../figures/" + dir + "/y_isolation.png")

% plotting the Z inertial position of the deck and end effector as well as
% the isolation percentage
fig = figure;
subplot(2,1,1)
plot(sol.x,zEE,"LineWidth",1.2)
xlabel("Time [s]",'FontWeight','bold')
ylabel("Z-axis [m]",'FontWeight','bold')
title("Inertial Z Position of the Deck & End-Effector")
hold on 
plot(sol.x,zB,"LineWidth",1.2)
hold off
legend("End Effector","Deck")

subplot(2,1,2)
plot(timeIsolation,IsolationPercentz,"x",'MarkerSize',6,"LineWidth",2)
xlabel("Time [seconds]",'FontWeight','bold')
ylabel("Isolation Percentage [%]",'FontWeight','bold')
hold on
yline(15,'r--',"LineWidth",2)
hold off
legend("Isolation Points","Required Percentage")
title("Isolation Percentage (in Z-Axis) at Various Time Points to Finish Time")

saveas(fig, "../figures/" + dir + "/z_isolation.png")

% torque values
tau = zeros(3,length(sol.x));
for i=1:length(sol.x)
[~,tau(:,i)] = EOM3D(sol.x(i),sol.y(:,i),platform,Gains,References, a);
end

% plotting the torques of each joint
fig = figure;
subplot(3,1,1)
plot(sol.x,tau(1,:),"LineWidth",1.2)
xlabel('Time [s]','FontWeight','bold')
ylabel('Input Torque [N.m]','FontWeight','bold')
title('Input Torque of Joint 1')
subplot(3,1,2)
plot(sol.x,tau(2,:),"LineWidth",1.2)
xlabel('Time [s]','FontWeight','bold')
ylabel('Input Torque [N.m]','FontWeight','bold')
title('Input Torque of Joint 2')
subplot(3,1,3)
plot(sol.x,tau(3,:),"LineWidth",1.2)
xlabel('Time [s]','FontWeight','bold')
ylabel('Input Torque [N.m]','FontWeight','bold')
title('Input Torque of Joint 3')

saveas(fig, "../figures/" + dir + "/torque_values.png")

% animate
animate3D(sol.x,sol.y,platform,d,theta_D)