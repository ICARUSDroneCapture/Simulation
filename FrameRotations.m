close all; clear; clc;

rng(1,"twister");

% Use same random seed

set(groot,'DefaultLineLineWidth',1)

simulationParameters;

close all;

% Acceleration curves (and derived velocity/position) for deck motion (wave environment)

a.real_pos_xI = @(t) 0.4/(beta^2)*sin(beta*t/2);
a.real_pos_yI = @(t) 1.6/(beta^2)*sin(beta*t/4);
a.real_pos_zI = @(t) alpha*sin(beta*t) + hdeck;

a.real_vel_xI = @(t) 0.2/beta*cos(beta/2*t);
a.real_vel_yI = @(t) 0.4/beta*cos(beta/4*t);
a.real_vel_zI = @(t) beta*alpha*cos(beta*t);

a.real_accel_xI = @(t) -0.1*sin(beta/2*t); % [m*s^-2]
a.real_accel_yI = @(t) -0.1*sin(beta/4*t); % [m*s^-2]
a.real_accel_zI = @(t) -beta^2*alpha*sin(beta*t); % [m*s^-2]

a.theta = @(t) -atan(beta*alpha*cos(beta*t)); % [rad]
a.phi = @(t) atan(0.2/beta*cos(beta/2*t)); % [rad]
a.psi = @(t) atan(0.4/beta*cos(beta/4*t)); % [rad]

a.theta_dot = @(t) (-(alpha*beta^2*sin(beta*t))./(alpha^2*beta^2*(cos(beta*t).^2)+1)); % [rad/s]
a.phi_dot = @(t) ((-0.1*sin(beta/2*t))/(((0.04*(cos(beta*t/2)^2))/(beta^2))+1)); % [rad/s]
a.psi_dot = @(t) ((-0.1*sin(beta/4*t))/(((0.16*(cos(beta*t/4)^2))/(beta^2))+1)); % [rad/s]



% a.real_pos_xI = @(t) 0*t;
% a.real_pos_yI = @(t) 0*t;
% a.real_pos_zI = @(t) alpha*sin(beta*t) + hdeck;
% 
% a.real_vel_xI = @(t) 0*t;
% a.real_vel_yI = @(t) 0*t;
% a.real_vel_zI = @(t) beta*alpha*cos(beta*t);
% 
% a.real_accel_xI = @(t) 0*t; % [m*s^-2]
% a.real_accel_yI = @(t) 0*t; % [m*s^-2]
% a.real_accel_zI = @(t) -beta^2*alpha*sin(beta*t); % [m*s^-2]
% 
% a.theta = @(t) atan(beta*alpha*cos(beta*t)); % [rad]
% a.phi = @(t) atan(0.2/beta*cos(beta/2*t)); % [rad]
% a.psi = @(t) atan(0.4/beta*cos(beta/4*t)); % [rad]
% 
% a.theta_dot = @(t) (-(alpha*beta^2*sin(beta*t))./(alpha^2*beta^2*(cos(beta*t).^2)+1)); % [rad/s]
% a.phi_dot = @(t) ((-0.1*sin(beta/2*t))/(((0.04*(cos(beta*t/2)^2))/(beta^2))+1)); % [rad/s]
% a.psi_dot = @(t) ((-0.1*sin(beta/4*t))/(((0.16*(cos(beta*t/4)^2))/(beta^2))+1)); % [rad/s]

s = @(x) sin(x);
c = @(x) cos(x);

%% Sensor Model Aspects

% Simulation time
startTime = 0;
finishTime = 10;
tspan = [startTime finishTime]; % [s]

% dt = 1/imu_rate;  % [s]
dt = 0.1;
t = (tspan(1):dt:tspan(2))';
t_count = length(t);
indeces = @(t) floor(t/dt)+1;

% defineSignals
defineSignalsNoNoise

%% Plot 3D Wave Motion

figure
plot3(a.real_pos_xI(t), a.real_pos_yI(t), a.real_pos_zI(t))
xlim([-1 1])
ylim([0 2.5])
zlim([0.5 1.5])
xlabel('East')
ylabel('North')
zlabel('Up')

% figure
% xlabel('East')
% ylabel('North')
% zlabel('Up')
% 
% for i = 1:length(t)
% 
%     t_i = t(i);
% 
%     x = a.real_pos_xI(t_i);
%     y = a.real_pos_yI(t_i);
%     z = a.real_pos_zI(t_i);
% 
%     scatter3(x,y,z, 10, 'filled', "MarkerEdgeColor", "black", "MarkerFaceColor","black")
%     xlim([-1 1])
%     ylim([-2.5 0])
%     zlim([0.5 1.5])
%     drawnow
%     hold on
%     pause(dt)
% end

% figure
% xlabel('East')
% ylabel('North')
% zlabel('Up')
% 
% for i = 1:length(t)
% 
%     t_i = t(i);
% 
%     x = a.real_pos_xI(t_i);
%     y = a.real_pos_yI(t_i);
%     z = a.real_pos_zI(t_i);
% 
% 
%     scatter3(x, y, z, 10, 'filled', "MarkerEdgeColor", "black", "MarkerFaceColor","black")
%     hold on
% 
%     a_vec = [0 0 -9.8];
%     a_coord = a_vec/norm(a_vec)/10;
% 
%     frame_vec = [x y z] + a_coord;
% 
%     plot3([x frame_vec(1)], [y y], [z z], Color='magenta')
%     hold on
% 
%     plot3([x x], [y frame_vec(2)], [z z], Color='magenta')
%     hold on
% 
%     plot3([x x], [y y], [z frame_vec(3)], Color='magenta')
% 
%     xlim([-1 1])
%     ylim([-2.5 0])
%     zlim([0.5 1.5])
%     xlabel('East')
%     ylabel('North')
%     zlabel('Up')
%     drawnow
%     hold on
%     pause(dt)
% end

%% See position curves on all planes

figure

subplot(3, 1, 1)
plot(a.real_pos_xI(t), a.real_pos_zI(t))
title("Platform/Deck Position in each Plane")
xlabel('East (m)')
ylabel('Up (m)')

subplot(3, 1, 2)
plot(a.real_pos_yI(t), a.real_pos_zI(t))
xlabel('North (m)')
ylabel('Up (m)')

subplot(3, 1, 3)
plot(a.real_pos_xI(t), a.real_pos_yI(t))
xlabel('East (m)')
ylabel('North (m)')


figure

subplot(3, 1, 1)
plot(t, a.real_pos_zI(t))
title("Platform/Deck Position over Time")
xlabel('Time (sec)')
ylabel('Up (m)')

subplot(3, 1, 2)
plot(t, a.real_pos_yI(t))
xlabel('Time (sec)')
ylabel('North (m)')

subplot(3, 1, 3)
plot(t, a.real_pos_xI(t))
xlabel('Time (sec)')
ylabel('East (m)')


%% See acceleration curves on all planes

figure

subplot(3, 1, 1)
plot(a.real_accel_xI(t), a.real_accel_zI(t))
title("Platform/Deck Acceleration in each Plane")
xlabel('East (m/s^2)')
ylabel('Up (m/s^2)')

subplot(3, 1, 2)
plot(a.real_accel_yI(t), a.real_accel_zI(t))
xlabel('North (m/s^2)')
ylabel('Up (m/s^2)')

subplot(3, 1, 3)
plot(a.real_accel_xI(t), a.real_accel_yI(t))
xlabel('East (m/s^2)')
ylabel('North (m/s^2)')


figure

subplot(3, 1, 1)
plot(t, a.real_accel_zI(t))
title("Platform/Deck Acceleration over Time")
xlabel('Time (sec)')
ylabel('Up Acceleraiton (m/s^2)')

subplot(3, 1, 2)
plot(t, a.real_accel_yI(t))
xlabel('Time (sec)')
ylabel('North Acceleraiton (m/s^2)')

subplot(3, 1, 3)
plot(t, a.real_accel_xI(t))
xlabel('Time (sec)')
ylabel('East Acceleraiton (m/s^2)')

%% Plot Velocities over Time

figure

subplot(3, 1, 1)
plot(t, a.real_vel_zI(t))
title("Platform/Deck Velocity in each Plane")
xlabel('Time (sec)')
ylabel('Up Velocity (m/s)')

subplot(3, 1, 2)
plot(t, a.real_vel_yI(t))
ylabel('North Velocity (m/s)')
xlabel('Time (sec)')
ylabel('Angle (deg)')

subplot(3, 1, 3)
plot(t, a.real_vel_xI(t))
ylabel('East Velocity (m/s)')
xlabel('Time (sec)')
ylabel('Angle (deg)')

%% Plot Angles over Time

figure

subplot(3, 1, 1)
plot(t, a.theta(t) / pi*180)
title("Platform/Deck Angles in each Plane")
xlabel('Time (sec)')
ylabel('Pitch (theta) (deg)')

subplot(3, 1, 2)
plot(t, a.phi(t) / pi*180)
xlabel('Time (sec)')
ylabel('Roll (phi) (deg)')

subplot(3, 1, 3)
plot(t, a.psi(t) / pi*180)
xlabel('Time (sec)')
ylabel('Yaw (psi) (deg)')


%% Checko

theta = 90/180*pi;
phi = 90/180*pi;
psi = 0;

% R_I_S = [c(phi)*c(psi)                                  c(phi)*sin(psi)                             -sin(phi);
%              s(theta)*s(phi)*cos(psi)-c(theta)*s(psi)       s(theta)*s(phi)*s(psi)+c(theta)*c(psi)      s(theta)*c(phi);
%              c(theta)*s(phi)*c(psi)+s(theta)*s(psi)         c(theta)*s(phi)*s(psi)-s(theta)*c(psi)      c(theta)*c(phi)];

a_I = [1; 0; 0];

a_S = Rotate_I_S(a_I, theta, phi, psi);

% a_S


% %% Drawing Unit Vectors
% 
% figure
% 
% for i = 1:length(t)
% 
%     t_i = t(i);
% 
%     x = a.real_pos_xI(t_i);
%     y = a.real_pos_yI(t_i);
%     z = a.real_pos_zI(t_i);
% 
%     scatter3(x, y, z, 10, 'filled', "MarkerEdgeColor", "black", "MarkerFaceColor","black")
%     hold on
% 
%     scale = 10;
% 
%     % Draw intertial unit vectors
%     a_I = [1; 1; 1];
% 
%     a_vec_I = a_I/norm(a_I)/scale; 
% 
%     frame_vec_I = [x y z] + a_vec_I';
% 
%     plot3([x frame_vec_I(1)], [y y], [z z], Color='black')
%     hold on
% 
%     plot3([x x], [y frame_vec_I(2)], [z z], Color='black')
%     hold on
% 
%     plot3([x x], [y y], [z frame_vec_I(3)], Color='black')
%     hold on
% 
%     % Convert acceleration vectors from inertial frame to sensor frame
%     theta = a.theta(t_i);
%     phi = a.phi(t_i);
%     psi = a.psi(t_i);
% 
%     a_S = Rotate_I_S(a_I, theta, phi, psi);
% 
%     % Drawn acceleration vectors in sensor frame
% 
%     eastChange = [a_S(1)*cos(psi) a_S(1)*sin(psi) a_S(1)*sin(theta)];
%     northChange = [-a_S(2)*sin(psi) a_S(2)*cos(psi) a_S(2)*sin(phi)];
%     upChange = [-a_S(3)*sin(theta) -a_S(3)*sin(phi) a_S(3)*cos(phi)];
% 
%     xChange_I = eastChange/norm(eastChange)/scale;
%     yChange_I = northChange/norm(northChange)/scale;
%     zChange_I = upChange/norm(upChange)/scale;
%     sensorEast = [x y z] + xChange_I;
%     sensorNorth = [x y z] + yChange_I;
%     sensorUp = [x y z] + zChange_I;
% 
%     plot3([x sensorEast(1)], [y sensorEast(2)], [z sensorEast(3)], Color='magenta')
%     hold on
% 
%     plot3([x sensorNorth(1)], [y sensorNorth(2)], [z sensorNorth(3)], Color='magenta')
%     hold on
% 
%     plot3([x sensorUp(1)], [y sensorUp(2)], [z sensorUp(3)], Color='magenta')
% 
%     xlim([-1 1])
%     ylim([0 2.5])
%     zlim([0.5 1.5])
%     xlabel('East')
%     ylabel('North')
%     zlabel('Up')
%     drawnow
%     hold on
%     pause(dt)
% end
% 
% %% Drawing Only Sensor Frame
% 
% figure
% 
% for i = 1:length(t)
% 
%     t_i = t(i);
% 
%     x = a.real_pos_xI(t_i);
%     y = a.real_pos_yI(t_i);
%     z = a.real_pos_zI(t_i);
% 
%     scatter3(x, y, z, 10, 'filled', "MarkerEdgeColor", "black", "MarkerFaceColor","black")
%     hold on
% 
%     scale = 10;
% 
%     % Draw intertial unit vectors
%     a_I = [1; 1; 1];
% 
%     a_vec_I = a_I/norm(a_I)/scale; 
% 
%     frame_vec_I = [x y z] + a_vec_I';
% 
%     % Convert acceleration vectors from inertial frame to sensor frame
%     theta = a.theta(t_i);
%     phi = a.phi(t_i);
%     psi = a.psi(t_i);
% 
%     a_S = Rotate_I_S(a_I, theta, phi, psi);
% 
%     % Drawn acceleration vectors in sensor frame
% 
%     eastChange = [a_S(1)*cos(psi) a_S(1)*sin(psi) a_S(1)*sin(theta)];
%     northChange = [-a_S(2)*sin(psi) a_S(2)*cos(psi) a_S(2)*sin(phi)];
%     upChange = [-a_S(3)*sin(theta) -a_S(3)*sin(phi) a_S(3)*cos(phi)];
% 
%     xChange_I = eastChange/norm(eastChange)/scale;
%     yChange_I = northChange/norm(northChange)/scale;
%     zChange_I = upChange/norm(upChange)/scale;
%     sensorEast = [x y z] + xChange_I;
%     sensorNorth = [x y z] + yChange_I;
%     sensorUp = [x y z] + zChange_I;
% 
%     plot3([x sensorEast(1)], [y sensorEast(2)], [z sensorEast(3)], Color='blue')
%     hold on
% 
%     plot3([x sensorNorth(1)], [y sensorNorth(2)], [z sensorNorth(3)], Color='green')
%     hold on
% 
%     plot3([x sensorUp(1)], [y sensorUp(2)], [z sensorUp(3)], Color='magenta')
% 
%     xlim([-1 1])
%     ylim([0 2.5])
%     zlim([0.5 1.5])
%     xlabel('East')
%     ylabel('North')
%     zlabel('Up')
%     drawnow
%     hold on
%     pause(dt)
% end
% 
% 
% 
% %% Rotate from Inertial to Sensor Frame
% 
% figure
% 
% for i = 1:length(t)
% 
%     t_i = t(i);
% 
%     x = a.real_pos_xI(t_i);
%     y = a.real_pos_yI(t_i);
%     z = a.real_pos_zI(t_i);
% 
%     scatter3(x, y, z, 10, 'filled', "MarkerEdgeColor", "black", "MarkerFaceColor","black")
%     hold on
% 
%     scale = 10;
% 
%     % Draw intertial unit vectors
%     a_I = [a.real_accel_xI(t_i); a.real_accel_yI(t_i); a.real_accel_zI(t_i)];
% 
%     a_vec_I = a_I/norm(a_I)/scale;
% 
%     frame_vec_I = [x y z] + a_vec_I';
% 
%     % plot3([x frame_vec_I(1)], [y y], [z z], Color='black')
%     % hold on
%     % 
%     % plot3([x x], [y frame_vec_I(2)], [z z], Color='black')
%     % hold on
%     % 
%     % plot3([x x], [y y], [z frame_vec_I(3)], Color='black')
%     % hold on
% 
%     % Convert acceleration vectors from inertial frame to sensor frame
%     theta = a.theta(t_i);
%     phi = a.phi(t_i);
%     psi = a.psi(t_i);
% 
%     a_S = Rotate_I_S(a_I, theta, phi, psi);
% 
%     % Drawn acceleration vectors in sensor frame
% 
%     eastChange = [a_S(1)*cos(psi) a_S(1)*sin(psi) a_S(1)*sin(theta)];
%     northChange = [-a_S(2)*sin(psi) a_S(2)*cos(psi) a_S(2)*sin(phi)];
%     upChange = [-a_S(3)*sin(theta) -a_S(3)*sin(phi) a_S(3)*cos(phi)];
% 
%     xChange_I = eastChange/norm(eastChange)/scale;
%     yChange_I = northChange/norm(northChange)/scale;
%     zChange_I = upChange/norm(upChange)/scale;
%     sensorEast = [x y z] + xChange_I;
%     sensorNorth = [x y z] + yChange_I;
%     sensorUp = [x y z] + zChange_I;
% 
%     plot3([x sensorEast(1)], [y sensorEast(2)], [z sensorEast(3)], Color='magenta')
%     hold on
% 
%     plot3([x sensorNorth(1)], [y sensorNorth(2)], [z sensorNorth(3)], Color='magenta')
%     hold on
% 
%     plot3([x sensorUp(1)], [y sensorUp(2)], [z sensorUp(3)], Color='magenta')
% 
%     xlim([-1 1])
%     ylim([0 2.5])
%     zlim([0.5 1.5])
%     xlabel('East')
%     ylabel('North')
%     zlabel('Up')
%     drawnow
%     hold on
%     pause(dt)
% end


%% Draw Inertial Vectors in Planar View

scale = 10;

figure

subplot(3, 1, 1)

for i = 1:length(t)

    t_i = t(i);

    x = a.real_pos_xI(t_i);
    y = a.real_pos_yI(t_i);
    z = a.real_pos_zI(t_i);

    % Draw intertial acceleration vectors
    a_I = [a.real_accel_xI(t_i); a.real_accel_yI(t_i); a.real_accel_zI(t_i)];

    a_vec_I = a_I/norm(a_I)/scale; 

    frame_vec_I = [x y z] + a_vec_I';
    
    % Plotting only east and up inertial vectors
    plot([x frame_vec_I(1)], [z z], Color='red')
    hold on

    plot([x x], [z frame_vec_I(3)], Color='black')
    hold on

end

title("Platform/Deck Position in each Plane (Inertial Frame)")
xlabel('East (m)')
ylabel('Up (m)')
legend('East', 'Up')
plot(a.real_pos_xI(t), a.real_pos_zI(t), HandleVisibility="off", Color='#0072BD')

subplot(3, 1, 2)

for i = 1:length(t)

    t_i = t(i);

    x = a.real_pos_xI(t_i);
    y = a.real_pos_yI(t_i);
    z = a.real_pos_zI(t_i);
   
    % Draw intertial acceleration vectors
    a_I = [a.real_accel_xI(t_i); a.real_accel_yI(t_i); a.real_accel_zI(t_i)];

    a_vec_I = a_I/norm(a_I)/scale; 

    frame_vec_I = [x y z] + a_vec_I';
    
    % Plotting only north and up inertial vectors
    plot([y frame_vec_I(2)], [z z], Color='#EDB120')
    hold on

    plot([y y], [z frame_vec_I(3)], Color='black')
    hold on

end


xlabel('North (m)')
ylabel('Up (m)')
legend('North', 'Up')
plot(a.real_pos_yI(t), a.real_pos_zI(t), HandleVisibility="off", Color='#0072BD')

subplot(3, 1, 3)

for i = 1:length(t)

    t_i = t(i);

    x = a.real_pos_xI(t_i);
    y = a.real_pos_yI(t_i);
    z = a.real_pos_zI(t_i);

    % Draw intertial acceleration vectors
    a_I = [a.real_accel_xI(t_i); a.real_accel_yI(t_i); a.real_accel_zI(t_i)];

    a_vec_I = a_I/norm(a_I)/scale; 

    frame_vec_I = [x y z] + a_vec_I';
    
    % Plotting only east and north inertial vectors
    plot([x frame_vec_I(1)], [y y], Color='red')
    hold on

    plot([x x], [y frame_vec_I(2)], Color='#EDB120')
    hold on

end

xlabel('East (m)')
ylabel('North (m)')
legend('East', 'North')
plot(a.real_pos_xI(t), a.real_pos_yI(t), HandleVisibility="off", Color='#0072BD')

%% Draw Rotated Vector in Planar View

scale = 20;

figure

subplot(3, 1, 1)

for i = 1:length(t)

    t_i = t(i);

    x = a.real_pos_xI(t_i);
    y = a.real_pos_yI(t_i);
    z = a.real_pos_zI(t_i);

    % Draw intertial acceleration vectors
    a_I = [a.real_accel_xI(t_i); a.real_accel_yI(t_i); a.real_accel_zI(t_i)];

    a_vec_I = a_I/norm(a_I)/scale; 

    frame_vec_I = [x y z] + a_vec_I';
    
    % Convert acceleration vectors from inertial frame to sensor frame
    theta = a.theta(t_i);
    phi = a.phi(t_i);
    psi = a.psi(t_i);

    a_S = Rotate_I_S(a_I, theta, phi, psi);
    
    % Drawn acceleration vectors in sensor frame
    
    eastChange = [a_S(1)*cos(psi) a_S(1)*sin(psi) a_S(1)*sin(theta)];
    northChange = [-a_S(2)*sin(psi) a_S(2)*cos(psi) a_S(2)*sin(phi)];
    upChange = [-a_S(3)*sin(theta) -a_S(3)*sin(phi) a_S(3)*cos(phi)];

    xChange_I = eastChange/norm(eastChange)/scale;
    yChange_I = northChange/norm(northChange)/scale;
    zChange_I = upChange/norm(upChange)/scale;
    sensorEast = [x y z] + xChange_I;
    sensorNorth = [x y z] + yChange_I;
    sensorUp = [x y z] + zChange_I;
    
    % Plotting only east and up sensor unit vectors
    plot([x sensorEast(1)], [z sensorEast(3)], Color='blue')
    hold on

    plot([x sensorUp(1)], [z sensorUp(3)], Color='magenta')
    hold on
    % 
    % xlim([-1 1])
    % zlim([0.5 1.5])

end

title("Platform/Deck Position in each Plane (Sensor Frame)")
xlabel('East (m)')
ylabel('Up (m)')
legend('East', 'Up')
plot(a.real_pos_xI(t), a.real_pos_zI(t), HandleVisibility="off")

subplot(3, 1, 2)

for i = 1:length(t)

    t_i = t(i);

    x = a.real_pos_xI(t_i);
    y = a.real_pos_yI(t_i);
    z = a.real_pos_zI(t_i);
   
    % Draw intertial acceleration vectors
    a_I = [a.real_accel_xI(t_i); a.real_accel_yI(t_i); a.real_accel_zI(t_i)];

    a_vec_I = a_I/norm(a_I)/scale; 

    frame_vec_I = [x y z] + a_vec_I';
    
    % Convert acceleration vectors from inertial frame to sensor frame
    theta = a.theta(t_i);
    phi = a.phi(t_i);
    psi = a.psi(t_i);

    a_S = Rotate_I_S(a_I, theta, phi, psi);
    
    % Drawn acceleration vectors in sensor frame
    
    eastChange = [a_S(1)*cos(psi) a_S(1)*sin(psi) a_S(1)*sin(theta)];
    northChange = [-a_S(2)*sin(psi) a_S(2)*cos(psi) a_S(2)*sin(phi)];
    upChange = [-a_S(3)*sin(theta) -a_S(3)*sin(phi) a_S(3)*cos(phi)];

    xChange_I = eastChange/norm(eastChange)/scale;
    yChange_I = northChange/norm(northChange)/scale;
    zChange_I = upChange/norm(upChange)/scale;
    sensorEast = [x y z] + xChange_I;
    sensorNorth = [x y z] + yChange_I;
    sensorUp = [x y z] + zChange_I;
    
    % Plotting only east and up sensor unit vectors
    plot([y sensorNorth(2)], [z sensorNorth(3)], Color='green')
    hold on

    plot([y sensorUp(2)], [z sensorUp(3)], Color='magenta')
    hold on

    % xlim([-1 1])
    % zlim([0.5 1.5])

end

xlabel('North (m)')
ylabel('Up (m)')
legend('North', 'Up')
plot(a.real_pos_yI(t), a.real_pos_zI(t), HandleVisibility="off")

subplot(3, 1, 3)

for i = 1:length(t)

    t_i = t(i);

    x = a.real_pos_xI(t_i);
    y = a.real_pos_yI(t_i);
    z = a.real_pos_zI(t_i);

    % Draw intertial acceleration vectors
    a_I = [a.real_accel_xI(t_i); a.real_accel_yI(t_i); a.real_accel_zI(t_i)];

    a_vec_I = a_I/norm(a_I)/scale; 

    frame_vec_I = [x y z] + a_vec_I';
    
    % Convert acceleration vectors from inertial frame to sensor frame
    theta = a.theta(t_i);
    phi = a.phi(t_i);
    psi = a.psi(t_i);

    a_S = Rotate_I_S(a_I, theta, phi, psi);
    
    % Drawn acceleration vectors in sensor frame
    
    eastChange = [a_S(1)*cos(psi) a_S(1)*sin(psi) a_S(1)*sin(theta)];
    northChange = [-a_S(2)*sin(psi) a_S(2)*cos(psi) a_S(2)*sin(phi)];
    upChange = [-a_S(3)*sin(theta) -a_S(3)*sin(phi) a_S(3)*cos(phi)];

    xChange_I = eastChange/norm(eastChange)/scale;
    yChange_I = northChange/norm(northChange)/scale;
    zChange_I = upChange/norm(upChange)/scale;
    sensorEast = [x y z] + xChange_I;
    sensorNorth = [x y z] + yChange_I;
    sensorUp = [x y z] + zChange_I;
    
    % Plotting only east and up sensor unit vectors
    plot([x sensorEast(1)], [y sensorEast(2)], Color='blue')
    hold on

    plot([x sensorNorth(1)], [y sensorNorth(2)], Color='green')
    hold on

    % xlim([-1 1])
    % zlim([0.5 1.5])

end

xlabel('East (m)')
ylabel('North (m)')
legend('East', 'North')
plot(a.real_pos_xI(t), a.real_pos_yI(t), HandleVisibility="off")

% %% Draw Over Time
% 
% scale = 500;
% 
% 
% figure
% 
% subplot(3, 1, 1)
% plot(t, a.real_pos_zI(t))
% title("Platform/Deck Position over Time")
% hold on
% 
% for i = 1:length(t)
% 
%     t_i = t(i);
% 
%     x = a.real_pos_xI(t_i);
%     y = a.real_pos_yI(t_i);
%     z = a.real_pos_zI(t_i);
% 
%     % Draw intertial acceleration vectors
%     a_I = [a.real_accel_xI(t_i); a.real_accel_yI(t_i); a.real_accel_zI(t_i)];
% 
%     a_vec_I = a_I/norm(a_I)/scale; 
% 
%     frame_vec_I = [x y z] + a_vec_I';
% 
%     % Convert acceleration vectors from inertial frame to sensor frame
%     theta = a.theta(t_i);
%     phi = a.phi(t_i);
%     psi = a.psi(t_i);
% 
%     a_S = Rotate_I_S(a_I, theta, phi, psi);
% 
%     % Drawn acceleration vectors in sensor frame
% 
%     eastChange = [a_S(1)*cos(psi) a_S(1)*sin(psi) a_S(1)*sin(theta)];
%     northChange = [-a_S(2)*sin(psi) a_S(2)*cos(psi) a_S(2)*sin(phi)];
%     upChange = [-a_S(3)*sin(theta) -a_S(3)*sin(phi) a_S(3)*cos(phi)];
% 
%     xChange_I = eastChange/norm(eastChange)/scale;
%     yChange_I = northChange/norm(northChange)/scale;
%     zChange_I = upChange/norm(upChange)/scale;
% 
%     % sensorEast = [x y z] + xChange_I;
%     % sensorNorth = [x y z] + yChange_I;
%     % sensorUp = [x y z] + zChange_I;
% 
%     sensorEast = xChange_I;
%     sensorNorth = yChange_I;
%     sensorUp = zChange_I;
% 
%     % Plotting only east and up sensor unit vectors
%     plot([t t+sensorEast(1)], [z z+sensorEast(3)], Color='blue')
%     hold on
% 
%     plot([t t+sensorUp(1)], [z z+sensorUp(3)], Color='magenta')
%     hold on
%     % 
%     % xlim([-1 1])
%     % zlim([0.5 1.5])
% 
% end
% 
% xlabel('Time (sec)')
% ylabel('Up (m)')
% 
% subplot(3, 1, 2)
% plot(t, a.real_pos_yI(t))
% 
% for i = 1:length(t)
% 
%     t_i = t(i);
% 
%     x = a.real_pos_xI(t_i);
%     y = a.real_pos_yI(t_i);
%     z = a.real_pos_zI(t_i);
% 
%     % Draw intertial acceleration vectors
%     a_I = [a.real_accel_xI(t_i); a.real_accel_yI(t_i); a.real_accel_zI(t_i)];
% 
%     a_vec_I = a_I/norm(a_I)/scale; 
% 
%     frame_vec_I = [x y z] + a_vec_I';
% 
%     % Convert acceleration vectors from inertial frame to sensor frame
%     theta = a.theta(t_i);
%     phi = a.phi(t_i);
%     psi = a.psi(t_i);
% 
%     a_S = Rotate_I_S(a_I, theta, phi, psi);
% 
%     % Drawn acceleration vectors in sensor frame
% 
%     eastChange = [a_S(1)*cos(psi) a_S(1)*sin(psi) a_S(1)*sin(theta)];
%     northChange = [-a_S(2)*sin(psi) a_S(2)*cos(psi) a_S(2)*sin(phi)];
%     upChange = [-a_S(3)*sin(theta) -a_S(3)*sin(phi) a_S(3)*cos(phi)];
% 
%     xChange_I = eastChange/norm(eastChange)/scale;
%     yChange_I = northChange/norm(northChange)/scale;
%     zChange_I = upChange/norm(upChange)/scale;
%     sensorEast = [x y z] + xChange_I;
%     sensorNorth = [x y z] + yChange_I;
%     sensorUp = [x y z] + zChange_I;
% 
%     % Plotting only east and up sensor unit vectors
%     plot([t t], [sensorNorth(2) sensorNorth(3)], Color='green')
%     hold on
% 
%     plot([t t], [sensorUp(2) sensorUp(3)], Color='magenta')
%     hold on
% 
%     % xlim([-1 1])
%     % zlim([0.5 1.5])
% 
% end
% 
% xlabel('Time (sec)')
% ylabel('North (m)')
% 
% subplot(3, 1, 3)
% plot(t, a.real_pos_xI(t))
% 
% 
% for i = 1:length(t)
% 
%     t_i = t(i);
% 
%     x = a.real_pos_xI(t_i);
%     y = a.real_pos_yI(t_i);
%     z = a.real_pos_zI(t_i);
% 
%     % Draw intertial acceleration vectors
%     a_I = [a.real_accel_xI(t_i); a.real_accel_yI(t_i); a.real_accel_zI(t_i)];
% 
%     a_vec_I = a_I/norm(a_I)/scale; 
% 
%     frame_vec_I = [x y z] + a_vec_I';
% 
%     % Plotting only east and up inertial vectors
%     plot([x frame_vec_I(1)], [y y], Color='black')
%     hold on
% 
%     plot([x x], [y frame_vec_I(2)], Color='black')
%     hold on
% 
%     % Convert acceleration vectors from inertial frame to sensor frame
%     theta = a.theta(t_i);
%     phi = a.phi(t_i);
%     psi = a.psi(t_i);
% 
%     a_S = Rotate_I_S(a_I, theta, phi, psi);
% 
%     % Drawn acceleration vectors in sensor frame
% 
%     eastChange = [a_S(1)*cos(psi) a_S(1)*sin(psi) a_S(1)*sin(theta)];
%     northChange = [-a_S(2)*sin(psi) a_S(2)*cos(psi) a_S(2)*sin(phi)];
%     upChange = [-a_S(3)*sin(theta) -a_S(3)*sin(phi) a_S(3)*cos(phi)];
% 
%     xChange_I = eastChange/norm(eastChange)/scale;
%     yChange_I = northChange/norm(northChange)/scale;
%     zChange_I = upChange/norm(upChange)/scale;
%     sensorEast = [x y z] + xChange_I;
%     sensorNorth = [x y z] + yChange_I;
%     sensorUp = [x y z] + zChange_I;
% 
%     % Plotting only east and up sensor unit vectors
%     plot([t t], [sensorEast(1) sensorEast(2)], Color='blue')
%     hold on
% 
%     plot([t t], [sensorNorth(1) sensorNorth(2)], Color='green')
%     hold on
% 
%     % xlim([-1 1])
%     % zlim([0.5 1.5])
% 
% end
% 
% xlabel('Time (sec)')
% ylabel('East (m)')
