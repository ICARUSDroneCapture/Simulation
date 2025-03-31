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
a.real_accel_zI = @(t) -beta^2*alpha*sin(beta*t) - 9.81; % [m*s^-2]

a.theta = @(t) -atan(beta*alpha*cos(beta*t)); % [rad]
a.phi = @(t) atan(0.2/beta*cos(beta/2*t)); % [rad]
a.psi = @(t) atan(0.4/beta*cos(beta/4*t)); % [rad]

% a.theta_dot = @(t) ((alpha*beta^2*sin(beta*t))./(alpha^2*beta^2*(cos(beta*t).^2)+1)); % [rad/s]
% a.phi_dot = @(t) ((-0.1*sin(beta/2*t))/(((0.04*(cos(beta*t/2).^2))/(beta^2))+1)); % [rad/s]
% a.psi_dot = @(t) ((-0.1*sin(beta/4*t))/(((0.16*(cos(beta*t/4).^2))/(beta^2))+1)); % [rad/s]

s = @(x) sin(x);
c = @(x) cos(x);

%% Sensor Model Aspects

% Simulation time
startTime = 0;
finishTime = 10;
tspan = [startTime finishTime]; % [s]

% dt = 1/imu_rate;  % [s]
dt = 0.001;
t = (tspan(1):dt:tspan(2))';
t_count = length(t);
indeces = @(t) floor(t/dt)+1;

defineSignals
% defineSignalsNoNoise

%% Take derivatives of angle functions manually


% Take derivatives of angle functions

theta_vals = a.theta(t);
phi_vals = a.phi(t);
psi_vals = a.psi(t);

theta_dot = zeros(1, length(t));
phi_dot = zeros(1, length(t));
psi_dot = zeros(1, length(t));

theta_dot(2:end) = diff(theta_vals)/dt;
phi_dot(2:end) = diff(phi_vals)/dt;
psi_dot(2:end) = diff(psi_vals)/dt;

figure

subplot(3, 1, 1)
plot(t, theta_dot/pi*180)
title("Platform/Deck Angular Rate (theta)")
xlabel('Time (sec)')
ylabel('Angular Rate (dps)')

subplot(3, 1, 2)
plot(t, phi_dot/pi*180)
title("Platform/Deck Angular Rate (phi)")
xlabel('Time (sec)')
ylabel('Angular Rate (dps)')

subplot(3, 1, 3)
plot(t, psi_dot/pi*180)
title("Platform/Deck Angular Rate (psi)")
xlabel('Time (sec)')
ylabel('Angular Rate (dps)')


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

% %% Draw Rotated Vector in Planar View
% 
% scale = 20;
% 
% figure
% 
% subplot(3, 1, 1)
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
%     plot([x sensorEast(1)], [z sensorEast(3)], Color='blue')
%     hold on
% 
%     plot([x sensorUp(1)], [z sensorUp(3)], Color='magenta')
%     hold on
%     % 
%     % xlim([-1 1])
%     % zlim([0.5 1.5])
% 
% end
% 
% title("Platform/Deck Position in each Plane (Sensor Frame)")
% xlabel('East (m)')
% ylabel('Up (m)')
% legend('East', 'Up')
% plot(a.real_pos_xI(t), a.real_pos_zI(t), HandleVisibility="off")
% 
% subplot(3, 1, 2)
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
%     plot([y sensorNorth(2)], [z sensorNorth(3)], Color='green')
%     hold on
% 
%     plot([y sensorUp(2)], [z sensorUp(3)], Color='magenta')
%     hold on
% 
%     % xlim([-1 1])
%     % zlim([0.5 1.5])
% 
% end
% 
% xlabel('North (m)')
% ylabel('Up (m)')
% legend('North', 'Up')
% plot(a.real_pos_yI(t), a.real_pos_zI(t), HandleVisibility="off")
% 
% subplot(3, 1, 3)
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
%     plot([x sensorEast(1)], [y sensorEast(2)], Color='blue')
%     hold on
% 
%     plot([x sensorNorth(1)], [y sensorNorth(2)], Color='green')
%     hold on
% 
%     % xlim([-1 1])
%     % zlim([0.5 1.5])
% 
% end
% 
% xlabel('East (m)')
% ylabel('North (m)')
% legend('East', 'North')
% plot(a.real_pos_xI(t), a.real_pos_yI(t), HandleVisibility="off")


%% Get Sensor Values over Time (without sensor error)

a_S_real = zeros(3, length(t));

for i = 1:length(t)

    t_i = t(i);

    theta = a.theta(t_i);
    phi = a.phi(t_i);
    psi = a.psi(t_i);

    a_I = [a.real_accel_xI(t_i); a.real_accel_yI(t_i); a.real_accel_zI(t_i)];

    a_S = Rotate_I_S(a_I, theta, phi, psi);

    a_S_real(:, i) = a_S;

end

%% Integrate Angular Rate with Sensor Error


% Initial Angle
p_theta0 = a.theta(tspan(1)); % Platform inertial angle [deg]
p_phi0 = a.phi(tspan(1)); % Platform inertial angle [deg]
p_psi0 = a.psi(tspan(1)); % Platform inertial angle [deg]

g_S_measured = zeros(3, length(t));

for i = 1:length(t)

    t_i = t(i);

    theta_d = theta_dot(i);
    phi_d = phi_dot(i);
    psi_d = psi_dot(i);

    theta_dot_m = a.measured_gyro(t_i, a.o_d_n_g_c, a.biasStabDistGyro, a.biasTempDistGyro, a.gyro_drift, a.noiseDistGyro, theta_d);
    phi_dot_m = a.measured_gyro(t_i, a.o_d_n_g_c, a.biasStabDistGyro, a.biasTempDistGyro, a.gyro_drift, a.noiseDistGyro, phi_d);
    psi_dot_m = a.measured_gyro(t_i, a.o_d_n_g_c, a.biasStabDistGyro, a.biasTempDistGyro, a.gyro_drift, a.noiseDistGyro, psi_d);

    g_S_measured(:, i) = [theta_dot_m; phi_dot_m; psi_dot_m];

end


a.theta_dot_eq = @(t, y) g_S_measured(1, floor(t./dt)+1);
a.phi_dot_eq = @(t, y) g_S_measured(2, floor(t./dt)+1);
a.psi_dot_eq = @(t, y) g_S_measured(3, floor(t./dt)+1);

[t, theta_err]= rk4_solver(a.theta_dot_eq, tspan, p_theta0, dt);
[t, phi_err]= rk4_solver(a.phi_dot_eq, tspan, p_phi0, dt);
[t, psi_err]= rk4_solver(a.psi_dot_eq, tspan, p_psi0, dt);

% Plot integrated angle from angular rate with error
figure
subplot(3,1,1)
plot(t, theta_err / pi * 180)
hold on
plot(t, a.theta(t) / pi * 180)
xlabel('Time (sec)')
ylabel('Angle (deg)')
title('Platform Angle')
legend('Measured', 'Real')

subplot(3,1,2)
plot(t, phi_err / pi * 180)
hold on
plot(t, a.phi(t) / pi * 180)
xlabel('Time (sec)')
ylabel('Angle (deg)')
title('Platform Angle')
legend('Measured', 'Real')

subplot(3,1,3)
plot(t, psi_err / pi * 180)
hold on
plot(t, a.psi(t) / pi * 180)
xlabel('Time (sec)')
ylabel('Angle (deg)')
title('Platform Angle')
legend('Measured', 'Real')

%% Integrate Angular Rate with Sensor Error

a_S_measured = zeros(3, length(t));

for i = 1:length(t)

    t_i = t(i);

    theta = theta_err(i);
    phi = phi_err(i);
    psi = psi_err(i);

    a_I = [a.real_accel_xI(t_i); a.real_accel_yI(t_i); a.real_accel_zI(t_i)];

    a_S = Rotate_I_S(a_I, theta, phi, psi);

    a_S_measured(:, i) = a_S;

end


%% Get Sensor Values over Time (with sensor error)

figure
subplot(3,1,1)
plot(t, g_S_measured(1, :))
hold on
plot(t, theta_dot)
xlabel('Time (sec)')
ylabel('Angular Rate (dps)')
title('Angular Rate (theta)')
legend('Measured', 'Real')

subplot(3,1,2)
plot(t, g_S_measured(2, :))
hold on
plot(t, phi_dot)
xlabel('Time (sec)')
ylabel('Angular Rate (dps)')
title('Angular Rate (phi)')
legend('Measured', 'Real')

subplot(3,1,3)
plot(t, g_S_measured(3, :))
hold on
plot(t, psi_dot)
xlabel('Time (sec)')
ylabel('Angular Rate (dps)')
title('Angular Rate (psi)')
legend('Measured', 'Real')

%% Plot acceleration difference

figure
subplot(3,1,1)
plot(t, a_S_measured(1, :))
hold on
plot(t, a_S_real(1, :))
xlabel('Time (sec)')
ylabel('Acceleration (m/s^2)')
title('Sensor Frame Acceleration X (m/s^2)')
legend('Measured', 'Real')

subplot(3,1,2)
plot(t, a_S_measured(2, :))
hold on
plot(t, a_S_real(2, :))
xlabel('Time (sec)')
ylabel('Acceleration (m/s^2)')
title('Sensor Frame Acceleration Y (m/s^2)')
legend('Measured', 'Real')

subplot(3,1,3)
plot(t, a_S_measured(3, :))
hold on
plot(t, a_S_real(3, :))
xlabel('Time (sec)')
ylabel('Acceleration (m/s^2)')
title('Sensor Frame Acceleration Z (m/s^2)')
legend('Measured', 'Real')

%% Rotate from Sensor back to Inertial

% Reintegrated for angle values with new random error numbers
g_S_measured = zeros(3, length(t));

for i = 1:length(t)

    t_i = t(i);

    theta_d = theta_dot(i);
    phi_d = phi_dot(i);
    psi_d = psi_dot(i);

    theta_dot_m = a.measured_gyro(t_i, a.o_d_n_g_c, a.biasStabDistGyro, a.biasTempDistGyro, a.gyro_drift, a.noiseDistGyro, theta_d);
    phi_dot_m = a.measured_gyro(t_i, a.o_d_n_g_c, a.biasStabDistGyro, a.biasTempDistGyro, a.gyro_drift, a.noiseDistGyro, phi_d);
    psi_dot_m = a.measured_gyro(t_i, a.o_d_n_g_c, a.biasStabDistGyro, a.biasTempDistGyro, a.gyro_drift, a.noiseDistGyro, psi_d);

    g_S_measured(:, i) = [theta_dot_m; phi_dot_m; psi_dot_m];

end


a.theta_dot_eq = @(t, y) g_S_measured(1, floor(t./dt)+1);
a.phi_dot_eq = @(t, y) g_S_measured(2, floor(t./dt)+1);
a.psi_dot_eq = @(t, y) g_S_measured(3, floor(t./dt)+1);

[t, theta_err]= rk4_solver(a.theta_dot_eq, tspan, p_theta0, dt);
[t, phi_err]= rk4_solver(a.phi_dot_eq, tspan, p_phi0, dt);
[t, psi_err]= rk4_solver(a.psi_dot_eq, tspan, p_psi0, dt);

a_I_time = zeros(3, length(t));

for i = 1:length(t)

    a_S = a_S_measured(:, i);

    t_i = t(i);

    theta = theta_err(i);
    phi = phi_err(i);
    psi = psi_err(i);

    a_I = Rotate_S_I(a_S, theta, phi, psi);

    a_I_time(:, i) = a_I;

end


figure
subplot(3,1,1)
plot(t, a_I_time(1, :))
hold on
plot(t, a.real_accel_xI(t))
xlabel('Time (sec)')
ylabel('Acceleration (m/s^2)')
title('Rotated Inertial Acceleration X (m/s^2)')
legend('Measured', 'Real')

subplot(3,1,2)
plot(t, a_I_time(2, :))
hold on
plot(t, a.real_accel_yI(t))
xlabel('Time (sec)')
ylabel('Acceleration (m/s^2)')
title('Rotated Inertial Acceleration Y (m/s^2)')
legend('Measured', 'Real')

subplot(3,1,3)
plot(t, a_I_time(3, :))
hold on
plot(t, a.real_accel_zI(t))
xlabel('Time (sec)')
ylabel('Acceleration (m/s^2)')
title('Rotated Inertial Acceleration Z (m/s^2)')
legend('Measured', 'Real')

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
    a_I = a_I_time(:, i);

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
    a_I = a_I_time(:, i);

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


%% Check Difference


a_I_og = [a.real_accel_xI(t); a.real_accel_yI(t); a.real_accel_zI(t)];
a_I_og(3, :) = a_I_og(3, :) + 9.81;

a_diff = a_I_time - a_I_og;

figure
scatter(t, a_diff)
xlabel('Time (sec)')
ylabel('Acceleration (m/s^2)')
title('Numerical Difference in Acceleration (Post Rotations)')

