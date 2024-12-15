%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%% Arm Parameters %%%%%%%%%%%%%%%%%%%%%%%

a.m = 1;    % Mass [kg]
a.g = 9.81; % Acceleration of gravity [m/s^2]
a.pr_d = 0.5; % Desired relative position [m]

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%% Environmental Model %%%%%%%%%%%%%%%%%%%%

alpha = 0.4; % wave amplitdue [m]
hdeck = 1;   % inertial reference deck hight [m] (arbitrary)

% Wave frequency
Tmax = 7.5;    % Maximum period [s]
k = 1;
T = Tmax / k;  % Period of deck disturbance [s]
beta = 2*pi/T; % wave frequency [rad/s]

% Inertial Position, Velocity, and Acceleration of Deck

% Sine Wave
a.d = @(t) alpha*sin(beta*t) + hdeck;      % [m]
a.d_dot = @(t) beta*alpha*cos(beta*t);    % [m/s]
a.d_ddot = @(t) -beta^2*alpha*sin(beta*t); % [m*s^-2]
 
% Cosine Wave
% a.d = @(t) alpha*cos(beta*t) + hdeck;      % [m]
% a.d_dot = @(t) -beta*alpha*sin(beta*t);    % [m/s]
% a.d_ddot = @(t) -beta^2*alpha*cos(beta*t); % [m*s^-2]

% Square Wave
% N = 3; % Number of terms in the Fourier series
% coefficients = 1:2:(2*N - 1); % Odd harmonics: 1, 3, 5, ..., (2*N-1)
% beta2 = 2*pi/15;
% a.sw = @(t) alpha*((4/pi) * ...
%     sum(arrayfun(@(n) sin(n*beta2*t)/n, coefficients))) + hdeck; % [m]
% a.d = @(t) arrayfun(a.sw, t);
% a.sw_dot = @(t) alpha*((4/pi) * ...
%     sum(arrayfun(@(n) beta2*cos(n*beta2*t), coefficients)));    % [m/s]
% a.d_dot = @(t) arrayfun(a.sw_dot, t);
% a.sw_ddot = @(t) alpha*((4/pi) * ...
%     sum(arrayfun(@(n) -n*beta2^2*sin(n*beta2*t), coefficients))); % [m*s^-2]
% a.d_ddot = @(t) arrayfun(a.sw_ddot, t);

% Stacked sine wave
% a.d = @(t) alpha*(1.5*sin(beta*t/6) + 0.75*sin(beta*t)) + hdeck; % [m]
% a.d_dot = @(t) alpha*(0.25*beta*cos(beta*t/6) ...
%                                     + 0.75*beta*cos(beta*t)); % [m/s]
% a.d_ddot = @(t) -alpha*(0.0417*beta^2*sin(beta*t/6) ...
%                                     + 0.75*beta^2*sin(beta*t)); % [m*s^-2]

a.a = @(t) atan(beta*alpha*cos(beta*t)); % [rad]
a.a_dot = @(t) -(alpha*beta^2*sin(beta*t))/(alpha^2*beta^2*(cos(beta*t)^2)+1); % [rad]

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%% Control Gains %%%%%%%%%%%%%%%%%%%%%%%

% Minimum Required Control Constants:
%   a.ka = 2800; a.kv = ?; a.ks = ?;
%   a.kp = 1000; a.kd = 10; a.ki = 100;

% Inertial Stabilization Control
a.ka =  1.8;  % Acceleration Control [kg]
a.kv = 90;  % Velocity Control [kg/s]
a.ks = 0;  % Position Control [kg*s^-2] 

% Relative Position Control
a.kp = 7;  % Proportional [kg*s^-2]
a.kd = 50;  % Derivative [kg/s]    
a.ki = 2;  % Integral [kg*s^-3]

% Low-Pass filter on measured acceleration
a.f_c = 1000; % Cutoff frequency [Hz]
a.omega = 2*pi*a.f_c; % Angular frequency [rad/s]

% Progressively Increase inertial stability gains to full gains so initial
% large values of velocity and acceleration do not cause large control
% forces
s0 = 0;  % Inital proportion of gain values to apply
gain_rate = 0.5; % Rate at which gains are increased
a.initial_scale = @(t) s0 + (1 - s0) * (1 - exp(-gain_rate*t));

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%% Control Gain Mixing %%%%%%%%%%%%%%%%%%%%

% If |pr-pr_d| <= r_g, apply full inertial stability control
a.r_g = 0.41; % [m]

% If |pr-pr_d| <= r_k, apply "h_k" proportion of relative position control
a.r_k = 0.41; % [m]
a.h_k = 0.0; % Proportion of applied relative position control

% Proportion of applied inertial stabilitycontrol
n = 1; % Polynomial order
a.I = @(x) ...
    ((-1/(1-(a.pr_d+a.r_g))^n) * ...
        (sign(x-a.pr_d).*(x-(a.pr_d+sign(x-a.pr_d)*a.r_g))).^n + 1) .* ...
            (abs(x-a.pr_d)> a.r_g) ...
    + 1 * (abs(x-a.pr_d) <= a.r_g);

% Proportion of applied relative position control
n = 1; % Polynomial order
a.K = @(x) ...
    ((1/(1-(a.pr_d+a.r_k))^n) * ...
        (sign(x-a.pr_d).*(x-(a.pr_d+sign(x-a.pr_d)*a.r_k))).^n) .*...
            (abs(x-a.pr_d)> a.r_k);
a.K_h = @(x) ...
    (-((a.h_k-1)/(1-(a.pr_d+a.r_k))^n) * ...
        (sign(x-a.pr_d).*(x-(a.pr_d+sign(x-a.pr_d)*a.r_k))).^n + a.h_k) .*...
            (abs(x-a.pr_d)> a.r_k) ...
    + a.h_k * (abs(x-a.pr_d) <= a.r_k);

%%%               Gain Mixing Visualization               %%%

pr = linspace(0,1,1000); % relative position test values

figure;
plot(pr,a.I(pr)) % Proportion of inertial control gain
hold on
plot(pr, a.K(pr)) % Proportion of relative position control gain
plot(pr, a.K_h(pr)) % Proportion of non-zero relative position control gain
xline(a.pr_d,'--','Label','$p_{rd}$','Interpreter','latex','FontSize',15,...
    'LabelOrientation','horizontal','LabelVerticalAlignment','middle')
title('Inertial and Relative Positional Control Mixing')
xlabel('Relative Position (m)')
ylabel('Gain Proportion')
xlim([-0.1 1.1])
ylim([-0.1 1.1])
legend('Inertial','Relative','Non-zero Relative')
grid on

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%% Performance Parameters %%%%%%%%%%%%%%%%%%%

% Maximum acceleration metric
p_ddot_max = 0.005*beta^2;

% Settling time
t_s = 2;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%