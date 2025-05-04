%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%% Arm Parameters %%%%%%%%%%%%%%%%%%%%%%%

a.m = 3.08;    % Mass [kg]
a.g = 9.81; % Acceleration of gravity [m/s^2]
a.pr_d = [0; 0; 0.5]; % Desired relative position [m]

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%% Environmental Model %%%%%%%%%%%%%%%%%%%%


define_sim_environment


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%% Control Gains %%%%%%%%%%%%%%%%%%%%%%%

plot_gain = false;

% Minimum Required Control Constants:
%   a.ka = 2800; a.kv = ?; a.ks = ?;
%   a.kp = 1000; a.kd = 10; a.ki = 100;


% % ----------------------------- Non-Zero Gains ----------------------------

% % Inertial Stabilization Control
% a.ka = [2; 2; 8];  % Acceleration Control [kg]
% a.kv = [12; 12; 72];  % Velocity Control [kg]

% Relative Position Control at center
scale_pc = 3;
a.kp_c = scale_pc*[0.3; 0.3; 5];  % Proportional [kg*s^-2]
a.kd_c = scale_pc*[0.0000001; 0.0000001; 1];  % Derivative [kg/s]    
a.ki_c = scale_pc*[0.0000000001; 0.0000000001; 0.1];  % Integral [kg*s^-3]

% Relative Position Control at boundaries
scale_pb = 1;
a.kp_b = scale_pb*[5; 5; 50];  % Proportional [kg*s^-2]
a.kd_b = scale_pb*[1; 1; 30];  % Derivative [kg/s]    
a.ki_b = scale_pb*[0.05; 0.05; 0.5];  % Integral [kg*s^-3]
% 
% % -------------------------------------------------------------------------



% ------------------------------- Zero Gains ------------------------------

% Inertial Stabilization Control
a.ka = [0; 0; 0];  % Acceleration Control [kg]
a.kv = [0; 0; 0];  % Velocity Control [kg]

% % Relative Position Control at center
% scale = 1;
% a.kp_c = scale*[0; 0; 0];  % Proportional [kg*s^-2]
% a.kd_c = scale*[0; 0; 0];  % Derivative [kg/s]
% a.ki_c = scale*[0; 0; 0];  % Integral [kg*s^-3]
% 
% % Relative Position Control at boundaries
% a.kp_b = [0; 0; 0];  % Proportional [kg*s^-2]
% a.kd_b = [0; 0; 0];  % Derivative [kg/s]    
% a.ki_b = [0; 0; 0];  % Integral [kg*s^-3]

% -------------------------------------------------------------------------

% % Inertial Stabilization Control
% a.ka =  0.98;  % Acceleration Control [kg]
% a.kv = 2;  % Velocity Control [kg/s]
% a.ks = 0;  % Position Control [kg*s^-2]
% 
% % % Relative Position Control
% a.kp = 300;  % Proportional [kg*s^-2]
% a.kd = 50;  % Derivative [kg/s]    
% a.ki = 20;  % Integral [kg*s^-3]

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


if plot_gain
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
end

% -------------------------- Parameters -------------------------------- %

a.w = 1; % Range of inputs

% Piecewise radii
c = 0.5;
r_c = a.w/2*c; % Full isolation control radius

b = 0.5;
r_b = a.w/2*b; % Zero relative position control radius

% Polynomial order
n = 4;

% --------------------------- Mixing Functions ------------------------- %

% Center gain scale (inertial isolation)
a_c = -1 / abs(r_c - a.w/2)^n;
k_c = 1;
a.C = @(x, d) 1 .* (abs(x - d) <= r_c) ...
       + (a_c*abs(x-(d+sign(x-d)*r_c)).^n + k_c) .* (abs(x - d) > r_c & ...
       abs(x - d) <= a.w/2);

% Boundary gain scale (relative position control)
a_b = 1 / abs(r_b - a.w/2)^n;
k_b = 0;
a.B = @(x, d) (a_b*abs(x-(d+sign(x-d)*r_b)).^n + k_b) .* (abs(x - d) > r_b & ...
       abs(x - d) <= a.w/2) ...
       + 1 .* (abs(x - d) > a.w/2);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%% Performance Parameters %%%%%%%%%%%%%%%%%%%

% Maximum acceleration metric
p_ddot_max = 0.005*beta^2;

% Settling time
t_s = 2;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
