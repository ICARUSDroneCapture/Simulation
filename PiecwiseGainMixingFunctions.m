close all; clear; clc;

% syms a k h w r
% % syms n
% % assume(n > 0)
% % assumeAlso(n, "integer")
% % assumptions(n)
% 
% for i = 1:4
% 
%     n = i;
% 
%     % Solving for piecewise function for intertial control
%     % e1 = (0 == a*(h - 0)^n + k);
%     % e2 = (1 == a*(h - (w/2-r))^n + k);
%     % e3 = (0 == a*(w - h)^n + k);
%     % [a_c, k_c, h_c] = solve(e1, e2, e3, a, k, h, real=true);
% 
%     % Solving for piecewise function for relative position control
%     % e1 = (1 == a*(h - 0)^n + k);
%     % e2 = (0 == a*(h - (w/2-r))^n + k);
%     % e3 = (1 == a*(w - h)^n + k);
%     % [a_b, k_b, h_b] = solve(e1, e2, e3, a, k, h, real=true);
% 
% end

%% Gain Mixing

% Parameters
w = 90; % Range of inputs
h = w/2; % Center of input region

% Piecewise radii
c = 0.8;
r_c = 45*c; % Full isolation control radius
d = 0.8;
r_b = 45*d; % Zero relative position control radius

% Polynomial order
n = 8;

% Center gain scale (inertial isolation)
a_c = 2^n / ((2*r_c)^n - w^n);
k_c = -w^n / ((2*r_c)^n - w^n);
C = @(x) 1 .* (abs(x - h) <= r_c) ...
                + (a_c*abs(x-h).^n + k_c) .* (abs(x - h) > r_c);

% Boundary gain scale (relative position control)
a_b = -2^n / ((2*r_b)^n - w^n);
k_b = (2*r_b)^n / ((2*r_b)^n - w^n);
B = @(x) (a_b*abs(x-h).^n + k_b) .* (abs(x - h) > r_b);

%%%               Gain Mixing Visualization               %%%

theta = linspace(0,w,1000); % relative position test values
cent_gain = C(theta); % Proportion of center control gain
boundary_gain = B(theta); % Proportion of boundary control gain

% figure;
plot(theta, cent_gain)
hold on
plot(theta, boundary_gain)
xline(h,'--','Label','$\theta_{rd}$','Interpreter','latex','FontSize',15,...
    'LabelOrientation','horizontal','LabelVerticalAlignment','middle')
title('Center and Boundary Control Mixing')
xlabel('Relative Position (deg)')
ylabel('Gain Proportion')
xlim([-0.1*w 1.1*w])
ylim([-0.1 1.1])
legend('Center','Boundary')
grid on
