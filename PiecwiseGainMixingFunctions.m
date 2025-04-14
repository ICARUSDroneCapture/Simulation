close all; clear; clc;

%% Determining mixing functions

syms a k d w r

for i = 1:4

    n = i;

    % Solving for piecewise function for intertial control
    % e1 = (0 == a*abs((d-w/2) - (d-r))^n + k);
    % e2 = (1 == a*abs((d-r) - (d-r))^n + k);
    % [a_c, k_c] = solve(e1, e2, a, k, real=true);

    % Solving for piecewise function for relative position control
    e1 = (1 == a*abs((d-w/2) - (d-r))^n + k);
    e2 = (0 == a*abs((d-r) - (d-r))^n + k);
    [a_b, k_b] = solve(e1, e2, a, k, real=true);

end


%% Gain Mixing

% -------------------------- Parameters -------------------------------- %
w = 90; % Range of inputs
d = w/2; % Center of input region

% Piecewise radii
c = 0.8;
r_c = 45*c; % Full isolation control radius
b = 0.8;
r_b = 45*b; % Zero relative position control radius

% Polynomial order
n = 4;


% --------------------------- Mixing Functions ------------------------- %

% % Center gain scale (inertial isolation)
% a_c = -1 / abs(r_c - w/2)^n;
% k_c = 1;
% C = @(x) 1 .* (abs(x - d) <= r_c) ...
%        + (a_c*abs(x-(d+sign(x-d)*r_c)).^n + k_c) .* (abs(x - d) > r_c);
% 
% % Boundary gain scale (relative position control)
% a_b = 1 / abs(r_b - w/2)^n;
% k_b = 0;
% B = @(x) (a_b*abs(x-(d+sign(x-d)*r_b)).^n + k_b) .* (abs(x - d) > r_b);

% Center gain scale (inertial isolation)
a_c = -1 / abs(r_c - w/2)^n;
k_c = 1;
C = @(x) 1 .* (abs(x - d) <= r_c) ...
       + (a_c*abs(x-(d+sign(x-d)*r_c)).^n + k_c) .* (abs(x - d) > r_c & ...
       abs(x - d) <= w/2);

% Boundary gain scale (relative position control)
a_b = 1 / abs(r_b - w/2)^n;
k_b = 0;
B = @(x) (a_b*abs(x-(d+sign(x-d)*r_b)).^n + k_b) .* (abs(x - d) > r_b & ...
       abs(x - d) <= w/2) ...
       + 1 .* (abs(x - d) > w/2);


% ---------------------- Gain Mixing Visualization --------------------- %

theta = linspace(d - 1.5*w/2, d + 1.5*w/2,1000); % relative position test values
cent_gain = C(theta); % Proportion of center control gain
boundary_gain = B(theta); % Proportion of boundary control gain

% figure;
plot(theta, cent_gain)
hold on
plot(theta, boundary_gain)
xline(d,'--','Label','$x_{d}$','Interpreter','latex','FontSize',15,...
    'LabelOrientation','horizontal','LabelVerticalAlignment','middle')
title('Center and Boundary Control Mixing')
xlabel('Relative Position')
ylabel('Gain Proportion')
xlim([d - 1.1*w/2 d + 1.1*w/2])
ylim([-0.1 1.1])
% legend('Center','Boundary')
grid on
