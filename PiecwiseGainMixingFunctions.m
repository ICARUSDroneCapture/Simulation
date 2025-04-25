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
w = 1; % Range of inputs
d = w/2; % Center of input region

% Piecewise radii
c = 0.5;
r_c = 0.5*c; % Full isolation control radius
b = 0.5;
r_b = 0.5*b; % Zero relative position control radius

% Polynomial order
n = 1;


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

m_lim = 1.14;
theta = linspace(d - m_lim*w/2, d + m_lim*w/2,1000); % relative position test values
cent_gain = C(theta); % Proportion of center control gain
boundary_gain = B(theta); % Proportion of boundary control gain

figure;
fill([d-r_c d-r_c d+r_c d+r_c], [1.1 -0.1 -0.1 1.1], "r", 'FaceColor', "#b6d7a8", ...
    'FaceAlpha', 0.4, 'EdgeColor','none');
hold on
fill([d-w/2 d-w/2 d-r_c d-r_c], [1.1 -0.1 -0.1 1.1], "r", 'FaceColor', "#ffe59a", ...
    'FaceAlpha', 0.4, 'EdgeColor','none');
fill([d+w/2 d+w/2 d+r_c d+r_c], [1.1 -0.1 -0.1 1.1], "r", 'FaceColor', "#ffe59a", ...
    'FaceAlpha', 0.4, 'EdgeColor','none');
fill([d-m_lim*w/2 d-m_lim*w/2 d-w/2 d-w/2], [1.1 -0.1 -0.1 1.1], "r", 'FaceColor', "#ea9998", ...
    'FaceAlpha', 0.4, 'EdgeColor','none');
fill([d+m_lim*w/2 d+m_lim*w/2 d+w/2 d+w/2], [1.1 -0.1 -0.1 1.1], "r", 'FaceColor', "#ea9998", ...
    'FaceAlpha', 0.4, 'EdgeColor','none');
plot(theta, cent_gain, 'LineWidth', 1.5)
text(d+0.2*w, 0, 'Boundary Control Scale', 'FontSize', 11, ...
    'HorizontalAlignment', 'right', 'VerticalAlignment', 'top');
plot(theta, boundary_gain, 'LineWidth', 1.5)
text(d-0.18*w, 1, 'Inertial Control Scale', 'FontSize', 11, ...
    'HorizontalAlignment', 'left', 'VerticalAlignment', 'bottom');
xline(d,'--','Label',{'Desired'; 'Position'},...
    'LabelHorizontalAlignment','right',...
    'LabelVerticalAlignment','middle',...
    'LabelOrientation','horizontal',...
    'FontSize', 11, 'LineWidth', 1)
title('Inertial Stability and Boundary Control Mixing')
xlabel('Relative Position')
ylabel('Gain Proportion')
xlim([d - m_lim*w/2 d + m_lim*w/2])
ylim([-0.1 1.1])
% legend('Center Control Scale','Boundary Control Scale', 'Desired Position', ...
%        'Location','eastoutside')
grid on
