close all; clear; clc;

% Testing piecwise polynomial functions for mixed inertial and relative
% position control

pr_d = 0.5; % desired relative position (m)
r_g = 0.4; % radius from pr_d for full inertial control
h_k = 0.1; % proportion of relative position control gain inside abs(pr-(pr_d+r_k))
r_k = 0.4; % radius from pr_d for h_k proportion of relative position control

pr = linspace(0,1,1000); % relative position test values
inert_gain = zeros(size(pr)); % Proportion of inertial control gain
rel_gain = zeros(size(pr)); % Proportion of relative position control gain

n = 10; % Polynomial order
g = @(x) (-((h_k-1)/(1-(pr_d+r_k))^n) * (sign(x-pr_d)*(x-(pr_d+sign(x-pr_d)*r_k)))^n + h_k) * (abs(x-pr_d)> r_k)...
    + h_k * (abs(x-pr_d) <= r_k);

n = 10; % Polynomial order
k = @(x) ((-1/(1-(pr_d+r_g))^n) * (sign(x-pr_d)*(x-(pr_d+sign(x-pr_d)*r_g)))^n + 1) * (abs(x-pr_d)> r_g)...
    + 1 * (abs(x-pr_d) <= r_g);

for i = 1:length(pr)
    inert_gain(i) = g(pr(i));
    rel_gain(i) = k(pr(i));
end

figure;
plot(pr,inert_gain)
hold on
plot(pr, rel_gain)
title('Inertial and Relative Positional Control Mixing')
xlabel('Relative Position (m)')
ylabel('Gain Proportion')
xlim([-0.1 1.1])
ylim([-0.1 1.1])
legend('Inertial','Relative')
grid on
