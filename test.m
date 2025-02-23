close all; clear; clc;

% for i = 1:length(rel_deriv_gains)
%     disp(i)
%     max_isolation(:,:,i)
% end

% % Get indices of 1's and 0's
% [a1_idx, p1_idx, d1_idx] = ind2sub(size(good_control), find(good_control == 1)); % 1's
% [a0_idx, p0_idx, d0_idx] = ind2sub(size(good_control), find(good_control == 0)); % 0's
% 
% % Map indices to actual axis values
% x1 = acceleration_gains(a1_idx);
% y1 = rel_prop_gains(p1_idx);
% z1 = rel_deriv_gains(d1_idx);
% 
% x0 = acceleration_gains(a0_idx);
% y0 = rel_prop_gains(p0_idx);
% z0 = rel_deriv_gains(d0_idx);
% 
% % Plotting
% figure;
% hold on;
% scatter3(x1, y1, z1, 50, 'g', 'filled');  % Green for 1's
% scatter3(x0, y0, z0, 50, 'r', 'filled');  % Red for 0's
% xlabel('Acceleration Gain');
% ylabel('Proportional Gain');
% zlabel('Derivative Gain');
% grid on;
% axis tight;
% axis equal;
% view(3);  % 3D view
% hold off;

data_files = ["below_good_control_info.mat",...
              "above_good_control_info.mat",... 
              "good_control_info.mat",....
              "side_good_control_info.mat",...
              "very_below_good_control_info.mat",...
              "very_wide_good_control_info.mat"];

for d = 1:length(data_files)
    load(data_files(d))
    % Map indices to actual axis values
    x1 = acceleration_gains;
    y1 = rel_prop_gains;
    z1 = rel_deriv_gains;

    [M,I] = min(max_isolation(:));
    [i, j, k] = ind2sub(size(max_isolation), I);
    fprintf(['Minimum isolation of %.4f:\n' ...
             'Acceleration Gain: %d\n' ...
             'Proportional Gain: %d\n' ...
             'Dervative Gain: %d\n'], M, x1(i), y1(j), z1(k))

    [X, Y, Z] = ndgrid(x1, y1, z1);

    % Plotting
    figure(1);
    hold on;
    scatter3(X(max_isolation < 0.5), Y(max_isolation < 0.5), Z(max_isolation < 0.5), 50, max_isolation(max_isolation < 0.5), 'filled');
end

% % Converting erronous ouput to NaN
% max_isolation(max_isolation == -1) = NaN;
% 
% % Map indices to actual axis values
% x1 = acceleration_gains;
% y1 = rel_prop_gains;
% z1 = rel_deriv_gains;
% 
% [X, Y, Z] = ndgrid(x1, y1, z1);
% 
% % Plotting
% figure;
% hold on;
% scatter3(X(max_isolation < 0.5), Y(max_isolation < 0.5), Z(max_isolation < 0.5), 50, max_isolation(max_isolation < 0.5), 'filled');
colorbar;
clim([0, 0.5])
% colormap(spring);
colormap(hot);
xlabel('Acceleration Gain');
ylabel('Proportional Gain');
zlabel('Derivative Gain');
grid on;
axis tight;
% axis equal;
view(3);  % 3D view
hold off;

