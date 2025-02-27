close all;

% Converting erronous ouput to NaN
max_isolation(max_isolation <= 0) = NaN;

% Map indices to actual axis values
x1 = acceleration_gains;
y1 = rel_prop_gains;
z1 = rel_deriv_gains;
a1 = rel_int_gains;
b1 = r_g;
c1 = B_scale;
d1 = n_vals;

[M,I] = min(max_isolation(:));
[i, j, k, l, m, n, o] = ind2sub(size(max_isolation), I);
fprintf(['Minimum isolation of %.4f:\n' ...
         'Acceleration Gain: %.2f\n' ...
         'Proportional Gain: %.2f\n' ...
         'Dervative Gain: %.2f\n' ...
         'Integral Gain: %.2f\n'...
         'Radius of center control: %.3f\n' ...
         'Boundary Gan Scale: %.3f\n' ...
         'Polynomial Order: %d\n'], M, x1(i), y1(j), z1(k), a1(l), ...
                                      b1(m), c1(n), d1(o))

% [X, Y, Z] = ndgrid(x1, y1, z1);
% max_isolation_3D = max_isolation(:,:,:,l,m,n,o);

[X, Y, Z] = ndgrid(x1, z1, a1);
max_isolation_3D = reshape(max_isolation(:,j,:,:,m,n,o), ...
                                            [length(acceleration_gains), ...
                                            length(rel_deriv_gains), ...
                                            length(rel_int_gains)]);

% Plotting
figure;
hold on;
% scatter3(X(max_isolation_3D < 0.5), ...
%          Y(max_isolation_3D < 0.5), ...
%          Z(max_isolation_3D < 0.5), 50, ...
%          max_isolation_3D(max_isolation_3D < 0.5), 'filled');
validIdx = ~isnan(max_isolation_3D) & ~isinf(max_isolation_3D); % Keep only finite values
scatter3(X(validIdx), Y(validIdx), Z(validIdx), ...
                    50, max_isolation_3D(validIdx), 'filled');
colorbar;
clim([0, 1])
% colormap(spring);
colormap(hot);
xlabel('Acceleration Gain');
% ylabel('Proportional Gain');
% zlabel('Derivative Gain');
ylabel('Derivative Gain');
zlabel('Integral Gain');
grid on;
axis tight;
% axis equal;
view(3);  % 3D view
hold off;

