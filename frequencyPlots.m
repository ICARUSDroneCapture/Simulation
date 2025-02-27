close all; clear; clc;

simulationParameters;

screenSize = get(0, 'ScreenSize'); 
width = screenSize(3) / 4;  % Half the screen width
height = screenSize(4) / 4; % Half the screen height

beta_min = 0.5;
beta_max = 2;
n_max = 15;
betas = linspace(beta_min, beta_max, n_max);
isolation = zeros(1,n_max);
max_wave = 2*pi/Tmax;
beta_wave = max_wave;
[~, wave_idx] = min(abs(betas-beta_wave));

% Computation time: 12.8444 hrs
acceleration_gains = 4000:500:15000;
rel_prop_gains = 50:50:1500;
rel_deriv_gains = 2500:500:6000;
% B_scale = 0:0.2:1;
% r_g = 0.2:0.1:0.5; 
% n_vals = 1:4;
% acceleration_gains = 12500;
% rel_prop_gains = 800;
% rel_deriv_gains = 3000;
B_scale = 0.2;
r_g = 0.4; 
n_vals = 4;

good_control = -1*ones(length(acceleration_gains), ...
                       length(rel_prop_gains), ...
                       length(rel_deriv_gains), ...
                       length(r_g), ...
                       length(B_scale), ...
                       length(n_vals));
max_isolation = zeros(size(good_control));

originalState = warning;
warnId = 'MATLAB:ode45:IntegrationTolNotMet';
warning('error', warnId);

tic

for n_idx = 1:length(n_vals)

n = n_vals(n_idx);
a.B = @(x) ...
    (-((a.h_k-1)/(1-(a.pr_d+a.r_k))^n) * ...
        (sign(x-a.pr_d).*(x-(a.pr_d+sign(x-a.pr_d)*a.r_k))).^n + a.h_k) .*...
            (abs(x-a.pr_d)> a.r_k) ...
    + a.h_k * (abs(x-a.pr_d) <= a.r_k);

for B_idx = 1:length(B_scale)

% Relative Position Control at boundaries
a.kp_b = B_scale(B_idx)*3000;  % Proportional [kg*s^-2]
a.kd_b = B_scale(B_idx)*500;  % Derivative [kg/s]    
a.ki_b = B_scale(B_idx)*200;  % Integral [kg*s^-3] 

for rg_idx = 1:length(r_g)

a.r_g = r_g(rg_idx);
a.C = @(x) ...
    ((-1/(1-(a.pr_d+a.r_g))^n) * ...
        (sign(x-a.pr_d).*(x-(a.pr_d+sign(x-a.pr_d)*a.r_g))).^n + 1) .* ...
            (abs(x-a.pr_d)> a.r_g) ...
    + 1 * (abs(x-a.pr_d) <= a.r_g);

for d_idx = 1:length(rel_deriv_gains)

a.kd_c = rel_deriv_gains(d_idx);

for p_idx = 1:length(rel_prop_gains)

a.kp_c = rel_prop_gains(p_idx);

for a_idx = 1:length(acceleration_gains)

a.ka = acceleration_gains(a_idx);

for idx = 1:n_max
    beta = betas(idx);
    a.d = @(t) alpha*cos(beta*t) + hdeck;      % [m]
    a.d_dot = @(t) -beta*alpha*sin(beta*t);    % [m/s]
    a.d_ddot = @(t) -beta^2*alpha*cos(beta*t); % [m*s^-2]
    
    % Simulating model
    controlSystemsTesting;
    
    % if idx == wave_idx
    % % Position
    % figure('Position', [100, 100, width, height]);
    % plot(t,s(:,1))
    % hold on
    % plot(t,a.d(t))
    % title('Inertial Position vs Time')
    % xlabel('Time (s)')
    % ylabel('Position (m)')
    % legend('Platform', 'Deck','Location','northeast')
    % 
    % figure('Position', [100, height + 150, width, height]);
    % plot(t,s(:,1)-a.d(t))
    % title('Relative Position vs Time')
    % xlabel('Time (s)')
    % ylabel('Position (m)')
    % 
    % figure('Position', [width + 150, 100, width, height]);
    % plot(t,s(:,2))
    % hold on
    % plot(t,a.d_dot(t))
    % title('Inertial Velocity vs Time')
    % xlabel('Time (s)')
    % ylabel('Velocity (m/s)')
    % legend('Platform', 'Deck','Location','southeast')
    % 
    % figure('Position', [width + 150, height + 150, width, height]);
    % plot(t,s(:,2)-a.d_dot(t))
    % title('Relative Velocity vs Time')
    % xlabel('Time (s)')
    % ylabel('Velocity (m/s)')
    % 
    % figure('Position', [2*width + 150, 100, width, height]);
    % % Feeding states back through EOM to calculating inertial acceleration of
    % % the platfor
    % p_ddot = zeros(size(t));
    % for i = 1:length(t)
    %     s_dot = rigidArmControl(t(i),s(i,:),a);
    %     p_ddot(i) = s_dot(2);
    % end
    % plot(t,p_ddot)
    % hold on
    % plot(t,a.d_ddot(t))
    % yline(p_ddot_max,'--')
    % yline(-p_ddot_max,'--')
    % % ylim([min(a.d_ddot(t))*1.25 max(a.d_ddot(t))*1.25])
    % title('Inertial Acceleration vs Time')
    % xlabel('Time (s)')
    % ylabel('Acceleration (m/s^2)')
    % legend('Platform', 'Deck','Location','southeast')
    % end
    
    % Isolation
    t_test0 = 20;
    if t(end) == tspan(end)
        x_I = s(t>t_test0,1);
        x_D = a.d(t(t>t_test0));
        isolation(idx) = calculateAverageIsolation(x_I, x_D);
        if isolation(idx) < 0
            isolation(idx) = NaN;
        end
    else
        isolation(idx) = NaN;
    end

end

if sum(isnan(isolation)) == 0
    max_isolation(a_idx, p_idx, d_idx, rg_idx, B_idx, n_idx) = max(isolation(betas >= max_wave));
else
    max_isolation(a_idx, p_idx, d_idx, rg_idx, B_idx, n_idx) = NaN;
end
% if max_isolation(a_idx, p_idx, d_idx, rg_idx, B_idx) < 0.4
%     good_control(a_idx, p_idx, d_idx, rg_idx, B_idx) = 1;
% else
%     good_control(a_idx, p_idx, d_idx, rg_idx, B_idx) = 0;
% end

% figure('Position', [2*width + 150, height + 150, width, height]);
% plot(betas,isolation)
% hold on
% yline(0.1, '--', 'DisplayName','')
% title('Isolation vs Wave Frequency')
% xlabel('Angular Frequency (rad/s)')
% ylabel('Isolation Ratio x_I/x_D')
% legend();

end
toc
% disp('---- One proportional gain done -----')
end
fprintf(['One derivative gain done: number %d\n' ...
         'r_g number: %d\n' ...
         'B_scale number: %d\n' ...
         'Polynomial order number: %d\n'], d_idx, rg_idx, B_idx, n_idx)
end
end
end
end

warning(originalState);

% save('test_good_control_info', "max_isolation", "acceleration_gains", ...
%     "rel_prop_gains", "rel_deriv_gains", "r_g", "B_scale")








