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
beta_wave = 2*pi/Tmax;
[~, wave_idx] = min(abs(betas-beta_wave));

% acceleration_gains = 500:500:10000;
% rel_prop_gains = 50:50:1000;
% rel_deriv_gains = 500:500:3500;
acceleration_gains = 29000;
rel_prop_gains = 200;
rel_deriv_gains = 3500;

good_control = -1*ones(length(acceleration_gains), ...
                       length(rel_prop_gains), ...
                       length(rel_deriv_gains));
max_isolation = zeros(size(good_control));

% originalState = warning;
% warnId = 'MATLAB:ode45:IntegrationTolNotMet';
% warning('error', warnId);

tic

for d_idx = 1:length(rel_deriv_gains)

a.kd_c = rel_deriv_gains(d_idx);

for p_idx = 1:length(rel_prop_gains)

a.kp_c = rel_prop_gains(p_idx);

for a_idx = 1:length(acceleration_gains)

a.ka_c = acceleration_gains(a_idx);

for idx = 1:n_max
    beta = betas(idx);
    a.d = @(t) alpha*cos(beta*t) + hdeck;      % [m]
    a.d_dot = @(t) -beta*alpha*sin(beta*t);    % [m/s]
    a.d_ddot = @(t) -beta^2*alpha*cos(beta*t); % [m*s^-2]
    
    % Simulating model
    controlSystemsTesting;
    
    if idx == wave_idx
    % Position
    figure('Position', [100, 100, width, height]);
    plot(t,s(:,1))
    hold on
    plot(t,a.d(t))
    title('Inertial Position vs Time')
    xlabel('Time (s)')
    ylabel('Position (m)')
    legend('Platform', 'Deck','Location','southeast')

    figure('Position', [width + 150, 100, width, height]);
    plot(t,s(:,2))
    hold on
    plot(t,a.d_dot(t))
    title('Inertial Velocity vs Time')
    xlabel('Time (s)')
    ylabel('Velocity (m/s)')
    legend('Platform', 'Deck','Location','southeast')

    figure('Position', [width + 150, height + 150, width, height]);
    plot(t,s(:,2)-a.d_dot(t))
    title('Relative Velocity vs Time')
    xlabel('Time (s)')
    ylabel('Velocity (m/s)')

    figure('Position', [2*width + 150, 100, width, height]);
    % Feeding states back through EOM to calculating inertial acceleration of
    % the platfor
    p_ddot = zeros(size(t));
    f_comp = m0*a.g; % Gravity compensation force mass [N]
    for i = 1:length(t)
        s_dot = rigidArmControl(t(i),s(i,:),a);
        p_ddot(i) = s_dot(2);
    end
    plot(t,p_ddot)
    hold on
    plot(t,a.d_ddot(t))
    yline(p_ddot_max,'--')
    yline(-p_ddot_max,'--')
    % ylim([min(a.d_ddot(t))*1.25 max(a.d_ddot(t))*1.25])
    title('Inertial Acceleration vs Time')
    xlabel('Time (s)')
    ylabel('Acceleration (m/s^2)')
    legend('Platform', 'Deck','Location','southeast')
    end
    
    % Isolation
    t_test0 = 20;
    x_I = s(t>t_test0,1);
    x_D = a.d(t(t>t_test0));
    if sum(t>t_test0) > 0
    isolation(idx) = calculateAverageIsolation(x_I, x_D);
    else
        isolation(idx) = -1;
    end
end

max_isolation(a_idx, p_idx, d_idx) = max(isolation(betas >= beta_wave))
% if max_isolation(a_idx, p_idx, d_idx) < 0.4
%     good_control(a_idx, p_idx, d_idx) = 1;
% else
%     good_control(a_idx, p_idx, d_idx) = 0;
% end

figure;
plot(betas,isolation)
hold on
yline(0.1, '--', 'DisplayName','')
title('Isolation vs Wave Frequency')
xlabel('Angular Frequency (rad/s)')
ylabel('Isolation Ratio x_I/x_D')
legend();

end
toc
% disp('---- One proportional gain done -----')
end
fprintf('One derivative gain done: number %d \n',d_idx)
end

% warning(originalState);

% save('very_wide_good_control_info_V2', "max_isolation", "acceleration_gains", ...
%     "rel_prop_gains", "rel_deriv_gains")






