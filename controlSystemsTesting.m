
set(0, 'DefaultLineLineWidth', 1);

% Simulation time
tspan = [0 30]; % [s]

% Initial States
p0 =  a.d(tspan(1))+a.pr_d; % Platform position [m]
p_dot0 = a.d_dot(tspan(1)); % Platform velocity [m/s]
pr_err_accum0 = 0;          % Integral of error in relative position [m*s]
s0 = [p0; p_dot0; pr_err_accum0];

% Running Simulation
lastwarn('')
op = odeset('RelTol',1e-5,'AbsTol',1e-5); % Tolerance options
% [t, s] = ode45(@(t,s)rigidArmControl(t,s,a),tspan,s0,op);
try
[t, s] = ode45(@(t,s)rigidArmControl(t,s,a),tspan,s0,op);
catch ME
    if strcmp(ME.identifier, warnId)
        % disp(beta)
        % disp(['Solver stopped due to tolerance failure at t = ', num2str(ME.stack(1).line)]);
        t = 0;

    else
        % Something else has gone wrong: just re-throw the error
        throw(ME);
    end
end

% Plotting Position vs Time and Acceleration vs Time

% % Position
% figure;
% plot(t,s(:,1))
% hold on
% plot(t,a.d(t))
% title('Inertial Position vs Time')
% xlabel('Time (s)')
% ylabel('Position (m)')
% legend('Platform', 'Deck','Location','southeast')
% 
% % Velocity
% figure;
% plot(t,s(:,2))
% hold on
% plot(t,a.d_dot(t))
% title('Inertial Velocity vs Time')
% xlabel('Time (s)')
% ylabel('Velocity (m/s)')
% legend('Platform', 'Deck','Location','southeast')
% 
% % Acceleration
% figure;
% % Feeding states back through EOM to calculating inertial acceleration of
% % the platfor
% p_ddot = zeros(size(t));
% f_comp = m0*a.g; % Gravity compensation force mass [N]
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
% 
% 
% % Plotting relative position
% figure;
% plot(t,s(:,1)-a.d(t))
% hold on
% % Plotting inertial control region
% x = [tspan, flip(tspan)];
% yf = [a.pr_d-a.r_g, a.pr_d-a.r_g, a.pr_d+a.r_g, a.pr_d+a.r_g];
% fill(x,yf,'y','FaceAlpha',0.2,'EdgeColor','none')
% % Plotting Relative position control region
% x = [tspan, flip(tspan)];
% yta = [a.pr_d+a.r_k, a.pr_d+a.r_k, 1, 1];
% ytb = [0, 0, a.pr_d-a.r_k, a.pr_d-a.r_k];
% fill(x,yta,'b','FaceAlpha',0.2,'EdgeColor','none')
% fill(x,ytb,'b','FaceAlpha',0.2,'EdgeColor','none')
% yline(a.pr_d,'--','Label','$p_{rd}$','Interpreter','latex','FontSize',15)
% title('Relative Position vs Time')
% xlabel('Time (s)')
% ylabel('Position (m)')
% legend('','Full Inertial','Relative Position','')
