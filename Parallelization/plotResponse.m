function plotResponse(t, s, betas, isolation, a, max_isolation)
%PLOTRESPONSE Summary of this function goes here
%   Detailed explanation goes here

screenSize = get(0, 'ScreenSize'); 
width = screenSize(3) / 4;  % Half the screen width
height = screenSize(4) / 4; % Half the screen height

% Position
figure('Position', [100, 100, width, height]);
plot(t,s(:,1))
hold on
plot(t,a.d(t))
title('Inertial Position vs Time')
xlabel('Time (s)')
ylabel('Position (m)')
legend('Platform', 'Deck','Location','northeast')

figure('Position', [100, height + 200, width, height]);
plot(t,s(:,1)-a.d(t))
title('Relative Position vs Time')
xlabel('Time (s)')
ylabel('Position (m)')

figure('Position', [width + 150, 100, width, height]);
plot(t,s(:,2))
hold on
plot(t,a.d_dot(t))
title('Inertial Velocity vs Time')
xlabel('Time (s)')
ylabel('Velocity (m/s)')
legend('Platform', 'Deck','Location','southeast')

figure('Position', [width + 150, height + 200, width, height]);
plot(t,s(:,2)-a.d_dot(t))
title('Relative Velocity vs Time')
xlabel('Time (s)')
ylabel('Velocity (m/s)')

figure('Position', [2*width + 150, 100, width, height]);
% Feeding states back through EOM to calculating inertial acceleration of
% the platfor
p_ddot = zeros(size(t));
for i = 1:length(t)
    s_dot = rigidArmControl(t(i),s(i,:),a);
    p_ddot(i) = s_dot(2);
end
plot(t,p_ddot)
hold on
plot(t,a.d_ddot(t))
title('Inertial Acceleration vs Time')
xlabel('Time (s)')
ylabel('Acceleration (m/s^2)')
legend('Platform', 'Deck','Location','southeast')

figure('Position', [2*width + 150, height + 200, width, height]);
plot(betas,isolation)
hold on
yline(0.1, '--')
title('Isolation vs Wave Frequency')
xlabel('Angular Frequency (rad/s)')
ylabel('Isolation Ratio x_I/x_D')

fprintf('Maximum Isolation above minimum wave frequency is %.3f\n',max_isolation)

end

