
out = sim("Control1DOF.slx",'StopTime', '30');

t = out.tout;
p = out.p.Data;

% Position
figure;
plot(t,p)
hold on
plot(t,a.d(t))
title('Inertial Position vs Time')
xlabel('Time (s)')
ylabel('Position (m)')
legend('Platform', 'Deck','Location','southeast')