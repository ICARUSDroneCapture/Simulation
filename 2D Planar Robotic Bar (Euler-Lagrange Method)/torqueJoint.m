function torque = torqueJoint(time,y,ref_q,gains)
%SUMMARY

% definding the symbolic variables of the system
syms q1(t)

% reference angles
ref_q1 = ref_q(1);

% gains
Kp1 = gains(1);
Ki1 = gains(2);
Kd1 = gains(3);


% initializing the torque vector
torque = zeros(1,length(time));

% writting diff(ref_q1-q1(t),'t') symbolically
diff_q1Dif = diff(ref_q1-q1(t),'t');

% writting int(ref_q1-q1(t)) as a vector
int_q1Dif = cumtrapz(time,ref_q1-y(1,:));

for i = 1:length(time)
    % evaluate diff(ref_q1-q1(t),'t') at time i
    diff_q1Dif_eval = double(subs(diff_q1Dif,diff(q1(t),t),y(2,i)));

    torque(1,i) =  Kp1*(ref_q1-y(1,i))+Ki1*int_q1Dif(i)+Kd1*diff_q1Dif_eval;
end

% plotting the torques of each joint
figure()
plot(time,torque(1,:))
xlabel('Time [s]')
ylabel('Input Torque [N.m]')
title('Input Torque of Joint 1')

end

