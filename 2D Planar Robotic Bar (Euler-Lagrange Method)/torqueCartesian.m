function torque = torqueCartesian(time, y, r_0, r_0_ref,p_ddot, angle2, gains)
%SUMMARY

% symbolic variables
syms q1(t) theta2(t)

% parameters of the system
constants;

% gains
Ka = gains(1);
Kv = gains(2);
Kp = gains(3);
Ki = gains(4);
Kd = gains(5);

% initializing a torque vector
torque = zeros(1,length(time));

% forward kinematics
r_0_eval1 = platform.l1*cos(y(1,:)+double(subs(angle2,t,time)));
r_0_eval2 = zeros(size(time));
r_0_eval3 = -platform.l1*sin(y(1,:)+double(subs(angle2,t,time)));
r_0_eval = [r_0_eval1;r_0_eval2;r_0_eval3];

% the integral of the (negative) inertial acceleration
np_dot = int(-p_ddot);

% writting int(r_0_ref-r_0) as a vector
int_positionDif = cumtrapz(time,r_0_ref-r_0_eval,2);

% writting diff(r_0_ref-r_0) symbolically
diff_positionDif = diff(r_0_ref-r_0);

% the second derivative of q1 and q2
Ddq1 = gradient(y(2,:),time);

% the first and second derivative of angle2
Dangle2 = diff(angle2,'t');
Ddangle2 = diff(angle2,'t',2);

for i=1:length(time)
    % the jocabian
    J_eval = [-platform.l1*sin(y(1,i)+double(subs(angle2,t,time(i))));0;-platform.l1*cos(y(1,i)+double(subs(angle2,t,time(i))))];
    
    % evaluating the inertial acceleration vector at time i
    p_ddot_eval = subs(p_ddot,diff(q1(t),t,t),Ddq1(i));
    p_ddot_eval = subs(p_ddot_eval,diff(theta2(t),t,t),double(subs(Ddangle2,t,time(i))));
    p_ddot_eval = subs(p_ddot_eval,diff(q1(t),t),y(2,i));
    p_ddot_eval = subs(p_ddot_eval,diff(theta2(t),t),double(subs(Dangle2,t,time(i))));
    p_ddot_eval = subs(p_ddot_eval,q1(t),y(1,i));
    p_ddot_eval = subs(p_ddot_eval,theta2(t),double(subs(angle2,t,time(i))));
    p_ddot_eval = double(subs(p_ddot_eval,t,time(i)));

    % evaluating the integral of the  (negative) inertial acceleration
    % vector at time i
    np_dot_eval = subs(np_dot,diff(q1(t),t),y(2,i));
    np_dot_eval = subs(np_dot_eval,diff(theta2(t),t),double(subs(Dangle2,t,time(i))));
    np_dot_eval = subs(np_dot_eval,q1(t),y(1,i));
    np_dot_eval = subs(np_dot_eval,theta2(t),double(subs(angle2,t,time(i))));
    np_dot_eval = double(subs(np_dot_eval,t,time(i)));

    % evaluating the diff_positionDif at time i
    diff_positionDif_eval = subs(diff_positionDif,diff(q1(t),t),y(2,i));
    diff_positionDif_eval = subs(diff_positionDif_eval,diff(theta2(t),t),double(subs(Dangle2,t,time(i))));
    diff_positionDif_eval = subs(diff_positionDif_eval,q1(t),y(1,i));
    diff_positionDif_eval = subs(diff_positionDif_eval,theta2(t),double(subs(angle2,t,time(i))));
    diff_positionDif_eval = double(subs(diff_positionDif_eval,t,time(i)));
    
    % evaluating the control force on the end effector
    F_control_eval = Kp*(r_0_ref-r_0_eval(:,i))+Ki*int_positionDif(:,i)+Kd*diff_positionDif_eval-Ka*p_ddot_eval+Kv*np_dot_eval;
    torque(:,i) = (J_eval')*F_control_eval;
end

% plotting the torques of each joint
figure()
plot(time,torque(1,:))
xlabel('Time [s]')
ylabel('Input Torque [N.m]')
title('Input Torque of Joint 1')

end

