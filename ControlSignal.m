function s_dot = ControlSignal(time_i, a, prev_state)

    % Current states

    vel = prev_state(1);
    theta = prev_state(2);
    p = prev_state(3);
    theta_err_accum = prev_state(4);
    pr_err_accum = prev_state(5);

    
    % --------------------- Control Sensor Signals ------------------------

    curr_state = [vel; theta];
    state_err_accum = [p; theta_err_accum];

    specs = a.specs;

    kw = a.kw; % 
    kt = a.kt; % 

    state_0 = 0;
    state_dot_0 = 0;

    accel_i = a.real_accel(time_i);
    ang_rate_i = a.real_ang_rate(time_i);
    
    accel_dot_m_v = a.measured_accel_vert(time_i, a.o_d_n_a_c_v, a.biasStabDistAccel, a.biasTempDistAccel, a.accel_drift_vert, a.noiseDistAccel, accel_i, a.theta_err);
    gyro_dot_m = a.measured_gyro(time_i, a.o_d_n_g_c, a.biasStabDistGyro, a.biasTempDistGyro, a.gyro_drift, a.noiseDistGyro, ang_rate_i);
    
    % No signal error - test logic
    accel_dot_m_v = -accel_i;
    gyro_dot_m = -ang_rate_i;
    % ----------------------------

    state_dot_m = [accel_dot_m_v; gyro_dot_m];

    state_control = curr_state - state_0;

    state_dot_comp = kt * state_err_accum + kw * state_control;

    state_dot = state_dot_m + state_dot_0 - state_dot_comp;

    s_dot = zeros(length(prev_state),1);


    % --------------------------- Arm Control -----------------------------

    % Error in relative position (distance to center of operation region)
    % pr = p-a.d(time_i);

    pr = p-0;
    pr_err = p-a.pr_d;
    
    % Control gain proportions

    I = a.I(pr); % Proportion of inertial stability control to apply

    ka = a.ka*I; % Acceleration [kg]
    kv = a.kv*I; % Velocity     [kg/s]
    ks = a.ks*I; % Position     [kg*s^-2]

    k_h = a.K_h(pr);    % Proportion of relative position control to apply

    kp = a.kp*k_h;     % Proportional [kg*s^-2]
    kd = a.kd*k_h;     % Derivative   [kg/s]
    ki = a.ki*k_h;     % Integral     [kg*s^-3]

    % Control Law

    % Inertial stability control force
    % c_i = a.initial_scale(time_i); % Initial scale of gains
    c_i = 1;
    f_i = -(ka*accel_dot_m_v + kv*vel + ks*p)*c_i;
    % f_i = 0.1
    % Relative position control force
    f_pr = -(kp*pr_err + ki*pr_err_accum + kd*(vel-a.d_dot(time_i)));

    % Platform EOM
    p_ddot = (f_i+f_pr) / a.m;
    state_dot(1) = p_ddot;
    
    % ------------------------ Assign Derivatives -------------------------

    s_dot(1:2) = state_dot;
    s_dot(3:4) = state_control;

    % Error in relative position
    s_dot(5) = pr_err;

end