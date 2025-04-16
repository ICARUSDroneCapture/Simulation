function state_dot = rigidArmControl_3D(t, a, prev_state)
    % rigidArmControl is the EOM for the 1 DOF model of the inertially
    % stabilized platform. It uses inertial acceleration control when the
    % platform is close to the center of the operation region, and uses PID
    % control on the relative position as the platform goes closer to the
    % operational bounderies
    %
    % Inputs:   t    = current time
    %           s    = vector of states
    %                = [p; p_dot; pr_err_accum] where p is the inertial position 
    %                  of the platform, pdot is the inertial velocity of the 
    %                  platform, and pr_err_accum is the integral of the error 
    %                  in the relative position of the platform
    %           a    = structure containing environmental constants and gain
    %                  values
    % Outputs:  sdot = time derivative of input state vector
    %                = [p_dot; p_ddot; pr_err] where pdot is the inertial 
    %                  velocity of the platform, p_ddot is the inertial 
    %                  acceleration of the platform,and pr_err is the error in 
    %                  the relative position of the platform
    
    % Current states

    p_x = prev_state(1);
    p_y = prev_state(2);
    p_z = prev_state(3);

    p_dot_x = prev_state(4);
    p_dot_y = prev_state(5);
    p_dot_z = prev_state(6);

    pr_err_accum_x = prev_state(7);
    pr_err_accum_y = prev_state(8);
    pr_err_accum_z = prev_state(9);

    pm_x = prev_state(10);
    pm_y = prev_state(11);
    pm_z = prev_state(12);

    pm_dot_x = prev_state(13);
    pm_dot_y = prev_state(14);
    pm_dot_z = prev_state(15);

    pm_ddot_x = prev_state(16);
    pm_ddot_y = prev_state(17);
    pm_ddot_z = prev_state(18);

    theta = prev_state(19);
    phi = prev_state(20);
    psi = prev_state(21);

    p_theta_err_accum = prev_state(22);
    p_phi_err_accum = prev_state(23);
    p_psi_err_accum = prev_state(24);

    specs = a.specs;

    p = [p_x; p_y; p_z];
    p_dot = [p_dot_x; p_dot_y; p_dot_z];
    pr_err_accum = [pr_err_accum_x; pr_err_accum_y; pr_err_accum_z];

    pm = [pm_x; pm_y; pm_z];
    pm_dot = [pm_dot_x; pm_dot_y; pm_dot_z];
    pm_ddot = [pm_ddot_x; pm_ddot_y; pm_ddot_z];

    p_theta = [theta; phi; psi];
    p_theta_err_accum = [p_theta_err_accum; p_phi_err_accum; p_psi_err_accum];

    % ------------------------- Check Reference Frames --------------------

    % Error in relative position (distance to center of operation region)
    pr = pm-a.real_pos_zI(t);
    pr_err = pr-a.pr_d;
    
    % ----------------- Inserting measured accel manually -----------------

    % Get REAL angle
    theta_real = a.theta(t);
    phi_real = a.phi(t);
    psi_real = a.psi(t);
    
    % Get REAL inertial acceleration (this has g in it)
    a_I = [a.real_accel_xI(t); a.real_accel_yI(t); a.real_accel_zI(t)];
    
    % Rotate inertial acceleration to sensor frame base on real angle/acceleration
    a_S = Rotate_I_S(a_I, theta_real, phi_real, psi_real);

    % Add sensor error to acceleration measurements
    accel_S = a.measured_accel_3D(a, t, a_S);
    % accel_S = a_S;

    % Get REAL angular velocity
    angle_d = [a.theta_dot(t); a.phi_dot(t); a.psi_dot(t)];

    % Add sensor error to gyroscope measurements
    gyro = a.measured_gyro_3D(a, t, angle_d);
    % gyro = angle_d;

    % Rotate realistic sensor acceleration measurements back to inertial
    % frame (still has g), using our INTEGRATED angle (has integration error)
    accel_I = Rotate_S_I(accel_S, theta_real, phi_real, psi_real);
    % accel_I = Rotate_S_I(accel_S, theta_real, phi_real, psi_real);
    % accel_I = a_I;

    % Simply remove g from inertial z vector
    % g is positive in simulation parameters. When measured it would be
    % negative though, hence why we add it here
    accel_I = [accel_I(1) accel_I(2) accel_I(3)+a.g];

    % Assign our acceleration and gyroscope measurements, with gravity
    % removed, to our vector for sensor error correction
    measured_state = [accel_I, gyro'];

    % Compensate for constant error values
    corrected_state = compensateError(measured_state, specs, t);

    int_state = [-p_dot; p_theta; pr_err; p_theta_err_accum];
    
    % input: [velocity; theta; position; theta_err_accum]
    % output: [accel; theta_dot; vel; theta]
    state_control = DriftCorrection(a, int_state, corrected_state);

    pm_ddot = state_control(1:3);
    p_thetadot = state_control(4:6);

    pm_dot = state_control(7:9);
    p_theta_err = state_control(10:12);
    
    % ---------------------------------------------------------------------
    
    % Control gain proportions
    
    I = a.I(pr); % Proportion of inertial stability control to apply
    % I = 0.5;
    
    ka = a.ka*I; % Acceleration [kg]
    kv = a.kv*I; % Velocity     [kg/s]
    ks = a.ks*I; % Position     [kg*s^-2]
    
    k = a.K(pr);     % Proportion of relative position control to apply
    k_h = a.K_h(pr);
    % k_h =  0.5;
    
    kp = a.kp*k_h;     % Proportional [kg*s^-2]
    kd = a.kd*k_h;     % Derivative   [kg/s]
    ki = a.ki*k_h;     % Integral     [kg*s^-3]
    
    % Derivative of states
    state_dot = zeros(24,1);
    
    % Derivative of position
    state_dot(1:3) = p_dot;  % inertial velocity
    state_dot(10:12) = pm_dot; % measured inertial velocity
    
    % Control Law
    
    % Inertial stability control force
    c_i = a.initial_scale(t); % Initial scale of gains
    % c_i = 1;
    f_i = -(ka.*pm_ddot + kv.*pm_dot + ks.*pm).*c_i;
    % f_i = 0.1
    % Relative position control force
    f_pr = -(kp.*pr_err + ki.*pr_err_accum + kd.*(pm_dot-a.real_vel_zI(t)));
    
    % Platform EOM
    p_ddot = (f_i+f_pr) / a.m;
    
    % Derivative of velocity
    state_dot(4:6) = p_ddot; % Inertial acceleration
    state_dot(13:15) = pm_ddot; % Measured inertial acceleration
    
    % Derivative of measured inertial acceleration
    state_dot(16:18) = a.omega*(p_ddot - pm_ddot);
    
    % Error in relative position
    state_dot(7:9) = pr_err;
    
    % Derivative of angle
    state_dot(19:21) = p_thetadot; % Angular Rate
    
    % Theta error
    state_dot(22:24) = p_theta_err;

end