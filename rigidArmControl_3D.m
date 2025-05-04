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

    pi_x = prev_state(1);
    pi_y = prev_state(2);
    pi_z = prev_state(3);

    p_dot_x = prev_state(4);
    p_dot_y = prev_state(5);
    p_dot_z = prev_state(6);

    p_err_accum_x = prev_state(7);
    p_err_accum_y = prev_state(8);
    p_err_accum_z = prev_state(9);

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

    pi = [pi_x; pi_y; pi_z];
    p_dot = [p_dot_x; p_dot_y; p_dot_z];
    p_err_accum = [p_err_accum_x; p_err_accum_y; p_err_accum_z];

    pm = [pm_x; pm_y; pm_z];
    pm_dot = [pm_dot_x; pm_dot_y; pm_dot_z];
    pm_ddot = [pm_ddot_x; pm_ddot_y; pm_ddot_z];
    
    p_theta = [theta; phi; psi];
    p_theta_err_accum = [p_theta_err_accum; p_phi_err_accum; p_psi_err_accum];

    % ------------------------- Check Reference Frames --------------------
    
    real_pos = [a.real_pos_xI(t); a.real_pos_yI(t); a.real_pos_zI(t)];
    real_vel = [a.real_vel_xI(t); a.real_vel_yI(t); a.real_vel_zI(t)];

    % Error in relative position (distance to center of operation region)
    p = pi-real_pos;
    p_err = p-a.pr_d;
    
    % ----------------- Inserting measured accel manually -----------------

    % Get REAL angle
    theta_real = a.theta(t);
    phi_real = a.phi(t);
    psi_real = a.psi(t);
    
    % Get REAL inertial acceleration (this has g in it)
    accel_S = a.measured_accel_3D(a, t, pm_ddot);

    % Get REAL angular velocity
    ang_rate_i = [a.theta_dot(t); a.phi_dot(t); a.psi_dot(t)];
    
    % Add sensor error to gyroscope measurements
    gyro_m = a.measured_gyro_3D(a, t, ang_rate_i);

    % Rotate inertial acceleration to sensor frame base on real angle/acceleration
    a_S = Rotate_I_S(accel_S, theta_real, phi_real, psi_real);

    theta_use = p_theta(1);
    phi_use = p_theta(2);
    psi_use = p_theta(3);

    % Rotate realistic sensor acceleration measurements back to inertial
    % frame (still has g), using our INTEGRATED angle (has integration error)
    accel_I = Rotate_S_I(a_S, theta_use, phi_use, psi_use);

    % Assign our acceleration and gyroscope measurements, with gravity
    % removed, to our vector for sensor error correction
    measured_state = [accel_I; gyro_m]';

    % Compensate for constant error values
    corrected_state = compensateError(measured_state, a.specs, t);

    % input: [velocity; theta; position; theta_err_accum]
    % output: [accel; theta_dot; vel; theta]
    int_state = [pm_dot; p_theta; pm; p_theta_err_accum];
    state_control = DriftCorrection3D(a, int_state, corrected_state);

    pm_ddot = state_control(1:3);
    p_theta_dot = state_control(4:6);

    pm_dot = state_control(7:9);
    p_theta_err = state_control(12:12);
    
    % ---------------------------------------------------------------------
    
    % Control gain proportions
    
    a.q1_ref_x = a.pr_d(1);
    a.q1_ref_y = a.pr_d(2);
    a.q1_ref_z = a.pr_d(3);

    d_x = a.q1_ref_x; % Center of input region
    d_y = a.q1_ref_y; % Center of input region
    d_z = a.q1_ref_z; % Center of input region

    C = [a.C(p(1), d_x); a.C(p(2), d_y); a.C(p(3), d_z)];
    B = [a.B(p(1), d_x); a.B(p(2), d_y); a.B(p(3), d_z)];

    % Calculate gains with gain mixing
    ka = a.ka.*C; % Acceleration [kg]
    kv = a.kv.*C; % Acceleration [kg]
    
    % Proportion of relative position control
    kp = a.kp_c.*C + a.kp_b.*B; % Proportional 
    kd = a.kd_c.*C + a.kd_b.*B; % Derivative   
    ki = a.ki_c.*C + a.ki_b.*B; % Integral   
    
    % Derivative of states
    state_dot = zeros(24,1);
    
    % Derivative of position
    state_dot(1:3) = p_dot;  % inertial velocity
    state_dot(10:12) = pm_dot; % measured inertial velocity
    
    % Control Law
    
    % Inertial stability control force
    c_i = a.initial_scale(t); % Initial scale of gains
    f_i = -(ka.*pm_ddot + kv.*pm_dot)*c_i;
    % Relative position control force
    f_pr = -(kp.*p_err + ki.*p_err_accum + kd.*(pm_dot-real_vel));
    
    % Platform EOM
    p_ddot = (f_i+f_pr) / a.m;
    
    % Derivative of velocity
    state_dot(4:6) = p_ddot; % Inertial acceleration
    state_dot(13:15) = pm_ddot; % Measured inertial acceleration
    
    % Derivative of measured inertial acceleration
    state_dot(16:18) = a.omega*(p_ddot - pm_ddot);
    % state_dot(16:18) = pm_ddot / a.dt;
    
    % Error in relative position
    state_dot(7:9) = p_err;

    % Derivative of angle
    state_dot(19:21) = p_theta_dot; % Angular Rate
    
    % Theta error
    state_dot(22:24) = p_theta_err;

end