function state_dot = rigidArmControl_1D(t, a, prev_state, finishCalibrationTime)

    % rigidArmControl is the EOM for the 1 DOF model of the inertially
    % stabilized platform. It uses inertial acceleration control when the
    % platform is close to the center of the operation region, and uses PID
    % control on the relative position as the platform goes closer to the
    % operational bounderies
    %
    % Inputs:   t    = current time
    %           a    = structure containing environmental constants and gain
    %                  values
    %           s    = vector of states
    %                = [pi; p_dot; p_err_accum; pm; pm_dot; pm_ddot] 
    %                  - pi: platform inertial position 
    %                  - p_dot: platform inertial velocity
    %                  - p_err_accum: platform relative position correction accumulation
    %                  - pm: measured platform relative position
    %                  - pm_dot: measured platform inertial velocity 
    %                  - pm_ddot: measured platform intertial acceleration
    % Outputs:  sdot = time derivative of input state vector
    %                = [p_dot; p_ddot; pr_err] where pdot is the inertial 
    %                  velocity of the platform, p_ddot is the inertial 
    %                  acceleration of the platform,and pr_err is the error in 
    %                  the relative position of the platform
    
    % Current states
    pi = prev_state(1);
    p_dot = prev_state(2);
    p_err_accum = prev_state(3);
    pm = prev_state(4);
    pm_dot = prev_state(5);
    pm_ddot = prev_state(6);
    p_theta = prev_state(7);
    p_theta_err_accum = prev_state(8);
    
    % Error in relative position (distance to center of operation region)
    p = pi-a.real_pos_zI(t);
    p_err = p-a.pr_d(3);
    
    % ----------------- Inserting measured accel manually -----------------

    % Get REAL angle
    theta_real = a.theta(t);
    phi_real = a.phi(t);
    psi_real = a.psi(t);

    % Get REAL inertial acceleration (this has g in it)
    accel_i = a.measured_accel_3D(a, t, pm_ddot);
    accel_S = [0; 0; accel_i];

    % Get REAL angular velocity
    ang_rate_i = [a.theta_dot(t); a.phi_dot(t); a.psi_dot(t)];
    
    % Add sensor error to gyroscope measurements
    gyro_m = a.measured_gyro_3D(a, t, ang_rate_i);

    % Rotate inertial acceleration to sensor frame base on real angle/acceleration
    a_S = Rotate_I_S(accel_S, theta_real, phi_real, psi_real);

    theta_use = theta_real;
    phi_use = p_theta;
    psi_use = psi_real;

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
    state_control = DriftCorrection1D(a, int_state, corrected_state);

    pm_ddot = state_control(1);
    p_theta_dot = state_control(2);

    pm_dot = state_control(3);
    p_theta_err = state_control(4);

    % ---------------------------------------------------------------------
    

    % Control gain proportions
    a.q1_ref_z = a.pr_d(3);

    d_z = a.q1_ref_z; % Center of input region

    C = a.C(p, d_z);
    B = a.B(p, d_z);

    % Calculate gains with gain mixing
    ka = a.ka(3).*C; % Acceleration [kg]
    kv = a.kv(3).*C; % Acceleration [kg]
    
    % Proportion of relative position control
    kp = a.kp_c(3).*C + a.kp_b(3).*B; % Proportional 
    kd = a.kd_c(3).*C + a.kd_b(3).*B; % Derivative   
    ki = a.ki_c(3).*C + a.ki_b(3).*B; % Integral   

    % Derivative of states
    state_dot = zeros(6,1);
    
    % Derivative of position
    state_dot(1) = p_dot;  % inertial velocity
    state_dot(4) = pm_dot; % measured inertial velocity
    
    % Control Law
    
    % Inertial stability control force
    c_i = a.initial_scale(t); % Initial scale of gains
    f_i = -(ka*pm_ddot + kv*pm_dot)*c_i;
    % Relative position control force
    f_pr = -(kp*p_err + ki*p_err_accum + kd*(pm_dot-a.real_vel_zI(t)));
    
    % Platform EOM
    p_ddot = (f_i+f_pr) / a.m;
    
    % Derivative of velocity
    state_dot(2) = p_ddot; % Inertial acceleration
    state_dot(5) = pm_ddot; % Measured inertial acceleration
    
    % Derivative of measured inertial acceleration
    state_dot(6) = a.omega*(p_ddot - pm_ddot);
    
    % Error in relative position
    state_dot(3) = p_err;

    state_dot(7) = p_theta_dot;
    state_dot(8) = p_theta_err;

end