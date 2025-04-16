function state_dot = rigidArmControl_FixedInt(t, a, prev_state)
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
    p = prev_state(1);
    p_dot = prev_state(2);
    pr_err_accum = prev_state(3);
    pm = prev_state(4);
    pm_dot = prev_state(5);
    pm_ddot = prev_state(6);
    p_theta = prev_state(7);
    p_theta_err_accum = prev_state(8);
    % pm_ddot_err_accum = prev_state(9);

    specs = a.specs;

    % ------------------------- Check Reference Frames --------------------

    % Error in relative position (distance to center of operation region)
    pr = pm-a.d(t);
    pr_err = pr-a.pr_d;
    
    % ----------------- Inserting measured accel manually -----------------
    
    % Get REAL angle
    theta_real = a.real_ang(t);
    phi_real = a.real_ang(t);
    psi_real = a.real_ang(t);
    
    % Get REAL inertial acceleration (this has g in it)
    a_I = [a.real_accel(t); a.real_accel(t); a.real_accel(t)];
    
    % Rotate inertial acceleration to sensor frame base on real angle/acceleration
    a_S = Rotate_I_S(a_I, theta_real, phi_real, psi_real);

    % Add sensor error to acceleration measurements
    accel_S = a.measured_accel_3D(t, a.o_d_n_a_c_3D, a.biasStabDistAccel, a.biasTempDistAccel, a.accel_drift_3D, a.noiseDistAccel, a_S);
    % accel_S = a_S;

    % Get REAL angular velocity
    theta_d = a.real_ang_rate(t);
    phi_d = a.real_ang_rate(t);
    psi_d = a.real_ang_rate(t);

    angle_d = [theta_d; phi_d; psi_d];

    % Add sensor error to gyroscope measurements
    gyro = a.measured_gyro_3D(t, a.o_d_n_g_c_3D, a.biasStabDistGyro, a.biasTempDistGyro, a.gyro_drift_3D, a.noiseDistGyro, angle_d);
    % gyro = angle_d;

    % Rotate realistic sensor acceleration measurements back to inertial
    % frame (still has g), using our INTEGRATED angle (has integration error)
    accel_I = Rotate_S_I(accel_S, p_theta, p_theta, p_theta);
    % accel_I = Rotate_S_I(accel_S, theta_real, phi_real, psi_real);
    % accel_I = a_I;

    % Simply remove g from inertial z vector
    % g is positive in simulation parameters. When measured it would be
    % negative though, hence why we add it here
    accel_I(3) = accel_I(3) + a.g;

    % Assign our acceleration and gyroscope measurements, with gravity
    % removed, to our vector for sensor error correction
    measured_state = [accel_I; gyro]';

    % Compensate for constant error values
    corrected_state = compensateError(measured_state, specs, t);

    % input: [velocity; theta; position; theta_err_accum]
    % output: [accel; theta_dot; vel; theta]
    int_state = [-p_dot; p_theta; pr_err; p_theta_err_accum];
    state_control = DriftCorrection1D(a, int_state, corrected_state);

    pm_ddot = state_control(1);
    p_thetadot = state_control(2);

    pm_dot = state_control(3);
    p_theta_err = state_control(4);
    
    % ---------------------------------------------------------------------
    
    % For testing gains without mixing proportions
    % ka = a.ka; % Acceleration [kg]
    % kv = a.kv; % Velocity     [kg/s]
    % ks = a.ks; % Position     [kg*s^-2]
    % kp = a.kp; % Proportional [kg*s^-2]
    % kd = a.kd; % Derivative   [kg/s]
    % ki = a.ki; % Integral     [kg*s^-3]
    
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
    state_dot = zeros(7,1);
    
    % Derivative of position
    state_dot(1) = p_dot;  % inertial velocity
    state_dot(4) = pm_dot; % measured inertial velocity
    
    % Control Law
    
    % Inertial stability control force
    % c_i = a.initial_scale(t); % Initial scale of gains
    c_i = 1;
    f_i = -(ka*pm_ddot + kv*pm_dot + ks*pm)*c_i;
    % f_i = 0.1
    % Relative position control force
    f_pr = -(kp*pr_err + ki*pr_err_accum + kd*(pm_dot-a.d_dot(t)));
    
    % Platform EOM
    p_ddot = (f_i+f_pr) / a.m;
    % p_ddot = 0;
    
    % scatter(a.fi1, t, p_ddot, sz, "blue")
    % hold on

    % if rem(t, 0.01)
    %     fprintf("", f_i)
    % end
    
    % Derivative of velocity
    state_dot(2) = p_ddot; % Inertial acceleration
    state_dot(5) = pm_ddot; % Measured inertial acceleration
    
    % Derivative of measured inertial acceleration
    state_dot(6) = a.omega*(p_ddot - pm_ddot);
    
    % Error in relative position
    state_dot(3) = pr_err;
    
    % Derivative of angle
    state_dot(7) = p_thetadot; % Angular Rate
    
    % Theta error
    state_dot(8) = p_theta_err;
    
    err_v=abs(p_dot-pm_dot);
    err_a=abs(p_ddot-pm_ddot);

    % t

end