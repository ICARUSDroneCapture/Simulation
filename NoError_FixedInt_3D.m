function state_dot = NoError_FixedInt_3D(t, a, prev_state)

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

    pi = [pi_x; pi_y; pi_z];
    p_dot = [p_dot_x; p_dot_y; p_dot_z];
    p_err_accum = [p_err_accum_x; p_err_accum_y; p_err_accum_z];

    pm = [pm_x; pm_y; pm_z];
    pm_dot = [pm_dot_x; pm_dot_y; pm_dot_z];
    pm_ddot = [pm_ddot_x; pm_ddot_y; pm_ddot_z];

    % ------------------------- Check Reference Frames --------------------
    
    real_pos = [a.real_pos_xI(t); a.real_pos_yI(t); a.real_pos_zI(t)];
    real_vel = [a.real_vel_xI(t); a.real_vel_yI(t); a.real_vel_zI(t)];

    % Error in relative position (distance to center of operation region)
    p = pi-real_pos;
    p_err = p-a.pr_d;
    
    % ---------------------------------------------------------------------
    
    % For testing gains without mixing proportions
    % ka = a.ka; % Acceleration [kg]
    % kv = a.kv; % Velocity     [kg/s]
    % ks = a.ks; % Position     [kg*s^-2]
    % kp = a.kp; % Proportional [kg*s^-2]
    % kd = a.kd; % Derivative   [kg/s]
    % ki = a.ki; % Integral     [kg*s^-3]

    % Control gain proportions
    
    a.q1_ref_x = a.pr_d(1);
    a.q1_ref_y = a.pr_d(2);
    a.q1_ref_z = a.hdeck + a.pr_d(3);

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
    state_dot = zeros(18,1);
    
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
    
    % Error in relative position
    state_dot(7:9) = p_err;

end