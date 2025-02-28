function s_dot = rigidArmControl_ODEFunc(t, s, a)
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
    p = s(1);
    p_dot = s(2);
    pr_err_accum = s(3);
    pm = s(4);
    pm_dot = s(5);
    pm_ddot = s(6);
    p_theta = s(7);
    
    specs = a.specs;
    
    p_thetadot = 1/(pm_dot^2 + 1);
    
    
    
    % Inserting measured accel manually ---------------------------------------
    sz = 3;
    if t > specs.accel_resolution
    
        measured_a_x = a.measured_accel_horz(t, a.o_d_n_a_c_h, a.biasStabDistAccel, a.biasTempDistAccel, a.accel_drift_horz, a.noiseDistAccel, a.real_accel, a.theta_err);
        measured_a_y = measured_a_x;
        measured_a_z = a.measured_accel_vert(t, a.o_d_n_a_c_v, a.biasStabDistAccel, a.biasTempDistAccel, a.accel_drift_vert, a.noiseDistAccel, a.real_accel, a.theta_err);
    
        measured_g_x = a.measured_gyro(t, a.o_d_n_g_c, a.biasStabDistGyro, a.biasTempDistGyro, a.gyro_drift, a.noiseDistGyro, a.real_ang_rate);
        measured_g_y = measured_g_x;
        measured_g_z = measured_g_x;
    
        measuredState = [measured_a_x measured_a_y measured_a_z measured_g_x measured_g_y measured_g_z];
        corrected_a = compensateError(measuredState, specs, t);
    
        pm_ddot_error = pm_ddot - corrected_a(3);
        p_thetadot_error = p_thetadot - corrected_a(4);
    
        pm_ddot = pm_ddot + (pm_ddot_error/1);
        p_thetadot = p_thetadot + (p_thetadot_error/1);
    
        % pm_ddot = corrected_a(3);
        % p_thetadot = corrected_a(4);
    
        % scatter(a.fi1, t, a_error, "red")
        % hold on
        % 
        % scatter(a.fi1, t, pm_ddot, "blue")
        % hold on
        % 
        % scatter(a.fi1, t, corrected_a(3), "green")
        % hold on
    
    end
    
    
    % scatter(a.fi1, t, corrected_a(1), sz, "filled", color='blue')
    % hold on
    % scatter(a.fi1, t, corrected_a(3), sz, color='blue')
    % hold on
    
    % 
    % error = abs(corrected_a(3) - p_ddot);
    % scatter(a.fi1, t, error, sz, "red")
    % hold on
    % scatter(a.fi3, t, corrected_a(4), sz, "filled", color='blue')
    % hold on
    
    % ----------------------------------------------------------------------------------------------
    
    % Error in relative position (distance to center of operation region)
    pr = pm-a.d(t);
    pr_err = pr-a.pr_d;
    
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
    s_dot = zeros(7,1);
    
    % Derivative of position
    s_dot(1) = p_dot;  % inertial velocity
    s_dot(4) = pm_dot; % measured inertial velocity
    
    % Control Law
    
    % Inertial stability control force
    c_i = a.initial_scale(t); % Initial scale of gains
    f_i = -(ka*pm_ddot + kv*pm_dot + ks*pm)*c_i;
    % Relative position control force
    f_pr = -(kp*pr_err + ki*pr_err_accum + kd*(pm_dot-a.d_dot(t)));
    
    % Platform EOM
    p_ddot = (f_i+f_pr) / a.m;
    
    % scatter(a.fi1, t, p_ddot, sz, "blue")
    % hold on
    
    % Derivative of velocity
    s_dot(2) = p_ddot; % Inertial acceleration
    s_dot(5) = pm_ddot; % Measured inertial acceleration
    
    % Derivative of measured inertial acceleration
    s_dot(6) = a.omega*(p_ddot - pm_ddot);
    
    % Error in relative position
    s_dot(3) = pr_err;
    
    % Derivative of angle
    s_dot(7) = p_thetadot; % Angular Rate
    
    err_v=abs(p_dot-pm_dot);
    err_a=abs(p_ddot-pm_ddot);
    
    % t

end