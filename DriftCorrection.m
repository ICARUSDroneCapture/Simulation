function s_dot = DriftCorrection(time_i, a, prev_state)

    % Current states
    vel_x = prev_state(1);
    vel_y = prev_state(2);
    vel_z = prev_state(3);
    theta = prev_state(4);
    psi = prev_state(5);
    phi = prev_state(6);
    vel_err_accum_x = prev_state(7);
    vel_err_accum_y = prev_state(8);
    vel_err_accum_z = prev_state(9);
    theta_err_accum = prev_state(10);
    psi_err_accum = prev_state(11);
    phi_err_accum = prev_state(12);

    curr_state = [vel_x; vel_y; vel_z; theta; psi; phi];
    state_err_accum = [vel_err_accum_x; vel_err_accum_y; vel_err_accum_z; theta_err_accum; psi_err_accum; phi_err_accum];

    specs = a.specs;

    kw = a.kw; % 
    kt = a.kt; % 

    state_0 = 0;
    state_dot_0 = 0;

    accel_i = a.real_accel(time_i);
    ang_rate_i = a.real_ang_rate(time_i);
    
    accel_dot_m_v = a.measured_accel_vert(time_i, a.o_d_n_a_c_v, a.biasStabDistAccel, a.biasTempDistAccel, a.accel_drift_vert, a.noiseDistAccel, accel_i, a.theta_err);
    accel_dot_m_h = a.measured_accel_horz(time_i, a.o_d_n_a_c_h, a.biasStabDistAccel, a.biasTempDistAccel, a.accel_drift_horz, a.noiseDistAccel, accel_i, a.theta_err);
    
    gyro_dot_m = a.measured_gyro(time_i, a.o_d_n_g_c, a.biasStabDistGyro, a.biasTempDistGyro, a.gyro_drift, a.noiseDistGyro, ang_rate_i);

    state_dot_m = [accel_dot_m_h; accel_dot_m_h; accel_dot_m_v; gyro_dot_m; gyro_dot_m; gyro_dot_m];

    state_control = curr_state - state_0;

    state_dot_comp = kt * state_err_accum + kw * state_control;

    state_dot = state_dot_m + state_dot_0 - state_dot_comp;

    s_dot = zeros(length(prev_state),1);

    s_dot(1:6) = state_dot;
    s_dot(7:12) = state_control;

end