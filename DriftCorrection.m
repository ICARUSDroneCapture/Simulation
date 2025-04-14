function s_dot = DriftCorrection(time_i, a, state, signals)

    % Current states
    vel_x = prev_state(1);
    vel_y = prev_state(2);
    vel_z = prev_state(3);

    angle_theta = prev_state(4);
    angle_phi = prev_state(5);
    angle_psi = prev_state(6);

    p_x = prev_state(7);
    p_y = prev_state(8);
    p_z = prev_state(9);

    theta_err_accum = prev_state(10);
    phi_err_accum = prev_state(11);
    psi_err_accum = prev_state(12);

    specs = a.specs;

    state = [vel_x; vel_y; vel_z; angle_theta; angle_phi; angle_psi];
    state_err_accum = [p_x; p_y; p_z; theta_err_accum; phi_err_accum; psi_err_accum];

    kw = a.kw;
    kt = a.kt;

    state_0 = 0;
    state_dot_0 = 0;
    
    % Get real inertial accelerations
    a_I = [a.real_accel_xI(time_i); a.real_accel_yI(time_i); a.real_accel_zI(time_i)];
    ang_rate_i = [a.theta_dot(time_i) ,a.phi_dot(time_i), a.psi_dot(time_i)];
    
    % Get real angles
    theta_real = a.theta(time_i);
    phi_real = a.phi(time_i);
    psi_real = a.psi(time_i);
    
    % Get real sensor frame accelerations
    a_S = Rotate_I_S(a_I, theta_real, phi_real, psi_real);
    
    accel_m = a.measured_accel_3D(a, time_i, a_S);


    gyro_m = a.measured_gyro_3D(a, time_i, ang_rate_i);

    theta_use = theta_real;
    phi_use = phi_real;
    psi_use = psi_real;

    if time_i > finishCalibrationTime
        theta_use = angle_theta;
        phi_use = angle_phi;
        psi_use = angle_psi;
    end

    accel_I = Rotate_S_I(accel_m, theta_use, phi_use, psi_use);

    accel_state = [accel_I(1) accel_I(2) accel_I(3)+a.g];

    measuredState = [accel_state, gyro_m];
    corrected_state = compensateError(measuredState, specs, time_i);
    
    state_dot_m = corrected_state';

    state_control = state - state_0;

    state_dot_comp = kt .* state_err_accum + kw .* state_control;

    state_dot = state_dot_m + state_dot_0 - state_dot_comp;

    s_dot = zeros(12,1);

    s_dot(1:6) = state_dot;
    s_dot(7:12) = state_control;

end