function s_dot = DriftCorrection1D(a, state, signals)

    % Current states
    vel = state(1);
    theta = state(2);
    p = state(3);
    theta_err_accum = state(4);

    curr_state = [vel; theta];
    state_err_accum = [p; theta_err_accum];

    kw = a.kw(3);
    kt = a.kt(5);

    state_0 = [0; 0];
    state_dot_0 = [0; 0];

    accel_dot_m_v = signals(3);
    gyro_dot_m = signals(5);

    state_dot_m = [accel_dot_m_v; gyro_dot_m];

    state_control = curr_state - state_0;

    state_dot_comp = kt * state_err_accum + kw * state_control;

    state_dot = state_dot_m + state_dot_0 - state_dot_comp;

    s_dot = zeros(length(state),1);

    s_dot(1:2) = state_dot;
    s_dot(3:4) = state_control;

end